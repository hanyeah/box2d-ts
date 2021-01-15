/*
* Copyright (c) 2006-2009 Erin Catto http://www.box2d.org
*
* This software is provided 'as-is', without any express or implied
* warranty.  In no event will the authors be held liable for any damages
* arising from the use of this software.
* Permission is granted to anyone to use this software for any purpose,
* including commercial applications, and to alter it and redistribute it
* freely, subject to the following restrictions:
* 1. The origin of this software must not be misrepresented; you must not
* claim that you wrote the original software. If you use this software
* in a product, an acknowledgment in the product documentation would be
* appreciated but is not required.
* 2. Altered source versions must be plainly marked as such, and must not be
* misrepresented as being the original software.
* 3. This notice may not be removed or altered from any source distribution.
*/
namespace b2 {
  // Delegate of World.
  export class ContactManager {
    public readonly broadPhase: BroadPhase<FixtureProxy> = new BroadPhase<FixtureProxy>();
    public contactList: Contact = null;
    public contactCount: number = 0;
    public contactFilter: ContactFilter = ContactFilter.defaultFilter;
    public contactListener: ContactListener = ContactListener.defaultListener;

    public readonly contactFactory: ContactFactory = new ContactFactory();

    // Broad-phase callback.
    public addPair(proxyA: FixtureProxy, proxyB: FixtureProxy): void {
      // DEBUG: Assert(proxyA instanceof FixtureProxy);
      // DEBUG: Assert(proxyB instanceof FixtureProxy);

      let fixtureA: Fixture = proxyA.fixture;
      let fixtureB: Fixture = proxyB.fixture;

      let indexA: number = proxyA.childIndex;
      let indexB: number = proxyB.childIndex;

      let bodyA: Body = fixtureA.getBody();
      let bodyB: Body = fixtureB.getBody();

      // Are the fixtures on the same body?
      if (bodyA === bodyB) {
        return;
      }

      // TODO_ERIN use a hash table to remove a potential bottleneck when both
      // bodies have a lot of contacts.
      // Does a contact already exist?
      let edge: ContactEdge = bodyB.getContactList();
      while (edge) {
        if (edge.other === bodyA) {
          const fA: Fixture = edge.contact.getFixtureA();
          const fB: Fixture = edge.contact.getFixtureB();
          const iA: number = edge.contact.getChildIndexA();
          const iB: number = edge.contact.getChildIndexB();

          if (fA === fixtureA && fB === fixtureB && iA === indexA && iB === indexB) {
            // A contact already exists.
            return;
          }

          if (fA === fixtureB && fB === fixtureA && iA === indexB && iB === indexA) {
            // A contact already exists.
            return;
          }
        }

        edge = edge.next;
      }

      // Check user filtering.
      if (this.contactFilter && !this.contactFilter.ShouldCollide(fixtureA, fixtureB)) {
        return;
      }

      // Call the factory.
      const c: Contact = this.contactFactory.create(fixtureA, indexA, fixtureB, indexB);
      if (c === null) {
        return;
      }

      // Contact creation may swap fixtures.
      fixtureA = c.getFixtureA();
      fixtureB = c.getFixtureB();
      indexA = c.getChildIndexA();
      indexB = c.getChildIndexB();
      bodyA = fixtureA.body;
      bodyB = fixtureB.body;

      // Insert into the world.
      c.prev = null;
      c.next = this.contactList;
      if (this.contactList !== null) {
        this.contactList.prev = c;
      }
      this.contactList = c;

      // Connect to island graph.

      // Connect to body A
      c.nodeA.other = bodyB;

      c.nodeA.prev = null;
      c.nodeA.next = bodyA.contactList;
      if (bodyA.contactList !== null) {
        bodyA.contactList.prev = c.nodeA;
      }
      bodyA.contactList = c.nodeA;

      // Connect to body B
      c.nodeB.other = bodyA;

      c.nodeB.prev = null;
      c.nodeB.next = bodyB.contactList;
      if (bodyB.contactList !== null) {
        bodyB.contactList.prev = c.nodeB;
      }
      bodyB.contactList = c.nodeB;

      ++this.contactCount;
    }

    public findNewContacts(): void {
      this.broadPhase.updatePairs((proxyA: FixtureProxy, proxyB: FixtureProxy): void => {
        this.addPair(proxyA, proxyB);
      });
    }

    public destroy(c: Contact): void {
      const fixtureA: Fixture = c.getFixtureA();
      const fixtureB: Fixture = c.getFixtureB();
      const bodyA: Body = fixtureA.getBody();
      const bodyB: Body = fixtureB.getBody();

      if (this.contactListener && c.isTouching()) {
        this.contactListener.endContact(c);
      }

      // Remove from the world.
      if (c.prev) {
        c.prev.next = c.next;
      }

      if (c.next) {
        c.next.prev = c.prev;
      }

      if (c === this.contactList) {
        this.contactList = c.next;
      }

      // Remove from body 1
      if (c.nodeA.prev) {
        c.nodeA.prev.next = c.nodeA.next;
      }

      if (c.nodeA.next) {
        c.nodeA.next.prev = c.nodeA.prev;
      }

      if (c.nodeA === bodyA.contactList) {
        bodyA.contactList = c.nodeA.next;
      }

      // Remove from body 2
      if (c.nodeB.prev) {
        c.nodeB.prev.next = c.nodeB.next;
      }

      if (c.nodeB.next) {
        c.nodeB.next.prev = c.nodeB.prev;
      }

      if (c.nodeB === bodyB.contactList) {
        bodyB.contactList = c.nodeB.next;
      }

      // moved this from ContactFactory:Destroy
      if (c.manifold.pointCount > 0 &&
        !fixtureA.isSensor &&
        !fixtureB.isSensor) {
        fixtureA.getBody().setAwake(true);
        fixtureB.getBody().setAwake(true);
      }

      // Call the factory.
      this.contactFactory.destroy(c);
      --this.contactCount;
    }

    // This is the top level collision call for the time step. Here
    // all the narrow phase collision is processed for the world
    // contact list.
    public collide(): void {
      // Update awake contacts.
      let c: Contact = this.contactList;
      while (c) {
        const fixtureA: Fixture = c.getFixtureA();
        const fixtureB: Fixture = c.getFixtureB();
        const indexA: number = c.getChildIndexA();
        const indexB: number = c.getChildIndexB();
        const bodyA: Body = fixtureA.getBody();
        const bodyB: Body = fixtureB.getBody();

        // Is this contact flagged for filtering?
        if (c.filterFlag) {
          // Check user filtering.
          if (this.contactFilter && !this.contactFilter.ShouldCollide(fixtureA, fixtureB)) {
            const cNuke: Contact = c;
            c = cNuke.next;
            this.destroy(cNuke);
            continue;
          }

          // Clear the filtering flag.
          c.filterFlag = false;
        }

        const activeA: boolean = bodyA.isAwake() && bodyA.type !== BodyType.StaticBody;
        const activeB: boolean = bodyB.isAwake() && bodyB.type !== BodyType.StaticBody;

        // At least one body must be awake and it must be dynamic or kinematic.
        if (!activeA && !activeB) {
          c = c.next;
          continue;
        }

        const treeNodeA: TreeNode<FixtureProxy> = fixtureA.proxies[indexA].treeNode;
        const treeNodeB: TreeNode<FixtureProxy> = fixtureB.proxies[indexB].treeNode;
        const overlap: boolean = testOverlapAABB(treeNodeA.aabb, treeNodeB.aabb);

        // Here we destroy contacts that cease to overlap in the broad-phase.
        if (!overlap) {
          const cNuke: Contact = c;
          c = cNuke.next;
          this.destroy(cNuke);
          continue;
        }

        // The contact persists.
        c.update(this.contactListener);
        c = c.next;
      }
    }
  }

}
