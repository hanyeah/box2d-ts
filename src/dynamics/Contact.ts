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
  /// Friction mixing law. The idea is to allow either fixture to drive the friction to zero.
/// For example, anything slides on ice.
  export function MixFriction(friction1: number, friction2: number): number {
    return Sqrt(friction1 * friction2);
  }

/// Restitution mixing law. The idea is allow for anything to bounce off an inelastic surface.
/// For example, a superball bounces on anything.
  export function MixRestitution(restitution1: number, restitution2: number): number {
    return restitution1 > restitution2 ? restitution1 : restitution2;
  }

/// Restitution mixing law. This picks the lowest value.
  export function MixRestitutionThreshold(threshold1: number, threshold2: number): number {
    return threshold1 < threshold2 ? threshold1 : threshold2;
  }

  export class ContactEdge {
    private _other: Body | null = null; ///< provides quick access to the other body attached.
    public get other(): Body {
      if (this._other === null) { throw new Error(); }
      return this._other;
    }
    public set other(value: Body) {
      if (this._other !== null) { throw new Error(); }
      this._other = value;
    }
    public readonly contact: Contact; ///< the contact
    public prev: ContactEdge | null = null; ///< the previous contact edge in the body's contact list
    public next: ContactEdge | null = null; ///< the next contact edge in the body's contact list
    constructor(contact: Contact) {
      this.contact = contact;
    }
    public Reset(): void {
      this._other = null;
      this.prev = null;
      this.next = null;
    }
  }

  export abstract class Contact<A extends Shape = Shape, B extends Shape = Shape> {
    public islandFlag: boolean = false; /// Used when crawling contact graph when forming islands.
    public touchingFlag: boolean = false; /// Set when the shapes are touching.
    public enabledFlag: boolean = false; /// This contact can be disabled (by user)
    public filterFlag: boolean = false; /// This contact needs filtering because a fixture filter was changed.
    public bulletHitFlag: boolean = false; /// This bullet contact had a TOI event
    public toiFlag: boolean = false; /// This contact has a valid TOI in toi

    public prev: Contact | null = null;
    public next: Contact | null = null;

    public readonly nodeA: ContactEdge = new ContactEdge(this);
    public readonly nodeB: ContactEdge = new ContactEdge(this);

    public fixtureA!: Fixture;
    public fixtureB!: Fixture;

    public indexA: number = 0;
    public indexB: number = 0;

    public manifold: Manifold = new Manifold(); // TODO: readonly

    public toiCount: number = 0;
    public toi: number = 0;

    public friction: number = 0;
    public restitution: number = 0;
    public restitutionThreshold: number = 0;

    public tangentSpeed: number = 0;

    public oldManifold: Manifold = new Manifold(); // TODO: readonly

    public GetManifold() {
      return this.manifold;
    }

    public GetWorldManifold(worldManifold: WorldManifold): void {
      const bodyA: Body = this.fixtureA.GetBody();
      const bodyB: Body = this.fixtureB.GetBody();
      const shapeA: A = this.GetShapeA();
      const shapeB: B = this.GetShapeB();
      worldManifold.Initialize(this.manifold, bodyA.GetTransform(), shapeA.radius, bodyB.GetTransform(), shapeB.radius);
    }

    public IsTouching(): boolean {
      return this.touchingFlag;
    }

    public SetEnabled(flag: boolean): void {
      this.enabledFlag = flag;
    }

    public IsEnabled(): boolean {
      return this.enabledFlag;
    }

    public GetNext(): Contact | null {
      return this.next;
    }

    public GetFixtureA(): Fixture {
      return this.fixtureA;
    }

    public GetChildIndexA(): number {
      return this.indexA;
    }

    public GetShapeA(): A {
      return this.fixtureA.GetShape() as A;
    }

    public GetFixtureB(): Fixture {
      return this.fixtureB;
    }

    public GetChildIndexB(): number {
      return this.indexB;
    }

    public GetShapeB(): B {
      return this.fixtureB.GetShape() as B;
    }

    public abstract Evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;

    public FlagForFiltering(): void {
      this.filterFlag = true;
    }

    public SetFriction(friction: number): void {
      this.friction = friction;
    }

    public GetFriction(): number {
      return this.friction;
    }

    public ResetFriction(): void {
      this.friction = MixFriction(this.fixtureA.friction, this.fixtureB.friction);
    }

    public SetRestitution(restitution: number): void {
      this.restitution = restitution;
    }

    public GetRestitution(): number {
      return this.restitution;
    }

    public ResetRestitution(): void {
      this.restitution = MixRestitution(this.fixtureA.restitution, this.fixtureB.restitution);
    }

    /// Override the default restitution velocity threshold mixture. You can call this in ContactListener::PreSolve.
    /// The value persists until you set or reset.
    public SetRestitutionThreshold(threshold: number): void {
      this.restitutionThreshold = threshold;
    }

    /// Get the restitution threshold.
    public GetRestitutionThreshold(): number {
      return this.restitutionThreshold;
    }

    /// Reset the restitution threshold to the default value.
    public ResetRestitutionThreshold(): void {
      this.restitutionThreshold = MixRestitutionThreshold(this.fixtureA.restitutionThreshold, this.fixtureB.restitutionThreshold);
    }

    public SetTangentSpeed(speed: number): void {
      this.tangentSpeed = speed;
    }

    public GetTangentSpeed(): number {
      return this.tangentSpeed;
    }

    public Reset(fixtureA: Fixture, indexA: number, fixtureB: Fixture, indexB: number): void {
      this.islandFlag = false;
      this.touchingFlag = false;
      this.enabledFlag = true;
      this.filterFlag = false;
      this.bulletHitFlag = false;
      this.toiFlag = false;

      this.fixtureA = fixtureA;
      this.fixtureB = fixtureB;

      this.indexA = indexA;
      this.indexB = indexB;

      this.manifold.pointCount = 0;

      this.prev = null;
      this.next = null;

      this.nodeA.Reset();
      this.nodeB.Reset();

      this.toiCount = 0;

      this.friction = MixFriction(this.fixtureA.friction, this.fixtureB.friction);
      this.restitution = MixRestitution(this.fixtureA.restitution, this.fixtureB.restitution);
      this.restitutionThreshold = MixRestitutionThreshold(this.fixtureA.restitutionThreshold, this.fixtureB.restitutionThreshold);
    }

    public Update(listener: ContactListener): void {
      const tManifold: Manifold = this.oldManifold;
      this.oldManifold = this.manifold;
      this.manifold = tManifold;

      // Re-enable this contact.
      this.enabledFlag = true;

      let touching: boolean = false;
      const wasTouching: boolean = this.touchingFlag;

      const sensorA: boolean = this.fixtureA.IsSensor();
      const sensorB: boolean = this.fixtureB.IsSensor();
      const sensor: boolean = sensorA || sensorB;

      const bodyA: Body = this.fixtureA.GetBody();
      const bodyB: Body = this.fixtureB.GetBody();
      const xfA: Transform = bodyA.GetTransform();
      const xfB: Transform = bodyB.GetTransform();

      // Is this contact a sensor?
      if (sensor) {
        const shapeA: A = this.GetShapeA();
        const shapeB: B = this.GetShapeB();
        touching = TestOverlapShape(shapeA, this.indexA, shapeB, this.indexB, xfA, xfB);

        // Sensors don't generate manifolds.
        this.manifold.pointCount = 0;
      } else {
        this.Evaluate(this.manifold, xfA, xfB);
        touching = this.manifold.pointCount > 0;

        // Match old contact ids to new contact ids and copy the
        // stored impulses to warm start the solver.
        for (let i: number = 0; i < this.manifold.pointCount; ++i) {
          const mp2: ManifoldPoint = this.manifold.points[i];
          mp2.normalImpulse = 0;
          mp2.tangentImpulse = 0;
          const id2: ContactID = mp2.id;

          for (let j: number = 0; j < this.oldManifold.pointCount; ++j) {
            const mp1: ManifoldPoint = this.oldManifold.points[j];

            if (mp1.id.key === id2.key) {
              mp2.normalImpulse = mp1.normalImpulse;
              mp2.tangentImpulse = mp1.tangentImpulse;
              break;
            }
          }
        }

        if (touching !== wasTouching) {
          bodyA.SetAwake(true);
          bodyB.SetAwake(true);
        }
      }

      this.touchingFlag = touching;

      if (!wasTouching && touching && listener) {
        listener.BeginContact(this);
      }

      if (wasTouching && !touching && listener) {
        listener.EndContact(this);
      }

      if (!sensor && touching && listener) {
        listener.PreSolve(this, this.oldManifold);
      }
    }

    private static ComputeTOI_s_input = new TOIInput();
    private static ComputeTOI_s_output = new TOIOutput();
    public ComputeTOI(sweepA: Sweep, sweepB: Sweep): number {
      const input: TOIInput = Contact.ComputeTOI_s_input;
      input.proxyA.SetShape(this.GetShapeA(), this.indexA);
      input.proxyB.SetShape(this.GetShapeB(), this.indexB);
      input.sweepA.Copy(sweepA);
      input.sweepB.Copy(sweepB);
      input.tMax = linearSlop;

      const output: TOIOutput = Contact.ComputeTOI_s_output;

      TimeOfImpact(output, input);

      return output.t;
    }
  }

}
