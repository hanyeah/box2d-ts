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
  /// This holds contact filtering data.
  export interface IFilter {
    /// The collision category bits. Normally you would just set one bit.
    categoryBits: number;

    /// The collision mask bits. This states the categories that this
    /// shape would accept for collision.
    maskBits: number;

    /// Collision groups allow a certain group of objects to never collide (negative)
    /// or always collide (positive). Zero means no collision group. Non-zero group
    /// filtering always wins against the mask bits.
    groupIndex?: number;
  }

/// This holds contact filtering data.
  export class Filter implements IFilter {
    public static readonly DEFAULT: Filter = new Filter();

    /// The collision category bits. Normally you would just set one bit.
    public categoryBits: number = 0x0001;

    /// The collision mask bits. This states the categories that this
    /// shape would accept for collision.
    public maskBits: number = 0xFFFF;

    /// Collision groups allow a certain group of objects to never collide (negative)
    /// or always collide (positive). Zero means no collision group. Non-zero group
    /// filtering always wins against the mask bits.
    public groupIndex: number = 0;

    public clone(): Filter {
      return new Filter().copy(this);
    }

    public copy(other: IFilter): this {
      // DEBUG: Assert(this !== other);
      this.categoryBits = other.categoryBits;
      this.maskBits = other.maskBits;
      this.groupIndex = other.groupIndex || 0;
      return this;
    }
  }

/// A fixture definition is used to create a fixture. This class defines an
/// abstract fixture definition. You can reuse fixture definitions safely.
  export interface IFixtureDef {
    /// The shape, this must be set. The shape will be cloned, so you
    /// can create the shape on the stack.
    shape: Shape;

    /// Use this to store application specific fixture data.
    userData?: any;

    /// The friction coefficient, usually in the range [0,1].
    friction?: number;

    /// The restitution (elasticity) usually in the range [0,1].
    restitution?: number;

    /// Restitution velocity threshold, usually in m/s. Collisions above this
    /// speed have restitution applied (will bounce).
    restitutionThreshold?: number;

    /// The density, usually in kg/m^2.
    density?: number;

    /// A sensor shape collects contact information but never generates a collision
    /// response.
    isSensor?: boolean;

    /// Contact filtering data.
    filter?: IFilter;
  }

/// A fixture definition is used to create a fixture. This class defines an
/// abstract fixture definition. You can reuse fixture definitions safely.
  export class FixtureDef implements IFixtureDef {
    /// The shape, this must be set. The shape will be cloned, so you
    /// can create the shape on the stack.
    public shape!: Shape;

    /// Use this to store application specific fixture data.
    public userData: any = null;

    /// The friction coefficient, usually in the range [0,1].
    public friction: number = 0.2;

    /// The restitution (elasticity) usually in the range [0,1].
    public restitution: number = 0;

    /// Restitution velocity threshold, usually in m/s. Collisions above this
    /// speed have restitution applied (will bounce).
    public restitutionThreshold: number = 1.0 * lengthUnitsPerMeter;

    /// The density, usually in kg/m^2.
    public density: number = 0;

    /// A sensor shape collects contact information but never generates a collision
    /// response.
    public isSensor: boolean = false;

    /// Contact filtering data.
    public readonly filter: Filter = new Filter();
  }

/// This proxy is used internally to connect fixtures to the broad-phase.
  export class FixtureProxy {
    public readonly aabb: AABB = new AABB();
    public readonly fixture: Fixture;
    public readonly childIndex: number = 0;
    public treeNode: TreeNode<FixtureProxy>;
    constructor(fixture: Fixture, childIndex: number) {
      this.fixture = fixture;
      this.childIndex = childIndex;
      this.fixture.shape.computeAABB(this.aabb, this.fixture.body.getTransform(), childIndex);
      this.treeNode = this.fixture.body.world.contactManager.broadPhase.createProxy(this.aabb, this);
    }
    public reset(): void {
      this.fixture.body.world.contactManager.broadPhase.destroyProxy(this.treeNode);
    }
    public touch(): void {
      this.fixture.body.world.contactManager.broadPhase.touchProxy(this.treeNode);
    }
    private static synchronize_s_aabb1 = new AABB();
    private static synchronize_s_aab = new AABB();
    private static synchronize_s_displacement = new Vec2();
    public synchronize(transform1: Transform, transform2: Transform): void {
      if (transform1 === transform2) {
        this.fixture.shape.computeAABB(this.aabb, transform1, this.childIndex);
        this.fixture.body.world.contactManager.broadPhase.moveProxy(this.treeNode, this.aabb, Vec2.ZERO);
      } else {
        // Compute an AABB that covers the swept shape (may miss some rotation effect).
        const aabb1: AABB = FixtureProxy.synchronize_s_aabb1;
        const aab: AABB = FixtureProxy.synchronize_s_aab;
        this.fixture.shape.computeAABB(aabb1, transform1, this.childIndex);
        this.fixture.shape.computeAABB(aab, transform2, this.childIndex);
        this.aabb.combine2(aabb1, aab);
        const displacement: Vec2 = FixtureProxy.synchronize_s_displacement;
        displacement.copy(aab.getCenter()).selfSub(aabb1.getCenter());
        this.fixture.body.world.contactManager.broadPhase.moveProxy(this.treeNode, this.aabb, displacement);
      }
    }
  }

/// A fixture is used to attach a shape to a body for collision detection. A fixture
/// inherits its transform from its parent. Fixtures hold additional non-geometric data
/// such as friction, collision filters, etc.
/// Fixtures are created via Body::CreateFixture.
/// @warning you cannot reuse fixtures.
  export class Fixture {
    public density: number = 0;

    public next: Fixture = null;
    public readonly body: Body;

    public readonly shape: Shape;

    public friction: number = 0;
    public restitution: number = 0;
    public restitutionThreshold: number = 1.0 * lengthUnitsPerMeter;

    public readonly proxies: FixtureProxy[] = [];
    public get proxyCount(): number { return this.proxies.length; }

    public readonly filter: Filter = new Filter();

    public isSensor: boolean = false;

    public userData: any = null;

    constructor(body: Body, def: IFixtureDef) {
      this.body = body;
      this.shape = def.shape.clone();
      this.userData = maybe(def.userData, null);
      this.friction = maybe(def.friction, 0.2);
      this.restitution = maybe(def.restitution, 0);
      this.restitutionThreshold = maybe(def.restitutionThreshold, 0);
      this.filter.copy(maybe(def.filter, Filter.DEFAULT));
      this.isSensor = maybe(def.isSensor, false);
      this.density = maybe(def.density, 0);
    }

    public reset(): void {
      // The proxies must be destroyed before calling this.
      // DEBUG: Assert(this.proxyCount === 0);
    }

    /// Get the type of the child shape. You can use this to down cast to the concrete shape.
    /// @return the shape type.
    public getType(): ShapeType {
      return this.shape.getType();
    }

    /// Get the child shape. You can modify the child shape, however you should not change the
    /// number of vertices because this will crash some collision caching mechanisms.
    /// Manipulating the shape may lead to non-physical behavior.
    public getShape(): Shape {
      return this.shape;
    }

    /// Set if this fixture is a sensor.
    public setSensor(sensor: boolean): void {
      if (sensor !== this.isSensor) {
        this.body.setAwake(true);
        this.isSensor = sensor;
      }
    }

    /// Set the contact filtering data. This will not update contacts until the next time
    /// step when either parent body is active and awake.
    /// This automatically calls Refilter.
    public setFilterData(filter: Filter): void {
      this.filter.copy(filter);

      this.refilter();
    }

    /// Get the contact filtering data.
    public getFilterData(): Filter {
      return this.filter;
    }

    /// Call this if you want to establish collision that was previously disabled by ContactFilter::ShouldCollide.
    public refilter(): void {
      // Flag associated contacts for filtering.
      let edge = this.body.getContactList();

      while (edge) {
        const contact = edge.contact;
        const fixtureA = contact.getFixtureA();
        const fixtureB = contact.getFixtureB();
        if (fixtureA === this || fixtureB === this) {
          contact.glagForFiltering();
        }

        edge = edge.next;
      }

      // Touch each proxy so that new pairs may be created
      this.touchProxies();
    }

    /// Get the parent body of this fixture. This is NULL if the fixture is not attached.
    /// @return the parent body.
    public getBody(): Body {
      return this.body;
    }

    /// Test a point for containment in this fixture.
    /// @param p a point in world coordinates.
    public testPoint(p: XY): boolean {
      return this.shape.testPoint(this.body.getTransform(), p);
    }

    // #if ENABLE_PARTICLE
    public computeDistance(p: Vec2, normal: Vec2, childIndex: number): number {
      return this.shape.computeDistance(this.body.getTransform(), p, normal, childIndex);
    }
    // #endif

    /// Cast a ray against this shape.
    /// @param output the ray-cast results.
    /// @param input the ray-cast input parameters.
    public rayCast(output: RayCastOutput, input: RayCastInput, childIndex: number): boolean {
      return this.shape.rayCast(output, input, this.body.getTransform(), childIndex);
    }

    /// Get the mass data for this fixture. The mass data is based on the density and
    /// the shape. The rotational inertia is about the shape's origin. This operation
    /// may be expensive.
    public getMassData(massData: MassData = new MassData()): MassData {
      this.shape.computeMass(massData, this.density);

      return massData;
    }

    /// Set the density of this fixture. This will _not_ automatically adjust the mass
    /// of the body. You must call Body::ResetMassData to update the body's mass.
    public setDensity(density: number): void {
      this.density = density;
    }

    /// Get the density of this fixture.
    public getDensity(): number {
      return this.density;
    }

    /// Get the coefficient of friction.
    public getFriction(): number {
      return this.friction;
    }

    /// Set the coefficient of friction. This will _not_ change the friction of
    /// existing contacts.
    public setFriction(friction: number): void {
      this.friction = friction;
    }

    /// Get the coefficient of restitution.
    public getRestitution(): number {
      return this.restitution;
    }

    /// Set the coefficient of restitution. This will _not_ change the restitution of
    /// existing contacts.
    public setRestitution(restitution: number): void {
      this.restitution = restitution;
    }

    /// Get the restitution velocity threshold.
    public getRestitutionThreshold(): number {
      return this.restitutionThreshold;
    }

    /// Set the restitution threshold. This will _not_ change the restitution threshold of
    /// existing contacts.
    public setRestitutionThreshold(threshold: number): void {
      this.restitutionThreshold = threshold;
    }

    /// Get the fixture's AABB. This AABB may be enlarge and/or stale.
    /// If you need a more accurate AABB, compute it using the shape and
    /// the body transform.
    public getAABB(childIndex: number): AABB {
      // DEBUG: Assert(0 <= childIndex && childIndex < this.proxyCount);
      return this.proxies[childIndex].aabb;
    }

    /// Dump this fixture to the log file.
    public dump(log: (format: string, ...args: any[]) => void, bodyIndex: number): void {
      log("    const fd: FixtureDef = new FixtureDef();\n");
      log("    fd.friction = %.15f;\n", this.friction);
      log("    fd.restitution = %.15f;\n", this.restitution);
      log("    fd.restitutionThreshold = %.15f;\n", this.restitutionThreshold);
      log("    fd.density = %.15f;\n", this.density);
      log("    fd.isSensor = %s;\n", (this.isSensor) ? ("true") : ("false"));
      log("    fd.filter.categoryBits = %d;\n", this.filter.categoryBits);
      log("    fd.filter.maskBits = %d;\n", this.filter.maskBits);
      log("    fd.filter.groupIndex = %d;\n", this.filter.groupIndex);

      this.shape.dump(log);

      log("\n");
      log("    fd.shape = shape;\n");
      log("\n");
      log("    bodies[%d].CreateFixture(fd);\n", bodyIndex);
    }

    // These support body activation/deactivation.
    public createProxies(): void {
      if (this.proxies.length !== 0) { throw new Error(); }
      // Create proxies in the broad-phase.
      for (let i: number = 0; i < this.shape.getChildCount(); ++i) {
        this.proxies[i] = new FixtureProxy(this, i);
      }
    }

    public destroyProxies(): void {
      // Destroy proxies in the broad-phase.
      for (const proxy of this.proxies) {
        proxy.reset();
      }
      this.proxies.length = 0;
    }

    public touchProxies(): void {
      for (const proxy of this.proxies) {
        proxy.touch();
      }
    }

    public synchronizeProxies(transform1: Transform, transform2: Transform): void {
      for (const proxy of this.proxies) {
        proxy.synchronize(transform1, transform2);
      }
    }
  }

}
