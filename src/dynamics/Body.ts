/*
* Copyright (c) 2006-2011 Erin Catto http://www.box2d.org
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
  /// The body type.
/// static: zero mass, zero velocity, may be manually moved
/// kinematic: zero mass, non-zero velocity set by user, moved by solver
/// dynamic: positive mass, non-zero velocity determined by forces, moved by solver
  export enum BodyType {
    Unknown = -1,
    StaticBody = 0,
    KinematicBody = 1,
    DynamicBody = 2,

    // TODO_ERIN
    // bulletBody = 3
  }

  export interface IBodyDef {
    /// The body type: static, kinematic, or dynamic.
    /// Note: if a dynamic body would have zero mass, the mass is set to one.
    type?: BodyType;

    /// The world position of the body. Avoid creating bodies at the origin
    /// since this can lead to many overlapping shapes.
    position?: XY;

    /// The world angle of the body in radians.
    angle?: number;

    /// The linear velocity of the body's origin in world co-ordinates.
    linearVelocity?: XY;

    /// The angular velocity of the body.
    angularVelocity?: number;

    /// Linear damping is use to reduce the linear velocity. The damping parameter
    /// can be larger than 1.0f but the damping effect becomes sensitive to the
    /// time step when the damping parameter is large.
    /// Units are 1/time
    linearDamping?: number;

    /// Angular damping is use to reduce the angular velocity. The damping parameter
    /// can be larger than 1.0f but the damping effect becomes sensitive to the
    /// time step when the damping parameter is large.
    /// Units are 1/time
    angularDamping?: number;

    /// Set this flag to false if this body should never fall asleep. Note that
    /// this increases CPU usage.
    allowSleep?: boolean;

    /// Is this body initially awake or sleeping?
    awake?: boolean;

    /// Should this body be prevented from rotating? Useful for characters.
    fixedRotation?: boolean;

    /// Is this a fast moving body that should be prevented from tunneling through
    /// other moving bodies? Note that all bodies are prevented from tunneling through
    /// kinematic and static bodies. This setting is only considered on dynamic bodies.
    /// @warning You should use this flag sparingly since it increases processing time.
    bullet?: boolean;

    /// Does this body start out enabled?
    enabled?: boolean;

    /// Use this to store application specific body data.
    userData?: any;

    /// Scale the gravity applied to this body.
    gravityScale?: number;
  }

/// A body definition holds all the data needed to construct a rigid body.
/// You can safely re-use body definitions. Shapes are added to a body after construction.
  export class BodyDef implements IBodyDef {
    /// The body type: static, kinematic, or dynamic.
    /// Note: if a dynamic body would have zero mass, the mass is set to one.
    public type: BodyType = BodyType.StaticBody;

    /// The world position of the body. Avoid creating bodies at the origin
    /// since this can lead to many overlapping shapes.
    public readonly position: Vec2 = new Vec2(0, 0);

    /// The world angle of the body in radians.
    public angle: number = 0;

    /// The linear velocity of the body's origin in world co-ordinates.
    public readonly linearVelocity: Vec2 = new Vec2(0, 0);

    /// The angular velocity of the body.
    public angularVelocity: number = 0;

    /// Linear damping is use to reduce the linear velocity. The damping parameter
    /// can be larger than 1.0f but the damping effect becomes sensitive to the
    /// time step when the damping parameter is large.
    public linearDamping: number = 0;

    /// Angular damping is use to reduce the angular velocity. The damping parameter
    /// can be larger than 1.0f but the damping effect becomes sensitive to the
    /// time step when the damping parameter is large.
    public angularDamping: number = 0;

    /// Set this flag to false if this body should never fall asleep. Note that
    /// this increases CPU usage.
    public allowSleep: boolean = true;

    /// Is this body initially awake or sleeping?
    public awake: boolean = true;

    /// Should this body be prevented from rotating? Useful for characters.
    public fixedRotation: boolean = false;

    /// Is this a fast moving body that should be prevented from tunneling through
    /// other moving bodies? Note that all bodies are prevented from tunneling through
    /// kinematic and static bodies. This setting is only considered on dynamic bodies.
    /// @warning You should use this flag sparingly since it increases processing time.
    public bullet: boolean = false;

    /// Does this body start out enabled?
    public enabled: boolean = true;

    /// Use this to store application specific body data.
    public userData: any = null;

    /// Scale the gravity applied to this body.
    public gravityScale: number = 1;
  }

/// A rigid body. These are created via World::CreateBody.
  export class Body {
    public type: BodyType = BodyType.StaticBody;

    public islandFlag: boolean = false;
    public awakeFlag: boolean = false;
    public autoSleepFlag: boolean = false;
    public bulletFlag: boolean = false;
    public fixedRotationFlag: boolean = false;
    public enabledFlag: boolean = false;
    public toiFlag: boolean = false;

    public islandIndex: number = 0;

    public readonly xf: Transform = new Transform();  // the body origin transform
    // #if ENABLE_PARTICLE
    public readonly xf0: Transform = new Transform();
    // #endif
    public readonly sweep: Sweep = new Sweep();    // the swept motion for CCD

    public readonly linearVelocity: Vec2 = new Vec2();
    public angularVelocity: number = 0;

    public readonly force: Vec2 = new Vec2();
    public torque: number = 0;

    public world: World;
    public prev: Body = null;
    public next: Body = null;

    public fixtureList: Fixture = null;
    public fixtureCount: number = 0;

    public jointList: JointEdge = null;
    public contactList: ContactEdge = null;

    public mass: number = 1;
    public invMass: number = 1;

    // Rotational inertia about the center of mass.
    public I: number = 0;
    public invI: number = 0;

    public linearDamping: number = 0;
    public angularDamping: number = 0;
    public gravityScale: number = 1;

    public sleepTime: number = 0;

    public userData: any = null;

    // #if ENABLE_CONTROLLER
    public controllerList: ControllerEdge = null;
    public controllerCount: number = 0;
    // #endif

    constructor(bd: IBodyDef, world: World) {
      this.bulletFlag = maybe(bd.bullet, false);
      this.fixedRotationFlag = maybe(bd.fixedRotation, false);
      this.autoSleepFlag = maybe(bd.allowSleep, true);
      // this.awakeFlag = Maybe(bd.awake, true);
      if (maybe(bd.awake, false) && maybe(bd.type, BodyType.StaticBody) !== BodyType.StaticBody) {
        this.awakeFlag = true;
      }
      this.enabledFlag = maybe(bd.enabled, true);

      this.world = world;

      this.xf.p.copy(maybe(bd.position, Vec2.ZERO));
      // DEBUG: Assert(this.xf.p.IsValid());
      this.xf.q.setAngle(maybe(bd.angle, 0));
      // DEBUG: Assert(IsValid(this.xf.q.GetAngle()));
      // #if ENABLE_PARTICLE
      this.xf0.copy(this.xf);
      // #endif

      this.sweep.localCenter.setZero();
      this.sweep.c0.copy(this.xf.p);
      this.sweep.c.copy(this.xf.p);
      this.sweep.a0 = this.sweep.a = this.xf.q.getAngle();
      this.sweep.alpha0 = 0;

      this.linearVelocity.copy(maybe(bd.linearVelocity, Vec2.ZERO));
      // DEBUG: Assert(this.linearVelocity.IsValid());
      this.angularVelocity = maybe(bd.angularVelocity, 0);
      // DEBUG: Assert(IsValid(this.angularVelocity));

      this.linearDamping = maybe(bd.linearDamping, 0);
      this.angularDamping = maybe(bd.angularDamping, 0);
      this.gravityScale = maybe(bd.gravityScale, 1);
      // DEBUG: Assert(IsValid(this.gravityScale) && this.gravityScale >= 0);
      // DEBUG: Assert(IsValid(this.angularDamping) && this.angularDamping >= 0);
      // DEBUG: Assert(IsValid(this.linearDamping) && this.linearDamping >= 0);

      this.force.setZero();
      this.torque = 0;

      this.sleepTime = 0;

      this.type = maybe(bd.type, BodyType.StaticBody);

      this.mass = 0;
      this.invMass = 0;

      this.I = 0;
      this.invI = 0;

      this.userData = bd.userData;

      this.fixtureList = null;
      this.fixtureCount = 0;

      // #if ENABLE_CONTROLLER
      this.controllerList = null;
      this.controllerCount = 0;
      // #endif
    }

    public createFixture(def: IFixtureDef): Fixture;
    public createFixture(shape: Shape): Fixture;
    public createFixture(shape: Shape, density: number): Fixture;
    public createFixture(a: IFixtureDef | Shape, b: number = 0): Fixture {
      if (a instanceof Shape) {
        return this.createFixtureShapeDensity(a, b);
      } else {
        return this.createFixtureDef(a);
      }
    }

    public addFixture(fixture: Fixture): void {
      if (this.enabledFlag) {
        fixture.createProxies();
      }

      fixture.next = this.fixtureList;
      this.fixtureList = fixture;
      ++this.fixtureCount;

      // fixture.body = this;

      // Adjust mass properties if needed.
      if (fixture.density > 0) {
        this.resetMassData();
      }

      // Let the world know we have a new fixture. This will cause new contacts
      // to be created at the beginning of the next time step.
      this.world.newContacts = true;
    }

    /// Creates a fixture and attach it to this body. Use this function if you need
    /// to set some fixture parameters, like friction. Otherwise you can create the
    /// fixture directly from a shape.
    /// If the density is non-zero, this function automatically updates the mass of the body.
    /// Contacts are not created until the next time step.
    /// @param def the fixture definition.
    /// @warning This function is locked during callbacks.
    public createFixtureDef(def: IFixtureDef): Fixture {
      if (this.world.isLocked()) { throw new Error(); }
      const fixture: Fixture = new Fixture(this, def);
      this.addFixture(fixture);
      return fixture;
    }

    /// Creates a fixture from a shape and attach it to this body.
    /// This is a convenience function. Use FixtureDef if you need to set parameters
    /// like friction, restitution, user data, or filtering.
    /// If the density is non-zero, this function automatically updates the mass of the body.
    /// @param shape the shape to be cloned.
    /// @param density the shape density (set to zero for static bodies).
    /// @warning This function is locked during callbacks.
    private static createFixtureShapeDensity_s_def: FixtureDef = new FixtureDef();
    public createFixtureShapeDensity(shape: Shape, density: number = 0): Fixture {
      const def: FixtureDef = Body.createFixtureShapeDensity_s_def;
      def.shape = shape;
      def.density = density;
      return this.createFixtureDef(def);
    }

    /// Destroy a fixture. This removes the fixture from the broad-phase and
    /// destroys all contacts associated with this fixture. This will
    /// automatically adjust the mass of the body if the body is dynamic and the
    /// fixture has positive density.
    /// All fixtures attached to a body are implicitly destroyed when the body is destroyed.
    /// @param fixture the fixture to be removed.
    /// @warning This function is locked during callbacks.
    public destroyFixture(fixture: Fixture): void {
      if (this.world.isLocked()) { throw new Error(); }

      // DEBUG: Assert(fixture.body === this);

      // Remove the fixture from this body's singly linked list.
      // DEBUG: Assert(this.fixtureCount > 0);
      let node: Fixture = this.fixtureList;
      let ppF: Fixture = null;
      // DEBUG: let found: boolean = false;
      while (node !== null) {
        if (node === fixture) {
          if (ppF) {
            ppF.next = fixture.next;
          } else {
            this.fixtureList = fixture.next;
          }
          // DEBUG: found = true;
          break;
        }

        ppF = node;
        node = node.next;
      }

      // You tried to remove a shape that is not attached to this body.
      // DEBUG: Assert(found);

      // Destroy any contacts associated with the fixture.
      let edge: ContactEdge = this.contactList;
      while (edge) {
        const c = edge.contact;
        edge = edge.next;

        const fixtureA: Fixture = c.getFixtureA();
        const fixtureB: Fixture = c.getFixtureB();

        if (fixture === fixtureA || fixture === fixtureB) {
          // This destroys the contact and removes it from
          // this body's contact list.
          this.world.contactManager.destroy(c);
        }
      }

      if (this.enabledFlag) {
        fixture.destroyProxies();
      }

      // fixture.body = null;
      fixture.next = null;
      fixture.reset();

      --this.fixtureCount;

      // Reset the mass data.
      this.resetMassData();
    }

    /// Set the position of the body's origin and rotation.
    /// This breaks any contacts and wakes the other bodies.
    /// Manipulating a body's transform may cause non-physical behavior.
    /// @param position the world position of the body's local origin.
    /// @param angle the world rotation in radians.
    public setTransformVec(position: XY, angle: number): void {
      this.setTransformXY(position.x, position.y, angle);
    }

    public setTransformXY(x: number, y: number, angle: number): void {
      if (this.world.isLocked()) { throw new Error(); }

      this.xf.q.setAngle(angle);
      this.xf.p.set(x, y);
      // #if ENABLE_PARTICLE
      this.xf0.copy(this.xf);
      // #endif

      Transform.mulXV(this.xf, this.sweep.localCenter, this.sweep.c);
      this.sweep.a = angle;

      this.sweep.c0.copy(this.sweep.c);
      this.sweep.a0 = angle;

      for (let f: Fixture = this.fixtureList; f; f = f.next) {
        f.synchronizeProxies(this.xf, this.xf);
      }

      // Check for new contacts the next step
      this.world.newContacts = true;
    }

    public setTransform(xf: Transform): void {
      this.setTransformVec(xf.p, xf.getAngle());
    }

    /// Get the body transform for the body's origin.
    /// @return the world transform of the body's origin.
    public getTransform(): Transform {
      return this.xf;
    }

    /// Get the world body origin position.
    /// @return the world position of the body's origin.
    public getPosition(): Vec2 {
      return this.xf.p;
    }

    public setPosition(position: XY): void {
      this.setTransformVec(position, this.getAngle());
    }

    public setPositionXY(x: number, y: number): void {
      this.setTransformXY(x, y, this.getAngle());
    }

    /// Get the angle in radians.
    /// @return the current world rotation angle in radians.
    public getAngle(): number {
      return this.sweep.a;
    }

    public setAngle(angle: number): void {
      this.setTransformVec(this.getPosition(), angle);
    }

    /// Get the world position of the center of mass.
    public getWorldCenter(): Vec2 {
      return this.sweep.c;
    }

    /// Get the local position of the center of mass.
    public getLocalCenter(): Vec2 {
      return this.sweep.localCenter;
    }

    /// Set the linear velocity of the center of mass.
    /// @param v the new linear velocity of the center of mass.
    public setLinearVelocity(v: XY): void {
      if (this.type === BodyType.StaticBody) {
        return;
      }

      if (Vec2.DotVV(v, v) > 0) {
        this.setAwake(true);
      }

      this.linearVelocity.copy(v);
    }

    /// Get the linear velocity of the center of mass.
    /// @return the linear velocity of the center of mass.
    public getLinearVelocity(): Vec2 {
      return this.linearVelocity;
    }

    /// Set the angular velocity.
    /// @param omega the new angular velocity in radians/second.
    public setAngularVelocity(w: number): void {
      if (this.type === BodyType.StaticBody) {
        return;
      }

      if (w * w > 0) {
        this.setAwake(true);
      }

      this.angularVelocity = w;
    }

    /// Get the angular velocity.
    /// @return the angular velocity in radians/second.
    public getAngularVelocity(): number {
      return this.angularVelocity;
    }

    public getDefinition(bd: BodyDef): BodyDef {
      bd.type = this.getType();
      bd.allowSleep = this.autoSleepFlag;
      bd.angle = this.getAngle();
      bd.angularDamping = this.angularDamping;
      bd.gravityScale = this.gravityScale;
      bd.angularVelocity = this.angularVelocity;
      bd.fixedRotation = this.fixedRotationFlag;
      bd.bullet = this.bulletFlag;
      bd.awake = this.awakeFlag;
      bd.linearDamping = this.linearDamping;
      bd.linearVelocity.copy(this.getLinearVelocity());
      bd.position.copy(this.getPosition());
      bd.userData = this.getUserData();
      return bd;
    }

    /// Apply a force at a world point. If the force is not
    /// applied at the center of mass, it will generate a torque and
    /// affect the angular velocity. This wakes up the body.
    /// @param force the world force vector, usually in Newtons (N).
    /// @param point the world position of the point of application.
    /// @param wake also wake up the body
    public applyForce(force: XY, point: XY, wake: boolean = true): void {
      if (this.type !== BodyType.DynamicBody) {
        return;
      }

      if (wake && !this.awakeFlag) {
        this.setAwake(true);
      }

      // Don't accumulate a force if the body is sleeping.
      if (this.awakeFlag) {
        this.force.x += force.x;
        this.force.y += force.y;
        this.torque += ((point.x - this.sweep.c.x) * force.y - (point.y - this.sweep.c.y) * force.x);
      }
    }

    /// Apply a force to the center of mass. This wakes up the body.
    /// @param force the world force vector, usually in Newtons (N).
    /// @param wake also wake up the body
    public applyForceToCenter(force: XY, wake: boolean = true): void {
      if (this.type !== BodyType.DynamicBody) {
        return;
      }

      if (wake && !this.awakeFlag) {
        this.setAwake(true);
      }

      // Don't accumulate a force if the body is sleeping.
      if (this.awakeFlag) {
        this.force.x += force.x;
        this.force.y += force.y;
      }
    }

    /// Apply a torque. This affects the angular velocity
    /// without affecting the linear velocity of the center of mass.
    /// @param torque about the z-axis (out of the screen), usually in N-m.
    /// @param wake also wake up the body
    public applyTorque(torque: number, wake: boolean = true): void {
      if (this.type !== BodyType.DynamicBody) {
        return;
      }

      if (wake && !this.awakeFlag) {
        this.setAwake(true);
      }

      // Don't accumulate a force if the body is sleeping.
      if (this.awakeFlag) {
        this.torque += torque;
      }
    }

    /// Apply an impulse at a point. This immediately modifies the velocity.
    /// It also modifies the angular velocity if the point of application
    /// is not at the center of mass. This wakes up the body.
    /// @param impulse the world impulse vector, usually in N-seconds or kg-m/s.
    /// @param point the world position of the point of application.
    /// @param wake also wake up the body
    public applyLinearImpulse(impulse: XY, point: XY, wake: boolean = true): void {
      if (this.type !== BodyType.DynamicBody) {
        return;
      }

      if (wake && !this.awakeFlag) {
        this.setAwake(true);
      }

      // Don't accumulate a force if the body is sleeping.
      if (this.awakeFlag) {
        this.linearVelocity.x += this.invMass * impulse.x;
        this.linearVelocity.y += this.invMass * impulse.y;
        this.angularVelocity += this.invI * ((point.x - this.sweep.c.x) * impulse.y - (point.y - this.sweep.c.y) * impulse.x);
      }
    }

    /// Apply an impulse at the center of gravity. This immediately modifies the velocity.
    /// @param impulse the world impulse vector, usually in N-seconds or kg-m/s.
    /// @param wake also wake up the body
    public applyLinearImpulseToCenter(impulse: XY, wake: boolean = true): void {
      if (this.type !== BodyType.DynamicBody) {
        return;
      }

      if (wake && !this.awakeFlag) {
        this.setAwake(true);
      }

      // Don't accumulate a force if the body is sleeping.
      if (this.awakeFlag) {
        this.linearVelocity.x += this.invMass * impulse.x;
        this.linearVelocity.y += this.invMass * impulse.y;
      }
    }

    /// Apply an angular impulse.
    /// @param impulse the angular impulse in units of kg*m*m/s
    /// @param wake also wake up the body
    public applyAngularImpulse(impulse: number, wake: boolean = true): void {
      if (this.type !== BodyType.DynamicBody) {
        return;
      }

      if (wake && !this.awakeFlag) {
        this.setAwake(true);
      }

      // Don't accumulate a force if the body is sleeping.
      if (this.awakeFlag) {
        this.angularVelocity += this.invI * impulse;
      }
    }

    /// Get the total mass of the body.
    /// @return the mass, usually in kilograms (kg).
    public getMass(): number {
      return this.mass;
    }

    /// Get the rotational inertia of the body about the local origin.
    /// @return the rotational inertia, usually in kg-m^2.
    public getInertia(): number {
      return this.I + this.mass * Vec2.DotVV(this.sweep.localCenter, this.sweep.localCenter);
    }

    /// Get the mass data of the body.
    /// @return a struct containing the mass, inertia and center of the body.
    public getMassData(data: MassData): MassData {
      data.mass = this.mass;
      data.I = this.I + this.mass * Vec2.DotVV(this.sweep.localCenter, this.sweep.localCenter);
      data.center.copy(this.sweep.localCenter);
      return data;
    }

    /// Set the mass properties to override the mass properties of the fixtures.
    /// Note that this changes the center of mass position.
    /// Note that creating or destroying fixtures can also alter the mass.
    /// This function has no effect if the body isn't dynamic.
    /// @param massData the mass properties.
    private static setMassData_s_oldCenter: Vec2 = new Vec2();
    public setMassData(massData: MassData): void {
      if (this.world.isLocked()) { throw new Error(); }

      if (this.type !== BodyType.DynamicBody) {
        return;
      }

      this.invMass = 0;
      this.I = 0;
      this.invI = 0;

      this.mass = massData.mass;
      if (this.mass <= 0) {
        this.mass = 1;
      }

      this.invMass = 1 / this.mass;

      if (massData.I > 0 && !this.fixedRotationFlag) {
        this.I = massData.I - this.mass * Vec2.DotVV(massData.center, massData.center);
        // DEBUG: Assert(this.I > 0);
        this.invI = 1 / this.I;
      }

      // Move center of mass.
      const oldCenter: Vec2 = Body.setMassData_s_oldCenter.copy(this.sweep.c);
      this.sweep.localCenter.copy(massData.center);
      Transform.mulXV(this.xf, this.sweep.localCenter, this.sweep.c);
      this.sweep.c0.copy(this.sweep.c);

      // Update center of mass velocity.
      Vec2.AddVCrossSV(this.linearVelocity, this.angularVelocity, Vec2.SubVV(this.sweep.c, oldCenter, Vec2.s_t0), this.linearVelocity);
    }

    /// This resets the mass properties to the sum of the mass properties of the fixtures.
    /// This normally does not need to be called unless you called SetMassData to override
    /// the mass and you later want to reset the mass.
    private static resetMassData_s_localCenter: Vec2 = new Vec2();
    private static resetMassData_s_oldCenter: Vec2 = new Vec2();
    private static resetMassData_s_massData: MassData = new MassData();
    public resetMassData(): void {
      // Compute mass data from shapes. Each shape has its own density.
      this.mass = 0;
      this.invMass = 0;
      this.I = 0;
      this.invI = 0;
      this.sweep.localCenter.setZero();

      // Static and kinematic bodies have zero mass.
      if (this.type === BodyType.StaticBody || this.type === BodyType.KinematicBody) {
        this.sweep.c0.copy(this.xf.p);
        this.sweep.c.copy(this.xf.p);
        this.sweep.a0 = this.sweep.a;
        return;
      }

      // DEBUG: Assert(this.type === BodyType.dynamicBody);

      // Accumulate mass over all fixtures.
      const localCenter: Vec2 = Body.resetMassData_s_localCenter.setZero();
      for (let f: Fixture = this.fixtureList; f; f = f.next) {
        if (f.density === 0) {
          continue;
        }

        const massData: MassData = f.getMassData(Body.resetMassData_s_massData);
        this.mass += massData.mass;
        localCenter.x += massData.center.x * massData.mass;
        localCenter.y += massData.center.y * massData.mass;
        this.I += massData.I;
      }

      // Compute center of mass.
      if (this.mass > 0) {
        this.invMass = 1 / this.mass;
        localCenter.x *= this.invMass;
        localCenter.y *= this.invMass;
      }

      if (this.I > 0 && !this.fixedRotationFlag) {
        // Center the inertia about the center of mass.
        this.I -= this.mass * Vec2.DotVV(localCenter, localCenter);
        // DEBUG: Assert(this.I > 0);
        this.invI = 1 / this.I;
      } else {
        this.I = 0;
        this.invI = 0;
      }

      // Move center of mass.
      const oldCenter: Vec2 = Body.resetMassData_s_oldCenter.copy(this.sweep.c);
      this.sweep.localCenter.copy(localCenter);
      Transform.mulXV(this.xf, this.sweep.localCenter, this.sweep.c);
      this.sweep.c0.copy(this.sweep.c);

      // Update center of mass velocity.
      Vec2.AddVCrossSV(this.linearVelocity, this.angularVelocity, Vec2.SubVV(this.sweep.c, oldCenter, Vec2.s_t0), this.linearVelocity);
    }

    /// Get the world coordinates of a point given the local coordinates.
    /// @param localPoint a point on the body measured relative the the body's origin.
    /// @return the same point expressed in world coordinates.
    public getWorldPoint<T extends XY>(localPoint: XY, out: T): T {
      return Transform.mulXV(this.xf, localPoint, out);
    }

    /// Get the world coordinates of a vector given the local coordinates.
    /// @param localVector a vector fixed in the body.
    /// @return the same vector expressed in world coordinates.
    public getWorldVector<T extends XY>(localVector: XY, out: T): T {
      return Rot.mulRV(this.xf.q, localVector, out);
    }

    /// Gets a local point relative to the body's origin given a world point.
    /// @param a point in world coordinates.
    /// @return the corresponding local point relative to the body's origin.
    public getLocalPoint<T extends XY>(worldPoint: XY, out: T): T {
      return Transform.mulTXV(this.xf, worldPoint, out);
    }

    /// Gets a local vector given a world vector.
    /// @param a vector in world coordinates.
    /// @return the corresponding local vector.
    public getLocalVector<T extends XY>(worldVector: XY, out: T): T {
      return Rot.mulTRV(this.xf.q, worldVector, out);
    }

    /// Get the world linear velocity of a world point attached to this body.
    /// @param a point in world coordinates.
    /// @return the world velocity of a point.
    public getLinearVelocityFromWorldPoint<T extends XY>(worldPoint: XY, out: T): T {
      return Vec2.AddVCrossSV(this.linearVelocity, this.angularVelocity, Vec2.SubVV(worldPoint, this.sweep.c, Vec2.s_t0), out);
    }

    /// Get the world velocity of a local point.
    /// @param a point in local coordinates.
    /// @return the world velocity of a point.
    public getLinearVelocityFromLocalPoint<T extends XY>(localPoint: XY, out: T): T {
      return this.getLinearVelocityFromWorldPoint(this.getWorldPoint(localPoint, out), out);
    }

    /// Get the linear damping of the body.
    public getLinearDamping(): number {
      return this.linearDamping;
    }

    /// Set the linear damping of the body.
    public setLinearDamping(linearDamping: number): void {
      this.linearDamping = linearDamping;
    }

    /// Get the angular damping of the body.
    public getAngularDamping(): number {
      return this.angularDamping;
    }

    /// Set the angular damping of the body.
    public setAngularDamping(angularDamping: number): void {
      this.angularDamping = angularDamping;
    }

    /// Get the gravity scale of the body.
    public getGravityScale(): number {
      return this.gravityScale;
    }

    /// Set the gravity scale of the body.
    public setGravityScale(scale: number): void {
      this.gravityScale = scale;
    }

    /// Set the type of this body. This may alter the mass and velocity.
    public setType(type: BodyType): void {
      if (this.world.isLocked()) { throw new Error(); }

      if (this.type === type) {
        return;
      }

      this.type = type;

      this.resetMassData();

      if (this.type === BodyType.StaticBody) {
        this.linearVelocity.setZero();
        this.angularVelocity = 0;
        this.sweep.a0 = this.sweep.a;
        this.sweep.c0.copy(this.sweep.c);
        this.awakeFlag = false;
        this.synchronizeFixtures();
      }

      this.setAwake(true);

      this.force.setZero();
      this.torque = 0;

      // Delete the attached contacts.
      let ce: ContactEdge = this.contactList;
      while (ce) {
        const ce0: ContactEdge = ce;
        ce = ce.next;
        this.world.contactManager.destroy(ce0.contact);
      }
      this.contactList = null;

      // Touch the proxies so that new contacts will be created (when appropriate)
      for (let f: Fixture = this.fixtureList; f; f = f.next) {
        f.touchProxies();
      }
    }

    /// Get the type of this body.
    public getType(): BodyType {
      return this.type;
    }

    /// Should this body be treated like a bullet for continuous collision detection?
    public setBullet(flag: boolean): void {
      this.bulletFlag = flag;
    }

    /// Is this body treated like a bullet for continuous collision detection?
    public isBullet(): boolean {
      return this.bulletFlag;
    }

    /// You can disable sleeping on this body. If you disable sleeping, the
    /// body will be woken.
    public setSleepingAllowed(flag: boolean): void {
      this.autoSleepFlag = flag;
      if (!flag) {
        this.setAwake(true);
      }
    }

    /// Is this body allowed to sleep
    public isSleepingAllowed(): boolean {
      return this.autoSleepFlag;
    }

    /// Set the sleep state of the body. A sleeping body has very
    /// low CPU cost.
    /// @param flag set to true to wake the body, false to put it to sleep.
    public setAwake(flag: boolean): void {
      if (this.type === BodyType.StaticBody) {
        return;
      }
      if (flag) {
        this.awakeFlag = true;
        this.sleepTime = 0;
      } else {
        this.awakeFlag = false;
        this.sleepTime = 0;
        this.linearVelocity.setZero();
        this.angularVelocity = 0;
        this.force.setZero();
        this.torque = 0;
      }
    }

    /// Get the sleeping state of this body.
    /// @return true if the body is sleeping.
    public isAwake(): boolean {
      return this.awakeFlag;
    }

    /// Allow a body to be disabled. A disabled body is not simulated and cannot
    /// be collided with or woken up.
    /// If you pass a flag of true, all fixtures will be added to the broad-phase.
    /// If you pass a flag of false, all fixtures will be removed from the
    /// broad-phase and all contacts will be destroyed.
    /// Fixtures and joints are otherwise unaffected. You may continue
    /// to create/destroy fixtures and joints on disabled bodies.
    /// Fixtures on a disabled body are implicitly disabled and will
    /// not participate in collisions, ray-casts, or queries.
    /// Joints connected to a disabled body are implicitly disabled.
    /// An diabled body is still owned by a World object and remains
    /// in the body list.
    public setEnabled(flag: boolean): void {
      if (this.world.isLocked()) { throw new Error(); }

      if (flag === this.isEnabled()) {
        return;
      }

      this.enabledFlag = flag;

      if (flag) {
        // Create all proxies.
        for (let f: Fixture = this.fixtureList; f; f = f.next) {
          f.createProxies();
        }
        // Contacts are created at the beginning of the next
        this.world.newContacts = true;
      } else {
        // Destroy all proxies.
        for (let f: Fixture = this.fixtureList; f; f = f.next) {
          f.destroyProxies();
        }
        // Destroy the attached contacts.
        let ce: ContactEdge = this.contactList;
        while (ce) {
          const ce0: ContactEdge = ce;
          ce = ce.next;
          this.world.contactManager.destroy(ce0.contact);
        }
        this.contactList = null;
      }
    }

    /// Get the active state of the body.
    public isEnabled(): boolean {
      return this.enabledFlag;
    }

    /// Set this body to have fixed rotation. This causes the mass
    /// to be reset.
    public setFixedRotation(flag: boolean): void {
      if (this.fixedRotationFlag === flag) {
        return;
      }

      this.fixedRotationFlag = flag;

      this.angularVelocity = 0;

      this.resetMassData();
    }

    /// Does this body have fixed rotation?
    public isFixedRotation(): boolean {
      return this.fixedRotationFlag;
    }

    /// Get the list of all fixtures attached to this body.
    public getFixtureList(): Fixture {
      return this.fixtureList;
    }

    /// Get the list of all joints attached to this body.
    public getJointList(): JointEdge {
      return this.jointList;
    }

    /// Get the list of all contacts attached to this body.
    /// @warning this list changes during the time step and you may
    /// miss some collisions if you don't use ContactListener.
    public getContactList(): ContactEdge {
      return this.contactList;
    }

    /// Get the next body in the world's body list.
    public getNext(): Body {
      return this.next;
    }

    /// Get the user data pointer that was provided in the body definition.
    public getUserData(): any {
      return this.userData;
    }

    /// Set the user data. Use this to store your application specific data.
    public setUserData(data: any): void {
      this.userData = data;
    }

    /// Get the parent world of this body.
    public getWorld(): World {
      return this.world;
    }

    /// Dump this body to a file
    public dump(log: (format: string, ...args: any[]) => void): void {
      const bodyIndex: number = this.islandIndex;

      log("{\n");
      log("  const bd: BodyDef = new BodyDef();\n");
      let type_str: string = "";
      switch (this.type) {
        case BodyType.StaticBody:
          type_str = "BodyType.staticBody";
          break;
        case BodyType.KinematicBody:
          type_str = "BodyType.kinematicBody";
          break;
        case BodyType.DynamicBody:
          type_str = "BodyType.dynamicBody";
          break;
        default:
          // DEBUG: Assert(false);
          break;
      }
      log("  bd.type = %s;\n", type_str);
      log("  bd.position.Set(%.15f, %.15f);\n", this.xf.p.x, this.xf.p.y);
      log("  bd.angle = %.15f;\n", this.sweep.a);
      log("  bd.linearVelocity.Set(%.15f, %.15f);\n", this.linearVelocity.x, this.linearVelocity.y);
      log("  bd.angularVelocity = %.15f;\n", this.angularVelocity);
      log("  bd.linearDamping = %.15f;\n", this.linearDamping);
      log("  bd.angularDamping = %.15f;\n", this.angularDamping);
      log("  bd.allowSleep = %s;\n", (this.autoSleepFlag) ? ("true") : ("false"));
      log("  bd.awake = %s;\n", (this.awakeFlag) ? ("true") : ("false"));
      log("  bd.fixedRotation = %s;\n", (this.fixedRotationFlag) ? ("true") : ("false"));
      log("  bd.bullet = %s;\n", (this.bulletFlag) ? ("true") : ("false"));
      log("  bd.active = %s;\n", (this.enabledFlag) ? ("true") : ("false"));
      log("  bd.gravityScale = %.15f;\n", this.gravityScale);
      log("\n");
      log("  bodies[%d] = this.world.CreateBody(bd);\n", this.islandIndex);
      log("\n");
      for (let f: Fixture = this.fixtureList; f; f = f.next) {
        log("  {\n");
        f.dump(log, bodyIndex);
        log("  }\n");
      }
      log("}\n");
    }

    private static synchronizeFixtures_s_xf1: Transform = new Transform();
    public synchronizeFixtures(): void {
      if (this.awakeFlag) {
        const xf1: Transform = Body.synchronizeFixtures_s_xf1;
        xf1.q.setAngle(this.sweep.a0);
        Rot.mulRV(xf1.q, this.sweep.localCenter, xf1.p);
        Vec2.SubVV(this.sweep.c0, xf1.p, xf1.p);

        for (let f: Fixture = this.fixtureList; f; f = f.next) {
          f.synchronizeProxies(xf1, this.xf);
        }
      } else {
        for (let f: Fixture = this.fixtureList; f; f = f.next) {
          f.synchronizeProxies(this.xf, this.xf);
        }
      }
    }

    public synchronizeTransform(): void {
      this.xf.q.setAngle(this.sweep.a);
      Rot.mulRV(this.xf.q, this.sweep.localCenter, this.xf.p);
      Vec2.SubVV(this.sweep.c, this.xf.p, this.xf.p);
    }

    // This is used to prevent connected bodies from colliding.
    // It may lie, depending on the collideConnected flag.
    public shouldCollide(other: Body): boolean {
      // At least one body should be dynamic or kinematic.
      if (this.type === BodyType.StaticBody && other.type === BodyType.StaticBody) {
        return false;
      }
      return this.shouldCollideConnected(other);
    }

    public shouldCollideConnected(other: Body): boolean {
      // Does a joint prevent collision?
      for (let jn: JointEdge = this.jointList; jn; jn = jn.next) {
        if (jn.other === other) {
          if (!jn.joint.collideConnected) {
            return false;
          }
        }
      }

      return true;
    }

    public advance(alpha: number): void {
      // Advance to the new safe time. This doesn't sync the broad-phase.
      this.sweep.advance(alpha);
      this.sweep.c.copy(this.sweep.c0);
      this.sweep.a = this.sweep.a0;
      this.xf.q.setAngle(this.sweep.a);
      Rot.mulRV(this.xf.q, this.sweep.localCenter, this.xf.p);
      Vec2.SubVV(this.sweep.c, this.xf.p, this.xf.p);
    }

    // #if ENABLE_CONTROLLER
    public getControllerList(): ControllerEdge {
      return this.controllerList;
    }

    public getControllerCount(): number {
      return this.controllerCount;
    }
    // #endif
  }

}
