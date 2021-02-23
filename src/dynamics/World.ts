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
/// The world class manages all physics entities, dynamic simulation,
/// and asynchronous queries. The world also contains efficient memory
/// management facilities.
namespace b2 {
  export class World {
    public readonly contactManager: ContactManager = new ContactManager();

    public bodyList: Body = null;
    public jointList: Joint = null;

    // #if ENABLE_PARTICLE
    public particleSystemList: ParticleSystem = null;
    // #endif

    public bodyCount: number = 0;
    public jointCount: number = 0;

    public readonly gravity: Vec2 = new Vec2();
    public allowSleep: boolean = true;

    public destructionListener: DestructionListener;
    public debugDrawInstance: Draw;

    // This is used to compute the time step ratio to
    // support a variable time step.
    public inv_dt0: number = 0;

    public newContacts: boolean = false;
    public locked: boolean = false;
    public clearForcesFlag: boolean = true;

    // These are for debugging the solver.
    public warmStarting: boolean = true;
    public continuousPhysics: boolean = true;
    public subStepping: boolean = false;

    public stepComplete: boolean = true;

    public readonly profile: Profile = new Profile();

    public readonly island: Island = new Island();

    public readonly s_stack: Array<Body> = [];

    // #if ENABLE_CONTROLLER
    public controllerList: Controller = null;
    public controllerCount: number = 0;
    // #endif

    /// Construct a world object.
    /// @param gravity the world gravity vector.
    constructor(gravity: XY) {
      this.gravity.copy(gravity);
    }

    /// Register a destruction listener. The listener is owned by you and must
    /// remain in scope.
    public setDestructionListener(listener: DestructionListener): void {
      this.destructionListener = listener;
    }

    /// Register a contact filter to provide specific control over collision.
    /// Otherwise the default filter is used (defaultFilter). The listener is
    /// owned by you and must remain in scope.
    public setContactFilter(filter: ContactFilter): void {
      this.contactManager.contactFilter = filter;
    }

    /// Register a contact event listener. The listener is owned by you and must
    /// remain in scope.
    public setContactListener(listener: ContactListener): void {
      this.contactManager.contactListener = listener;
    }

    /// Register a routine for debug drawing. The debug draw functions are called
    /// inside with World::DebugDraw method. The debug draw object is owned
    /// by you and must remain in scope.
    public setDebugDraw(debugDraw: Draw): void {
      this.debugDrawInstance = debugDraw;
    }

    public addBody(b: Body): void {
      b.prev = null;
      b.next = this.bodyList;
      if (this.bodyList) {
        this.bodyList.prev = b;
      }
      this.bodyList = b;
      ++this.bodyCount;
    }

    /// Create a rigid body given a definition. No reference to the definition
    /// is retained.
    /// @warning This function is locked during callbacks.
    public createBody(def: IBodyDef = {}): Body {
      if (this.isLocked()) { throw new Error(); }
      const b: Body = new Body(def, this);
      // Add to world doubly linked list.
      this.addBody(b);
      return b;
    }

    /// Destroy a rigid body given a definition. No reference to the definition
    /// is retained. This function is locked during callbacks.
    /// @warning This automatically deletes all associated shapes and joints.
    /// @warning This function is locked during callbacks.
    public destroyBody(b: Body): void {
      // DEBUG: Assert(this.bodyCount > 0);
      if (this.isLocked()) { throw new Error(); }

      // Delete the attached joints.
      let je: JointEdge = b.jointList;
      while (je) {
        const je0: JointEdge = je;
        je = je.next;

        if (this.destructionListener) {
          this.destructionListener.sayGoodbyeJoint(je0.joint);
        }

        this.destroyJoint(je0.joint);

        b.jointList = je;
      }
      b.jointList = null;

      // #if ENABLE_CONTROLLER
      // @see Controller list
      let coe: ControllerEdge = b.controllerList;
      while (coe) {
        const coe0: ControllerEdge = coe;
        coe = coe.nextController;
        coe0.controller.removeBody(b);
      }
      // #endif

      // Delete the attached contacts.
      let ce: ContactEdge = b.contactList;
      while (ce) {
        const ce0: ContactEdge = ce;
        ce = ce.next;
        this.contactManager.destroy(ce0.contact);
      }
      b.contactList = null;

      // Delete the attached fixtures. This destroys broad-phase proxies.
      let f: Fixture = b.fixtureList;
      while (f) {
        const f0: Fixture = f;
        f = f.next;

        if (this.destructionListener) {
          this.destructionListener.sayGoodbyeFixture(f0);
        }

        f0.destroyProxies();
        f0.reset();

        b.fixtureList = f;
        b.fixtureCount -= 1;
      }
      b.fixtureList = null;
      b.fixtureCount = 0;

      // Remove world body list.
      if (b.prev) {
        b.prev.next = b.next;
      }

      if (b.next) {
        b.next.prev = b.prev;
      }

      if (b === this.bodyList) {
        this.bodyList = b.next;
      }

      --this.bodyCount;
    }

    private static _createJoint(def: IJointDef): Joint {
      switch (def.type) {
        case JointType.DistanceJoint: return new DistanceJoint(def as IDistanceJointDef);
        case JointType.MouseJoint: return new MouseJoint(def as IMouseJointDef);
        case JointType.PrismaticJoint: return new PrismaticJoint(def as IPrismaticJointDef);
        case JointType.RevoluteJoint: return new RevoluteJoint(def as IRevoluteJointDef);
        case JointType.PulleyJoint: return new PulleyJoint(def as IPulleyJointDef);
        case JointType.GearJoint: return new GearJoint(def as IGearJointDef);
        case JointType.WheelJoint: return new WheelJoint(def as IWheelJointDef);
        case JointType.WeldJoint: return new WeldJoint(def as IWeldJointDef);
        case JointType.FrictionJoint: return new FrictionJoint(def as IFrictionJointDef);
        case JointType.MotorJoint: return new MotorJoint(def as IMotorJointDef);
        case JointType.AreaJoint: return new AreaJoint(def as IAreaJointDef);
      }
      throw new Error();
    }

    private static _destroyJoint(joint: Joint): void {
    }

    /// Create a joint to constrain bodies together. No reference to the definition
    /// is retained. This may cause the connected bodies to cease colliding.
    /// @warning This function is locked during callbacks.
    public createJoint(def: IAreaJointDef): AreaJoint;
    public createJoint(def: IDistanceJointDef): DistanceJoint;
    public createJoint(def: IFrictionJointDef): FrictionJoint;
    public createJoint(def: IGearJointDef): GearJoint;
    public createJoint(def: IMotorJointDef): MotorJoint;
    public createJoint(def: IMouseJointDef): MouseJoint;
    public createJoint(def: IPrismaticJointDef): PrismaticJoint;
    public createJoint(def: IPulleyJointDef): PulleyJoint;
    public createJoint(def: IRevoluteJointDef): RevoluteJoint;
    public createJoint(def: IWeldJointDef): WeldJoint;
    public createJoint(def: IWheelJointDef): WheelJoint;
    public createJoint(def: IJointDef): Joint {
      if (this.isLocked()) { throw new Error(); }
      const j: Joint = World._createJoint(def);
      this.addJoint(j);
      // Note: creating a joint doesn't wake the bodies.
      return j;
    }

    public addJoint(j: Joint): void {
      // Connect to the world list.
      j.prev = null;
      j.next = this.jointList;
      if (this.jointList) {
        this.jointList.prev = j;
      }
      this.jointList = j;
      ++this.jointCount;

      // Connect to the bodies' doubly linked lists.
      // j.edgeA.other = j.bodyB; // done in Joint constructor
      j.edgeA.prev = null;
      j.edgeA.next = j.bodyA.jointList;
      if (j.bodyA.jointList) { j.bodyA.jointList.prev = j.edgeA; }
      j.bodyA.jointList = j.edgeA;

      // j.edgeB.other = j.bodyA; // done in Joint constructor
      j.edgeB.prev = null;
      j.edgeB.next = j.bodyB.jointList;
      if (j.bodyB.jointList) { j.bodyB.jointList.prev = j.edgeB; }
      j.bodyB.jointList = j.edgeB;

      const bodyA: Body = j.bodyA;
      const bodyB: Body = j.bodyB;
      const collideConnected: boolean = j.collideConnected;

      // If the joint prevents collisions, then flag any contacts for filtering.
      if (!collideConnected) {
        let edge: ContactEdge = bodyB.getContactList();
        while (edge) {
          if (edge.other === bodyA) {
            // Flag the contact for filtering at the next time step (where either
            // body is awake).
            edge.contact.glagForFiltering();
          }

          edge = edge.next;
        }
      }
    }

    /// Destroy a joint. This may cause the connected bodies to begin colliding.
    /// @warning This function is locked during callbacks.
    public destroyJoint(j: Joint): void {
      if (this.isLocked()) { throw new Error(); }

      // Remove from the doubly linked list.
      if (j.prev) {
        j.prev.next = j.next;
      }

      if (j.next) {
        j.next.prev = j.prev;
      }

      if (j === this.jointList) {
        this.jointList = j.next;
      }

      // Disconnect from island graph.
      const bodyA: Body = j.bodyA;
      const bodyB: Body = j.bodyB;
      const collideConnected: boolean = j.collideConnected;

      // Wake up connected bodies.
      bodyA.setAwake(true);
      bodyB.setAwake(true);

      // Remove from body 1.
      if (j.edgeA.prev) {
        j.edgeA.prev.next = j.edgeA.next;
      }

      if (j.edgeA.next) {
        j.edgeA.next.prev = j.edgeA.prev;
      }

      if (j.edgeA === bodyA.jointList) {
        bodyA.jointList = j.edgeA.next;
      }

      j.edgeA.reset();

      // Remove from body 2
      if (j.edgeB.prev) {
        j.edgeB.prev.next = j.edgeB.next;
      }

      if (j.edgeB.next) {
        j.edgeB.next.prev = j.edgeB.prev;
      }

      if (j.edgeB === bodyB.jointList) {
        bodyB.jointList = j.edgeB.next;
      }

      j.edgeB.reset();

      World._destroyJoint(j);

      // DEBUG: Assert(this.jointCount > 0);
      --this.jointCount;

      // If the joint prevents collisions, then flag any contacts for filtering.
      if (!collideConnected) {
        let edge: ContactEdge = bodyB.getContactList();
        while (edge) {
          if (edge.other === bodyA) {
            // Flag the contact for filtering at the next time step (where either
            // body is awake).
            edge.contact.glagForFiltering();
          }

          edge = edge.next;
        }
      }
    }

    // #if ENABLE_PARTICLE

    public createParticleSystem(def: ParticleSystemDef): ParticleSystem {
      if (this.isLocked()) { throw new Error(); }

      const p = new ParticleSystem(def, this);

      // Add to world doubly linked list.
      p.prev = null;
      p.next = this.particleSystemList;
      if (this.particleSystemList) {
        this.particleSystemList.prev = p;
      }
      this.particleSystemList = p;

      return p;
    }

    public destroyParticleSystem(p: ParticleSystem): void {
      if (this.isLocked()) { throw new Error(); }

      // Remove world particleSystem list.
      if (p.prev) {
        p.prev.next = p.next;
      }

      if (p.next) {
        p.next.prev = p.prev;
      }

      if (p === this.particleSystemList) {
        this.particleSystemList = p.next;
      }
    }

    public calculateReasonableParticleIterations(timeStep: number): number {
      if (this.particleSystemList === null) {
        return 1;
      }

      function GetSmallestRadius(world: World): number {
        let smallestRadius = maxFloat;
        for (let system = world.getParticleSystemList(); system !== null; system = system.next) {
          smallestRadius = Min(smallestRadius, system.getRadius());
        }
        return smallestRadius;
      }

      // Use the smallest radius, since that represents the worst-case.
      return CalculateParticleIterations(this.gravity.length(), GetSmallestRadius(this), timeStep);
    }

    // #endif

    /// Take a time step. This performs collision detection, integration,
    /// and constraint solution.
    /// @param timeStep the amount of time to simulate, this should not vary.
    /// @param velocityIterations for the velocity constraint solver.
    /// @param positionIterations for the position constraint solver.
    private static step_s_step = new TimeStep();
    private static step_s_stepTimer = new Timer();
    private static step_s_timer = new Timer();
    // #if ENABLE_PARTICLE
    public step(dt: number, velocityIterations: number, positionIterations: number, particleIterations: number = this.calculateReasonableParticleIterations(dt)): void {
      // #else
      // public Step(dt: number, velocityIterations: number, positionIterations: number): void {
      // #endif
      const stepTimer: Timer = World.step_s_stepTimer.reset();

      // If new fixtures were added, we need to find the new contacts.
      if (this.newContacts) {
        this.contactManager.findNewContacts();
        this.newContacts = false;
      }

      this.locked = true;

      const step: TimeStep = World.step_s_step;
      step.dt = dt;
      step.velocityIterations = velocityIterations;
      step.positionIterations = positionIterations;
      // #if ENABLE_PARTICLE
      step.particleIterations = particleIterations;
      // #endif
      if (dt > 0) {
        step.inv_dt = 1 / dt;
      } else {
        step.inv_dt = 0;
      }

      step.dtRatio = this.inv_dt0 * dt;

      step.warmStarting = this.warmStarting;

      // Update contacts. This is where some contacts are destroyed.
      const timer: Timer = World.step_s_timer.reset();
      this.contactManager.collide();
      this.profile.collide = timer.getMilliseconds();

      // Integrate velocities, solve velocity constraints, and integrate positions.
      if (this.stepComplete && step.dt > 0) {
        const timer: Timer = World.step_s_timer.reset();
        // #if ENABLE_PARTICLE
        for (let p = this.particleSystemList; p; p = p.next) {
          p.solve(step); // Particle Simulation
        }
        // #endif
        this.solve(step);
        this.profile.solve = timer.getMilliseconds();
      }

      // Handle TOI events.
      if (this.continuousPhysics && step.dt > 0) {
        const timer: Timer = World.step_s_timer.reset();
        this.solveTOI(step);
        this.profile.solveTOI = timer.getMilliseconds();
      }

      if (step.dt > 0) {
        this.inv_dt0 = step.inv_dt;
      }

      if (this.clearForcesFlag) {
        this.clearForces();
      }

      this.locked = false;

      this.profile.step = stepTimer.getMilliseconds();
    }

    /// Manually clear the force buffer on all bodies. By default, forces are cleared automatically
    /// after each call to Step. The default behavior is modified by calling SetAutoClearForces.
    /// The purpose of this function is to support sub-stepping. Sub-stepping is often used to maintain
    /// a fixed sized time step under a variable frame-rate.
    /// When you perform sub-stepping you will disable auto clearing of forces and instead call
    /// ClearForces after all sub-steps are complete in one pass of your game loop.
    /// @see SetAutoClearForces
    public clearForces(): void {
      for (let body = this.bodyList; body; body = body.next) {
        body.force.setZero();
        body.torque = 0;
      }
    }

    // #if ENABLE_PARTICLE

    public drawParticleSystem(system: ParticleSystem): void {
      if (this.debugDrawInstance === null) {
        return;
      }
      const particleCount = system.getParticleCount();
      if (particleCount) {
        const radius = system.getRadius();
        const positionBuffer = system.getPositionBuffer();
        if (system.colorBuffer.data) {
          const colorBuffer = system.getColorBuffer();
          this.debugDrawInstance.drawParticles(positionBuffer, radius, colorBuffer, particleCount);
        } else {
          this.debugDrawInstance.drawParticles(positionBuffer, radius, null, particleCount);
        }
      }
    }

    // #endif

    /// Call this to draw shapes and other debug draw data.
    private static debugdraw_s_color = new Color(0, 0, 0);
    private static debugdraw_s_vs = Vec2.MakeArray(4);
    private static debugdraw_s_xf = new Transform();
    public debugDraw(): void {
      if (this.debugDrawInstance === null) {
        return;
      }

      const flags: number = this.debugDrawInstance.getFlags();
      const color: Color = World.debugdraw_s_color.setRGB(0, 0, 0);

      if (flags & DrawFlags.ShapeBit) {
        for (let b: Body = this.bodyList; b; b = b.next) {
          const xf: Transform = b.xf;

          this.debugDrawInstance.pushTransform(xf);

          for (let f: Fixture = b.getFixtureList(); f; f = f.next) {
            if (b.getType() === BodyType.DynamicBody && b.mass === 0.0) {
              // Bad body
              this.drawShape(f, new Color(1.0, 0.0, 0.0));
            } else if (!b.isEnabled()) {
              color.setRGB(0.5, 0.5, 0.3);
              this.drawShape(f, color);
            } else if (b.getType() === BodyType.StaticBody) {
              color.setRGB(0.5, 0.9, 0.5);
              this.drawShape(f, color);
            } else if (b.getType() === BodyType.KinematicBody) {
              color.setRGB(0.5, 0.5, 0.9);
              this.drawShape(f, color);
            } else if (!b.isAwake()) {
              color.setRGB(0.6, 0.6, 0.6);
              this.drawShape(f, color);
            } else {
              color.setRGB(0.9, 0.7, 0.7);
              this.drawShape(f, color);
            }
          }

          this.debugDrawInstance.popTransform(xf);
        }
      }

      // #if ENABLE_PARTICLE
      if (flags & DrawFlags.ParticleBit) {
        for (let p = this.particleSystemList; p; p = p.next) {
          this.drawParticleSystem(p);
        }
      }
      // #endif

      if (flags & DrawFlags.JointBit) {
        for (let j: Joint = this.jointList; j; j = j.next) {
          j.draw(this.debugDrawInstance);
        }
      }

      if (flags & DrawFlags.PairBit) {
        color.setRGB(0.3, 0.9, 0.9);
        for (let contact = this.contactManager.contactList; contact; contact = contact.next) {
          const fixtureA = contact.getFixtureA();
          const fixtureB = contact.getFixtureB();
          const indexA = contact.getChildIndexA();
          const indexB = contact.getChildIndexB();
          const cA = fixtureA.getAABB(indexA).getCenter();
          const cB = fixtureB.getAABB(indexB).getCenter();

          this.debugDrawInstance.drawSegment(cA, cB, color);
        }
      }

      if (flags & DrawFlags.AABBBit) {
        color.setRGB(0.9, 0.3, 0.9);
        const vs: Vec2[] = World.debugdraw_s_vs;

        for (let b: Body = this.bodyList; b; b = b.next) {
          if (!b.isEnabled()) {
            continue;
          }

          for (let f: Fixture = b.getFixtureList(); f; f = f.next) {
            for (let i: number = 0; i < f.proxyCount; ++i) {
              const proxy: FixtureProxy = f.proxies[i];

              const aabb: AABB = proxy.treeNode.aabb;
              vs[0].set(aabb.lowerBound.x, aabb.lowerBound.y);
              vs[1].set(aabb.upperBound.x, aabb.lowerBound.y);
              vs[2].set(aabb.upperBound.x, aabb.upperBound.y);
              vs[3].set(aabb.lowerBound.x, aabb.upperBound.y);

              this.debugDrawInstance.drawPolygon(vs, 4, color);
            }
          }
        }
      }

      if (flags & DrawFlags.CenterOfMassBit) {
        for (let b: Body = this.bodyList; b; b = b.next) {
          const xf: Transform = World.debugdraw_s_xf;
          xf.q.copy(b.xf.q);
          xf.p.copy(b.getWorldCenter());
          this.debugDrawInstance.drawTransform(xf);
        }
      }

      // #if ENABLE_CONTROLLER
      // @see Controller list
      if (flags & DrawFlags.ControllerBit) {
        for (let c = this.controllerList; c; c = c.next) {
          c.draw(this.debugDrawInstance);
        }
      }
      // #endif
    }

    /// Query the world for all fixtures that potentially overlap the
    /// provided AABB.
    /// @param callback a user implemented callback class.
    /// @param aabb the query box.
    public queryAABB(callback: QueryCallback, aabb: AABB): void;
    public queryAABB(aabb: AABB, fn: QueryCallbackFunction): void;
    public queryAABB(...args: any[]): void {
      if (args[0] instanceof QueryCallback) {
        this._queryAABB(args[0], args[1]);
      } else {
        this._queryAABB(null, args[0], args[1]);
      }
    }
    private _queryAABB(callback: QueryCallback, aabb: AABB, fn?: QueryCallbackFunction): void {
      this.contactManager.broadPhase.query(aabb, (proxy: TreeNode<FixtureProxy>): boolean => {
        const fixture_proxy: FixtureProxy = proxy.userData;
        // DEBUG: Assert(fixture_proxy instanceof FixtureProxy);
        const fixture: Fixture = fixture_proxy.fixture;
        if (callback) {
          return callback.reportFixture(fixture);
        } else if (fn) {
          return fn(fixture);
        }
        return true;
      });
      // #if ENABLE_PARTICLE
      if (callback instanceof QueryCallback) {
        for (let p = this.particleSystemList; p; p = p.next) {
          if (callback.shouldQueryParticleSystem(p)) {
            p.queryAABB(callback, aabb);
          }
        }
      }
      // #endif
    }

    public queryAllAABB(aabb: AABB, out: Fixture[] = []): Fixture[] {
      this.queryAABB(aabb, (fixture: Fixture): boolean => { out.push(fixture); return true; });
      return out;
    }

    /// Query the world for all fixtures that potentially overlap the
    /// provided point.
    /// @param callback a user implemented callback class.
    /// @param point the query point.
    public queryPointAABB(callback: QueryCallback, point: XY): void;
    public queryPointAABB(point: XY, fn: QueryCallbackFunction): void;
    public queryPointAABB(...args: any[]): void {
      if (args[0] instanceof QueryCallback) {
        this._queryPointAABB(args[0], args[1]);
      } else {
        this._queryPointAABB(null, args[0], args[1]);
      }
    }
    private _queryPointAABB(callback: QueryCallback, point: XY, fn?: QueryCallbackFunction): void {
      this.contactManager.broadPhase.queryPoint(point, (proxy: TreeNode<FixtureProxy>): boolean => {
        const fixture_proxy: FixtureProxy = proxy.userData;
        // DEBUG: Assert(fixture_proxy instanceof FixtureProxy);
        const fixture: Fixture = fixture_proxy.fixture;
        if (callback) {
          return callback.reportFixture(fixture);
        } else if (fn) {
          return fn(fixture);
        }
        return true;
      });
      // #if ENABLE_PARTICLE
      if (callback instanceof QueryCallback) {
        for (let p = this.particleSystemList; p; p = p.next) {
          if (callback.shouldQueryParticleSystem(p)) {
            p.queryPointAABB(callback, point);
          }
        }
      }
      // #endif
    }

    public queryAllPointAABB(point: XY, out: Fixture[] = []): Fixture[] {
      this.queryPointAABB(point, (fixture: Fixture): boolean => { out.push(fixture); return true; });
      return out;
    }

    public queryFixtureShape(callback: QueryCallback, shape: Shape, index: number, transform: Transform): void;
    public queryFixtureShape(shape: Shape, index: number, transform: Transform, fn: QueryCallbackFunction): void;
    public queryFixtureShape(...args: any[]): void {
      if (args[0] instanceof QueryCallback) {
        this._queryFixtureShape(args[0], args[1], args[2], args[3]);
      } else {
        this._queryFixtureShape(null, args[0], args[1], args[2], args[3]);
      }
    }
    private static queryFixtureShape_s_aabb = new AABB();
    private _queryFixtureShape(callback: QueryCallback, shape: Shape, index: number, transform: Transform, fn?: QueryCallbackFunction): void {
      const aabb: AABB = World.queryFixtureShape_s_aabb;
      shape.computeAABB(aabb, transform, index);
      this.contactManager.broadPhase.query(aabb, (proxy: TreeNode<FixtureProxy>): boolean => {
        const fixture_proxy: FixtureProxy = proxy.userData;
        // DEBUG: Assert(fixture_proxy instanceof FixtureProxy);
        const fixture: Fixture = fixture_proxy.fixture;
        if (testOverlapShape(shape, index, fixture.getShape(), fixture_proxy.childIndex, transform, fixture.getBody().getTransform())) {
          if (callback) {
            return callback.reportFixture(fixture);
          } else if (fn) {
            return fn(fixture);
          }
        }
        return true;
      });
      // #if ENABLE_PARTICLE
      if (callback instanceof QueryCallback) {
        for (let p = this.particleSystemList; p; p = p.next) {
          if (callback.shouldQueryParticleSystem(p)) {
            p.queryAABB(callback, aabb);
          }
        }
      }
      // #endif
    }

    public queryAllFixtureShape(shape: Shape, index: number, transform: Transform, out: Fixture[] = []): Fixture[] {
      this.queryFixtureShape(shape, index, transform, (fixture: Fixture): boolean => { out.push(fixture); return true; });
      return out;
    }

    public queryFixturePoint(callback: QueryCallback, point: XY): void;
    public queryFixturePoint(point: XY, fn: QueryCallbackFunction): void;
    public queryFixturePoint(...args: any[]): void {
      if (args[0] instanceof QueryCallback) {
        this._queryFixturePoint(args[0], args[1]);
      } else {
        this._queryFixturePoint(null, args[0], args[1]);
      }
    }
    private _queryFixturePoint(callback: QueryCallback, point: XY, fn?: QueryCallbackFunction): void {
      this.contactManager.broadPhase.queryPoint(point, (proxy: TreeNode<FixtureProxy>): boolean => {
        const fixture_proxy: FixtureProxy = proxy.userData;
        // DEBUG: Assert(fixture_proxy instanceof FixtureProxy);
        const fixture: Fixture = fixture_proxy.fixture;
        if (fixture.testPoint(point)) {
          if (callback) {
            return callback.reportFixture(fixture);
          } else if (fn) {
            return fn(fixture);
          }
        }
        return true;
      });
      // #if ENABLE_PARTICLE
      if (callback) {
        for (let p = this.particleSystemList; p; p = p.next) {
          if (callback.shouldQueryParticleSystem(p)) {
            p.queryPointAABB(callback, point);
          }
        }
      }
      // #endif
    }

    public queryAllFixturePoint(point: XY, out: Fixture[] = []): Fixture[] {
      this.queryFixturePoint(point, (fixture: Fixture): boolean => { out.push(fixture); return true; });
      return out;
    }

    /// Ray-cast the world for all fixtures in the path of the ray. Your callback
    /// controls whether you get the closest point, any point, or n-points.
    /// The ray-cast ignores shapes that contain the starting point.
    /// @param callback a user implemented callback class.
    /// @param point1 the ray starting point
    /// @param point2 the ray ending point
    public rayCast(callback: RayCastCallback, point1: XY, point2: XY): void;
    public rayCast(point1: XY, point2: XY, fn: RayCastCallbackFunction): void;
    public rayCast(...args: any[]): void {
      if (args[0] instanceof RayCastCallback) {
        this._rayCast(args[0], args[1], args[2]);
      } else {
        this._rayCast(null, args[0], args[1], args[2]);
      }
    }
    private static rayCast_s_input = new RayCastInput();
    private static rayCast_s_output = new RayCastOutput();
    private static rayCast_s_point = new Vec2();
    private _rayCast(callback: RayCastCallback, point1: XY, point2: XY, fn?: RayCastCallbackFunction): void {
      const input: RayCastInput = World.rayCast_s_input;
      input.maxFraction = 1;
      input.p1.copy(point1);
      input.p2.copy(point2);
      this.contactManager.broadPhase.rayCast(input, (input: RayCastInput, proxy: TreeNode<FixtureProxy>): number => {
        const fixture_proxy: FixtureProxy = proxy.userData;
        // DEBUG: Assert(fixture_proxy instanceof FixtureProxy);
        const fixture: Fixture = fixture_proxy.fixture;
        const index: number = fixture_proxy.childIndex;
        const output: RayCastOutput = World.rayCast_s_output;
        const hit: boolean = fixture.rayCast(output, input, index);
        if (hit) {
          const fraction: number = output.fraction;
          const point: Vec2 = World.rayCast_s_point;
          point.set((1 - fraction) * point1.x + fraction * point2.x, (1 - fraction) * point1.y + fraction * point2.y);
          if (callback) {
            return callback.reportFixture(fixture, point, output.normal, fraction);
          } else if (fn) {
            return fn(fixture, point, output.normal, fraction);
          }
        }
        return input.maxFraction;
      });
      // #if ENABLE_PARTICLE
      if (callback) {
        for (let p = this.particleSystemList; p; p = p.next) {
          if (callback.shouldQueryParticleSystem(p)) {
            p.rayCast(callback, point1, point2);
          }
        }
      }
      // #endif
    }

    public rayCastOne(point1: XY, point2: XY): Fixture {
      let result: Fixture = null;
      let min_fraction: number = 1;
      this.rayCast(point1, point2, (fixture: Fixture, point: Vec2, normal: Vec2, fraction: number): number => {
        if (fraction < min_fraction) {
          min_fraction = fraction;
          result = fixture;
        }
        return min_fraction;
      });
      return result;
    }

    public rayCastAll(point1: XY, point2: XY, out: Fixture[] = []): Fixture[] {
      this.rayCast(point1, point2, (fixture: Fixture, point: Vec2, normal: Vec2, fraction: number): number => {
        out.push(fixture);
        return 1;
      });
      return out;
    }

    // #if ENABLE_PARTICLE
    public getParticleSystemList(): ParticleSystem {
      return this.particleSystemList;
    }
    // #endif

    /// Get the world contact list. With the returned contact, use Contact::GetNext to get
    /// the next contact in the world list. A NULL contact indicates the end of the list.
    /// @return the head of the world contact list.
    /// @warning contacts are created and destroyed in the middle of a time step.
    /// Use ContactListener to avoid missing contacts.
    public getContactList(): Contact {
      return this.contactManager.contactList;
    }

    /// Enable/disable sleep.
    public setAllowSleeping(flag: boolean): void {
      if (flag === this.allowSleep) {
        return;
      }

      this.allowSleep = flag;
      if (!this.allowSleep) {
        for (let b = this.bodyList; b; b = b.next) {
          b.setAwake(true);
        }
      }
    }

    public getAllowSleeping(): boolean {
      return this.allowSleep;
    }

    /// Enable/disable warm starting. For testing.
    public setWarmStarting(flag: boolean): void {
      this.warmStarting = flag;
    }

    public getWarmStarting(): boolean {
      return this.warmStarting;
    }

    /// Enable/disable continuous physics. For testing.
    public setContinuousPhysics(flag: boolean): void {
      this.continuousPhysics = flag;
    }

    public getContinuousPhysics(): boolean {
      return this.continuousPhysics;
    }

    /// Enable/disable single stepped continuous physics. For testing.
    public setSubStepping(flag: boolean): void {
      this.subStepping = flag;
    }

    public getSubStepping(): boolean {
      return this.subStepping;
    }

    /// Get the number of broad-phase proxies.
    public getProxyCount(): number {
      return this.contactManager.broadPhase.getProxyCount();
    }

    /// Get the number of bodies.
    public getBodyCount(): number {
      return this.bodyCount;
    }

    /// Get the number of joints.
    public getJointCount(): number {
      return this.jointCount;
    }

    /// Get the number of contacts (each may have 0 or more contact points).
    public getContactCount(): number {
      return this.contactManager.contactCount;
    }

    /// Get the height of the dynamic tree.
    public getTreeHeight(): number {
      return this.contactManager.broadPhase.getTreeHeight();
    }

    /// Get the balance of the dynamic tree.
    public getTreeBalance(): number {
      return this.contactManager.broadPhase.getTreeBalance();
    }

    /// Get the quality metric of the dynamic tree. The smaller the better.
    /// The minimum is 1.
    public getTreeQuality(): number {
      return this.contactManager.broadPhase.getTreeQuality();
    }

    /// Change the global gravity vector.
    public setGravity(gravity: XY, wake: boolean = true) {
      if (!Vec2.IsEqualToV(this.gravity, gravity)) {
        this.gravity.copy(gravity);

        if (wake) {
          for (let b: Body = this.bodyList; b; b = b.next) {
            b.setAwake(true);
          }
        }
      }
    }

    /// Is the world locked (in the middle of a time step).
    public isLocked(): boolean {
      return this.locked;
    }

    /// Set flag to control automatic clearing of forces after each time step.
    public setAutoClearForces(flag: boolean): void {
      this.clearForcesFlag = flag;
    }

    /// Get the flag that controls automatic clearing of forces after each time step.
    public getAutoClearForces(): boolean {
      return this.clearForcesFlag;
    }

    /// Shift the world origin. Useful for large worlds.
    /// The body shift formula is: position -= newOrigin
    /// @param newOrigin the new origin with respect to the old origin
    public shiftOrigin(newOrigin: XY): void {
      if (this.isLocked()) { throw new Error(); }

      for (let b: Body = this.bodyList; b; b = b.next) {
        b.xf.p.selfSub(newOrigin);
        b.sweep.c0.selfSub(newOrigin);
        b.sweep.c.selfSub(newOrigin);
      }

      for (let j: Joint = this.jointList; j; j = j.next) {
        j.shiftOrigin(newOrigin);
      }

      this.contactManager.broadPhase.shiftOrigin(newOrigin);
    }

    /// Get the contact manager for testing.
    public getContactManager(): ContactManager {
      return this.contactManager;
    }

    /// Get the current profile.
    public getProfile(): Profile {
      return this.profile;
    }

    /// Dump the world into the log file.
    /// @warning this should be called outside of a time step.
    public dump(log: (format: string, ...args: any[]) => void): void {
      if (this.locked) {
        return;
      }

      // OpenDump("box2d_dump.inl");

      log("const g: Vec2 = new Vec2(%.15f, %.15f);\n", this.gravity.x, this.gravity.y);
      log("this.world.SetGravity(g);\n");

      log("const bodies: Body[] = [];\n");
      log("const joints: Joint[] = [];\n");
      let i: number = 0;
      for (let b: Body = this.bodyList; b; b = b.next) {
        b.islandIndex = i;
        b.dump(log);
        ++i;
      }

      i = 0;
      for (let j: Joint = this.jointList; j; j = j.next) {
        j.index = i;
        ++i;
      }

      // First pass on joints, skip gear joints.
      for (let j: Joint = this.jointList; j; j = j.next) {
        if (j.type === JointType.GearJoint) {
          continue;
        }

        log("{\n");
        j.dump(log);
        log("}\n");
      }

      // Second pass on joints, only gear joints.
      for (let j: Joint = this.jointList; j; j = j.next) {
        if (j.type !== JointType.GearJoint) {
          continue;
        }

        log("{\n");
        j.dump(log);
        log("}\n");
      }

      // CloseDump();
    }

    public drawShape(fixture: Fixture, color: Color): void {
      if (this.debugDrawInstance === null) {
        return;
      }
      const shape: Shape = fixture.getShape();

      switch (shape.type) {
        case ShapeType.CircleShape: {
          const circle: CircleShape = shape as CircleShape;
          const center: Vec2 = circle.p;
          const radius: number = circle.radius;
          const axis: Vec2 = Vec2.UNITX;
          this.debugDrawInstance.drawSolidCircle(center, radius, axis, color);
          break;
        }

        case ShapeType.EdgeShape: {
          const edge: EdgeShape = shape as EdgeShape;
          const v1: Vec2 = edge.vertex1;
          const v2: Vec2 = edge.vertex2;
          this.debugDrawInstance.drawSegment(v1, v2, color);

          if (edge.oneSided === false) {
            this.debugDrawInstance.drawPoint(v1, 4.0, color);
            this.debugDrawInstance.drawPoint(v2, 4.0, color);
          }
          break;
        }

        case ShapeType.ChainShape: {
          const chain: ChainShape = shape as ChainShape;
          const count: number = chain.count;
          const vertices: Vec2[] = chain.vertices;
          let v1: Vec2 = vertices[0];
          for (let i: number = 1; i < count; ++i) {
            const v2: Vec2 = vertices[i];
            this.debugDrawInstance.drawSegment(v1, v2, color);
            v1 = v2;
          }

          break;
        }

        case ShapeType.PolygonShape: {
          const poly: PolygonShape = shape as PolygonShape;
          const vertexCount: number = poly.count;
          const vertices: Vec2[] = poly.vertices;
          this.debugDrawInstance.drawSolidPolygon(vertices, vertexCount, color);
          break;
        }
      }
    }

    public solve(step: TimeStep): void {
      // #if ENABLE_PARTICLE
      // update previous transforms
      for (let b = this.bodyList; b; b = b.next) {
        b.xf0.copy(b.xf);
      }
      // #endif

      // #if ENABLE_CONTROLLER
      // @see Controller list
      for (let controller = this.controllerList; controller; controller = controller.next) {
        controller.step(step);
      }
      // #endif

      this.profile.solveInit = 0;
      this.profile.solveVelocity = 0;
      this.profile.solvePosition = 0;

      // Size the island for the worst case.
      const island: Island = this.island;
      island.initialize(this.bodyCount,
        this.contactManager.contactCount,
        this.jointCount,
        this.contactManager.contactListener);

      // Clear all the island flags.
      for (let b: Body = this.bodyList; b; b = b.next) {
        b.islandFlag = false;
      }
      for (let c: Contact = this.contactManager.contactList; c; c = c.next) {
        c.islandFlag = false;
      }
      for (let j: Joint = this.jointList; j; j = j.next) {
        j.islandFlag = false;
      }

      // Build and simulate all awake islands.
      // DEBUG: const stackSize: number = this.bodyCount;
      const stack: Array<Body> = this.s_stack;
      for (let seed: Body = this.bodyList; seed; seed = seed.next) {
        if (seed.islandFlag) {
          continue;
        }

        if (!seed.isAwake() || !seed.isEnabled()) {
          continue;
        }

        // The seed can be dynamic or kinematic.
        if (seed.getType() === BodyType.StaticBody) {
          continue;
        }

        // Reset island and stack.
        island.clear();
        let stackCount: number = 0;
        stack[stackCount++] = seed;
        seed.islandFlag = true;

        // Perform a depth first search (DFS) on the constraint graph.
        while (stackCount > 0) {
          // Grab the next body off the stack and add it to the island.
          const b: Body = stack[--stackCount];
          if (!b) { throw new Error(); }
          // DEBUG: Assert(b.IsEnabled());
          island.addBody(b);

          // To keep islands as small as possible, we don't
          // propagate islands across static bodies.
          if (b.getType() === BodyType.StaticBody) {
            continue;
          }

          // Make sure the body is awake. (without resetting sleep timer).
          b.awakeFlag = true;

          // Search all contacts connected to this body.
          for (let ce: ContactEdge = b.contactList; ce; ce = ce.next) {
            const contact: Contact = ce.contact;

            // Has this contact already been added to an island?
            if (contact.islandFlag) {
              continue;
            }

            // Is this contact solid and touching?
            if (!contact.isEnabled() || !contact.isTouching()) {
              continue;
            }

            // Skip sensors.
            const sensorA: boolean = contact.fixtureA.isSensor;
            const sensorB: boolean = contact.fixtureB.isSensor;
            if (sensorA || sensorB) {
              continue;
            }

            island.addContact(contact);
            contact.islandFlag = true;

            const other: Body = ce.other;

            // Was the other body already added to this island?
            if (other.islandFlag) {
              continue;
            }

            // DEBUG: Assert(stackCount < stackSize);
            stack[stackCount++] = other;
            other.islandFlag = true;
          }

          // Search all joints connect to this body.
          for (let je: JointEdge = b.jointList; je; je = je.next) {
            if (je.joint.islandFlag) {
              continue;
            }

            const other: Body = je.other;

            // Don't simulate joints connected to disabled bodies.
            if (!other.isEnabled()) {
              continue;
            }

            island.addJoint(je.joint);
            je.joint.islandFlag = true;

            if (other.islandFlag) {
              continue;
            }

            // DEBUG: Assert(stackCount < stackSize);
            stack[stackCount++] = other;
            other.islandFlag = true;
          }
        }

        const profile: Profile = new Profile();
        island.solve(profile, step, this.gravity, this.allowSleep);
        this.profile.solveInit += profile.solveInit;
        this.profile.solveVelocity += profile.solveVelocity;
        this.profile.solvePosition += profile.solvePosition;

        // Post solve cleanup.
        for (let i: number = 0; i < island.bodyCount; ++i) {
          // Allow static bodies to participate in other islands.
          const b: Body = island.bodies[i];
          if (b.getType() === BodyType.StaticBody) {
            b.islandFlag = false;
          }
        }
      }

      for (let i: number = 0; i < stack.length; ++i) {
        if (!stack[i]) { break; }
        stack[i] = null;
      }

      const timer: Timer = new Timer();

      // Synchronize fixtures, check for out of range bodies.
      for (let b = this.bodyList; b; b = b.next) {
        // If a body was not in an island then it did not move.
        if (!b.islandFlag) {
          continue;
        }

        if (b.getType() === BodyType.StaticBody) {
          continue;
        }

        // Update fixtures (for broad-phase).
        b.synchronizeFixtures();
      }

      // Look for new contacts.
      this.contactManager.findNewContacts();
      this.profile.broadphase = timer.getMilliseconds();
    }

    private static solveTOI_s_subStep = new TimeStep();
    private static solveTOI_s_backup = new Sweep();
    private static solveTOI_s_backup1 = new Sweep();
    private static solveTOI_s_backup2 = new Sweep();
    private static solveTOI_s_toi_input = new TOIInput();
    private static solveTOI_s_toi_output = new TOIOutput();
    public solveTOI(step: TimeStep): void {
      const island: Island = this.island;
      island.initialize(2 * maxTOIContacts, maxTOIContacts, 0, this.contactManager.contactListener);

      if (this.stepComplete) {
        for (let b: Body = this.bodyList; b; b = b.next) {
          b.islandFlag = false;
          b.sweep.alpha0 = 0;
        }

        for (let c: Contact = this.contactManager.contactList; c; c = c.next) {
          // Invalidate TOI
          c.toiFlag = false;
          c.islandFlag = false;
          c.toiCount = 0;
          c.toi = 1;
        }
      }

      // Find TOI events and solve them.
      for (; ;) {
        // Find the first TOI.
        let minContact: Contact = null;
        let minAlpha: number = 1;

        for (let c: Contact = this.contactManager.contactList; c; c = c.next) {
          // Is this contact disabled?
          if (!c.isEnabled()) {
            continue;
          }

          // Prevent excessive sub-stepping.
          if (c.toiCount > maxSubSteps) {
            continue;
          }

          let alpha: number = 1;
          if (c.toiFlag) {
            // This contact has a valid cached TOI.
            alpha = c.toi;
          } else {
            const fA: Fixture = c.getFixtureA();
            const fB: Fixture = c.getFixtureB();

            // Is there a sensor?
            if (fA.isSensor || fB.isSensor) {
              continue;
            }

            const bA: Body = fA.getBody();
            const bB: Body = fB.getBody();

            const typeA: BodyType = bA.type;
            const typeB: BodyType = bB.type;
            // DEBUG: Assert(typeA !== BodyType.staticBody || typeB !== BodyType.staticBody);

            const activeA: boolean = bA.isAwake() && typeA !== BodyType.StaticBody;
            const activeB: boolean = bB.isAwake() && typeB !== BodyType.StaticBody;

            // Is at least one body active (awake and dynamic or kinematic)?
            if (!activeA && !activeB) {
              continue;
            }

            const collideA: boolean = bA.isBullet() || typeA !== BodyType.DynamicBody;
            const collideB: boolean = bB.isBullet() || typeB !== BodyType.DynamicBody;

            // Are these two non-bullet dynamic bodies?
            if (!collideA && !collideB) {
              continue;
            }

            // Compute the TOI for this contact.
            // Put the sweeps onto the same time interval.
            let alpha0: number = bA.sweep.alpha0;

            if (bA.sweep.alpha0 < bB.sweep.alpha0) {
              alpha0 = bB.sweep.alpha0;
              bA.sweep.advance(alpha0);
            } else if (bB.sweep.alpha0 < bA.sweep.alpha0) {
              alpha0 = bA.sweep.alpha0;
              bB.sweep.advance(alpha0);
            }

            // DEBUG: Assert(alpha0 < 1);

            const indexA: number = c.getChildIndexA();
            const indexB: number = c.getChildIndexB();

            // Compute the time of impact in interval [0, minTOI]
            const input: TOIInput = World.solveTOI_s_toi_input;
            input.proxyA.setShape(fA.getShape(), indexA);
            input.proxyB.setShape(fB.getShape(), indexB);
            input.sweepA.copy(bA.sweep);
            input.sweepB.copy(bB.sweep);
            input.tMax = 1;

            const output: TOIOutput = World.solveTOI_s_toi_output;
            timeOfImpact(output, input);

            // Beta is the fraction of the remaining portion of the .
            const beta: number = output.t;
            if (output.state === TOIOutputState.Touching) {
              alpha = Min(alpha0 + (1 - alpha0) * beta, 1);
            } else {
              alpha = 1;
            }

            c.toi = alpha;
            c.toiFlag = true;
          }

          if (alpha < minAlpha) {
            // This is the minimum TOI found so far.
            minContact = c;
            minAlpha = alpha;
          }
        }

        if (minContact === null || 1 - 10 * epsilon < minAlpha) {
          // No more TOI events. Done!
          this.stepComplete = true;
          break;
        }

        // Advance the bodies to the TOI.
        const fA: Fixture = minContact.getFixtureA();
        const fB: Fixture = minContact.getFixtureB();
        const bA: Body = fA.getBody();
        const bB: Body = fB.getBody();

        const backup1: Sweep = World.solveTOI_s_backup1.copy(bA.sweep);
        const backup2: Sweep = World.solveTOI_s_backup2.copy(bB.sweep);

        bA.advance(minAlpha);
        bB.advance(minAlpha);

        // The TOI contact likely has some new contact points.
        minContact.update(this.contactManager.contactListener);
        minContact.toiFlag = false;
        ++minContact.toiCount;

        // Is the contact solid?
        if (!minContact.isEnabled() || !minContact.isTouching()) {
          // Restore the sweeps.
          minContact.setEnabled(false);
          bA.sweep.copy(backup1);
          bB.sweep.copy(backup2);
          bA.synchronizeTransform();
          bB.synchronizeTransform();
          continue;
        }

        bA.setAwake(true);
        bB.setAwake(true);

        // Build the island
        island.clear();
        island.addBody(bA);
        island.addBody(bB);
        island.addContact(minContact);

        bA.islandFlag = true;
        bB.islandFlag = true;
        minContact.islandFlag = true;

        // Get contacts on bodyA and bodyB.
        // const bodies: Body[] = [bA, bB];
        for (let i: number = 0; i < 2; ++i) {
          const body: Body = (i === 0) ? (bA) : (bB); // bodies[i];
          if (body.type === BodyType.DynamicBody) {
            for (let ce: ContactEdge = body.contactList; ce; ce = ce.next) {
              if (island.bodyCount === island.bodyCapacity) {
                break;
              }

              if (island.contactCount === island.contactCapacity) {
                break;
              }

              const contact: Contact = ce.contact;

              // Has this contact already been added to the island?
              if (contact.islandFlag) {
                continue;
              }

              // Only add static, kinematic, or bullet bodies.
              const other: Body = ce.other;
              if (other.type === BodyType.DynamicBody &&
                !body.isBullet() && !other.isBullet()) {
                continue;
              }

              // Skip sensors.
              const sensorA: boolean = contact.fixtureA.isSensor;
              const sensorB: boolean = contact.fixtureB.isSensor;
              if (sensorA || sensorB) {
                continue;
              }

              // Tentatively advance the body to the TOI.
              const backup: Sweep = World.solveTOI_s_backup.copy(other.sweep);
              if (!other.islandFlag) {
                other.advance(minAlpha);
              }

              // Update the contact points
              contact.update(this.contactManager.contactListener);

              // Was the contact disabled by the user?
              if (!contact.isEnabled()) {
                other.sweep.copy(backup);
                other.synchronizeTransform();
                continue;
              }

              // Are there contact points?
              if (!contact.isTouching()) {
                other.sweep.copy(backup);
                other.synchronizeTransform();
                continue;
              }

              // Add the contact to the island
              contact.islandFlag = true;
              island.addContact(contact);

              // Has the other body already been added to the island?
              if (other.islandFlag) {
                continue;
              }

              // Add the other body to the island.
              other.islandFlag = true;

              if (other.type !== BodyType.StaticBody) {
                other.setAwake(true);
              }

              island.addBody(other);
            }
          }
        }

        const subStep: TimeStep = World.solveTOI_s_subStep;
        subStep.dt = (1 - minAlpha) * step.dt;
        subStep.inv_dt = 1 / subStep.dt;
        subStep.dtRatio = 1;
        subStep.positionIterations = 20;
        subStep.velocityIterations = step.velocityIterations;
        // #if ENABLE_PARTICLE
        subStep.particleIterations = step.particleIterations;
        // #endif
        subStep.warmStarting = false;
        island.solveTOI(subStep, bA.islandIndex, bB.islandIndex);

        // Reset island flags and synchronize broad-phase proxies.
        for (let i: number = 0; i < island.bodyCount; ++i) {
          const body: Body = island.bodies[i];
          body.islandFlag = false;

          if (body.type !== BodyType.DynamicBody) {
            continue;
          }

          body.synchronizeFixtures();

          // Invalidate all contact TOIs on this displaced body.
          for (let ce: ContactEdge = body.contactList; ce; ce = ce.next) {
            ce.contact.toiFlag = false;
            ce.contact.islandFlag = false;
          }
        }

        // Commit fixture proxy movements to the broad-phase so that new contacts are created.
        // Also, some contacts can be destroyed.
        this.contactManager.findNewContacts();

        if (this.subStepping) {
          this.stepComplete = false;
          break;
        }
      }
    }

    // #if ENABLE_CONTROLLER
    public addController(controller: Controller): Controller {
      // Assert(controller.world === null, "Controller can only be a member of one world");
      // controller.world = this;
      controller.next = this.controllerList;
      controller.prev = null;
      if (this.controllerList) {
        this.controllerList.prev = controller;
      }
      this.controllerList = controller;
      ++this.controllerCount;
      return controller;
    }

    public removeController(controller: Controller): Controller {
      // Assert(controller.world === this, "Controller is not a member of this world");
      if (controller.prev) {
        controller.prev.next = controller.next;
      }
      if (controller.next) {
        controller.next.prev = controller.prev;
      }
      if (this.controllerList === controller) {
        this.controllerList = controller.next;
      }
      --this.controllerCount;
      controller.prev = null;
      controller.next = null;
      // delete controller.world; // = null;
      return controller;
    }
    // #endif
  }

}
