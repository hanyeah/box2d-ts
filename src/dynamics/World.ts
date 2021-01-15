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

    public bodyList: Body | null = null;
    public jointList: Joint | null = null;

    // #if ENABLE_PARTICLE
    public particleSystemList: ParticleSystem | null = null;
    // #endif

    public bodyCount: number = 0;
    public jointCount: number = 0;

    public readonly gravity: Vec2 = new Vec2();
    public allowSleep: boolean = true;

    public destructionListener: DestructionListener;
    public debugDraw: Draw;

    // This is used to compute the time step ratio to
    // support a variable time step.
    public inv_dt0: number = 0;

    public newContacts: boolean = false;
    public locked: boolean = false;
    public clearForces: boolean = true;

    // These are for debugging the solver.
    public warmStarting: boolean = true;
    public continuousPhysics: boolean = true;
    public subStepping: boolean = false;

    public stepComplete: boolean = true;

    public readonly profile: Profile = new Profile();

    public readonly island: Island = new Island();

    public readonly s_stack: Array<Body | null> = [];

    // #if ENABLE_CONTROLLER
    public controllerList: Controller | null = null;
    public controllerCount: number = 0;
    // #endif

    /// Construct a world object.
    /// @param gravity the world gravity vector.
    constructor(gravity: XY) {
      this.gravity.Copy(gravity);
    }

    /// Register a destruction listener. The listener is owned by you and must
    /// remain in scope.
    public SetDestructionListener(listener: DestructionListener | null): void {
      this.destructionListener = listener;
    }

    /// Register a contact filter to provide specific control over collision.
    /// Otherwise the default filter is used (defaultFilter). The listener is
    /// owned by you and must remain in scope.
    public SetContactFilter(filter: ContactFilter): void {
      this.contactManager.contactFilter = filter;
    }

    /// Register a contact event listener. The listener is owned by you and must
    /// remain in scope.
    public SetContactListener(listener: ContactListener): void {
      this.contactManager.contactListener = listener;
    }

    /// Register a routine for debug drawing. The debug draw functions are called
    /// inside with World::DebugDraw method. The debug draw object is owned
    /// by you and must remain in scope.
    public SetDebugDraw(debugDraw: Draw): void {
      this.debugDraw = debugDraw;
    }

    /// Create a rigid body given a definition. No reference to the definition
    /// is retained.
    /// @warning This function is locked during callbacks.
    public CreateBody(def: IBodyDef = {}): Body {
      if (this.IsLocked()) { throw new Error(); }

      const b: Body = new Body(def, this);

      // Add to world doubly linked list.
      b.prev = null;
      b.next = this.bodyList;
      if (this.bodyList) {
        this.bodyList.prev = b;
      }
      this.bodyList = b;
      ++this.bodyCount;

      return b;
    }

    /// Destroy a rigid body given a definition. No reference to the definition
    /// is retained. This function is locked during callbacks.
    /// @warning This automatically deletes all associated shapes and joints.
    /// @warning This function is locked during callbacks.
    public DestroyBody(b: Body): void {
      // DEBUG: Assert(this.bodyCount > 0);
      if (this.IsLocked()) { throw new Error(); }

      // Delete the attached joints.
      let je: JointEdge | null = b.jointList;
      while (je) {
        const je0: JointEdge = je;
        je = je.next;

        if (this.destructionListener) {
          this.destructionListener.SayGoodbyeJoint(je0.joint);
        }

        this.DestroyJoint(je0.joint);

        b.jointList = je;
      }
      b.jointList = null;

      // #if ENABLE_CONTROLLER
      // @see Controller list
      let coe: ControllerEdge | null = b.controllerList;
      while (coe) {
        const coe0: ControllerEdge = coe;
        coe = coe.nextController;
        coe0.controller.RemoveBody(b);
      }
      // #endif

      // Delete the attached contacts.
      let ce: ContactEdge | null = b.contactList;
      while (ce) {
        const ce0: ContactEdge = ce;
        ce = ce.next;
        this.contactManager.Destroy(ce0.contact);
      }
      b.contactList = null;

      // Delete the attached fixtures. This destroys broad-phase proxies.
      let f: Fixture | null = b.fixtureList;
      while (f) {
        const f0: Fixture = f;
        f = f.next;

        if (this.destructionListener) {
          this.destructionListener.SayGoodbyeFixture(f0);
        }

        f0.DestroyProxies();
        f0.Reset();

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

    private static _Joint_Create(def: IJointDef): Joint {
      switch (def.type) {
        case JointType.e_distanceJoint: return new DistanceJoint(def as IDistanceJointDef);
        case JointType.e_mouseJoint: return new MouseJoint(def as IMouseJointDef);
        case JointType.e_prismaticJoint: return new PrismaticJoint(def as IPrismaticJointDef);
        case JointType.e_revoluteJoint: return new RevoluteJoint(def as IRevoluteJointDef);
        case JointType.e_pulleyJoint: return new PulleyJoint(def as IPulleyJointDef);
        case JointType.e_gearJoint: return new GearJoint(def as IGearJointDef);
        case JointType.e_wheelJoint: return new WheelJoint(def as IWheelJointDef);
        case JointType.e_weldJoint: return new WeldJoint(def as IWeldJointDef);
        case JointType.e_frictionJoint: return new FrictionJoint(def as IFrictionJointDef);
        case JointType.e_motorJoint: return new MotorJoint(def as IMotorJointDef);
        case JointType.e_areaJoint: return new AreaJoint(def as IAreaJointDef);
      }
      throw new Error();
    }

    private static _Joint_Destroy(joint: Joint): void {
    }

    /// Create a joint to constrain bodies together. No reference to the definition
    /// is retained. This may cause the connected bodies to cease colliding.
    /// @warning This function is locked during callbacks.
    public CreateJoint(def: IAreaJointDef): AreaJoint;
    public CreateJoint(def: IDistanceJointDef): DistanceJoint;
    public CreateJoint(def: IFrictionJointDef): FrictionJoint;
    public CreateJoint(def: IGearJointDef): GearJoint;
    public CreateJoint(def: IMotorJointDef): MotorJoint;
    public CreateJoint(def: IMouseJointDef): MouseJoint;
    public CreateJoint(def: IPrismaticJointDef): PrismaticJoint;
    public CreateJoint(def: IPulleyJointDef): PulleyJoint;
    public CreateJoint(def: IRevoluteJointDef): RevoluteJoint;
    public CreateJoint(def: IWeldJointDef): WeldJoint;
    public CreateJoint(def: IWheelJointDef): WheelJoint;
    public CreateJoint(def: IJointDef): Joint {
      if (this.IsLocked()) { throw new Error(); }

      const j: Joint = World._Joint_Create(def);

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
        let edge: ContactEdge | null = bodyB.GetContactList();
        while (edge) {
          if (edge.other === bodyA) {
            // Flag the contact for filtering at the next time step (where either
            // body is awake).
            edge.contact.FlagForFiltering();
          }

          edge = edge.next;
        }
      }

      // Note: creating a joint doesn't wake the bodies.

      return j;
    }

    /// Destroy a joint. This may cause the connected bodies to begin colliding.
    /// @warning This function is locked during callbacks.
    public DestroyJoint(j: Joint): void {
      if (this.IsLocked()) { throw new Error(); }

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
      bodyA.SetAwake(true);
      bodyB.SetAwake(true);

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

      j.edgeA.Reset();

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

      j.edgeB.Reset();

      World._Joint_Destroy(j);

      // DEBUG: Assert(this.jointCount > 0);
      --this.jointCount;

      // If the joint prevents collisions, then flag any contacts for filtering.
      if (!collideConnected) {
        let edge: ContactEdge | null = bodyB.GetContactList();
        while (edge) {
          if (edge.other === bodyA) {
            // Flag the contact for filtering at the next time step (where either
            // body is awake).
            edge.contact.FlagForFiltering();
          }

          edge = edge.next;
        }
      }
    }

    // #if ENABLE_PARTICLE

    public CreateParticleSystem(def: ParticleSystemDef): ParticleSystem {
      if (this.IsLocked()) { throw new Error(); }

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

    public DestroyParticleSystem(p: ParticleSystem): void {
      if (this.IsLocked()) { throw new Error(); }

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

    public CalculateReasonableParticleIterations(timeStep: number): number {
      if (this.particleSystemList === null) {
        return 1;
      }

      function GetSmallestRadius(world: World): number {
        let smallestRadius = maxFloat;
        for (let system = world.GetParticleSystemList(); system !== null; system = system.next) {
          smallestRadius = Min(smallestRadius, system.GetRadius());
        }
        return smallestRadius;
      }

      // Use the smallest radius, since that represents the worst-case.
      return CalculateParticleIterations(this.gravity.Length(), GetSmallestRadius(this), timeStep);
    }

    // #endif

    /// Take a time step. This performs collision detection, integration,
    /// and constraint solution.
    /// @param timeStep the amount of time to simulate, this should not vary.
    /// @param velocityIterations for the velocity constraint solver.
    /// @param positionIterations for the position constraint solver.
    private static Step_s_step = new TimeStep();
    private static Step_s_stepTimer = new Timer();
    private static Step_s_timer = new Timer();
    // #if ENABLE_PARTICLE
    public Step(dt: number, velocityIterations: number, positionIterations: number, particleIterations: number = this.CalculateReasonableParticleIterations(dt)): void {
      // #else
      // public Step(dt: number, velocityIterations: number, positionIterations: number): void {
      // #endif
      const stepTimer: Timer = World.Step_s_stepTimer.Reset();

      // If new fixtures were added, we need to find the new contacts.
      if (this.newContacts) {
        this.contactManager.FindNewContacts();
        this.newContacts = false;
      }

      this.locked = true;

      const step: TimeStep = World.Step_s_step;
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
      const timer: Timer = World.Step_s_timer.Reset();
      this.contactManager.Collide();
      this.profile.collide = timer.GetMilliseconds();

      // Integrate velocities, solve velocity constraints, and integrate positions.
      if (this.stepComplete && step.dt > 0) {
        const timer: Timer = World.Step_s_timer.Reset();
        // #if ENABLE_PARTICLE
        for (let p = this.particleSystemList; p; p = p.next) {
          p.Solve(step); // Particle Simulation
        }
        // #endif
        this.Solve(step);
        this.profile.solve = timer.GetMilliseconds();
      }

      // Handle TOI events.
      if (this.continuousPhysics && step.dt > 0) {
        const timer: Timer = World.Step_s_timer.Reset();
        this.SolveTOI(step);
        this.profile.solveTOI = timer.GetMilliseconds();
      }

      if (step.dt > 0) {
        this.inv_dt0 = step.inv_dt;
      }

      if (this.clearForces) {
        this.ClearForces();
      }

      this.locked = false;

      this.profile.step = stepTimer.GetMilliseconds();
    }

    /// Manually clear the force buffer on all bodies. By default, forces are cleared automatically
    /// after each call to Step. The default behavior is modified by calling SetAutoClearForces.
    /// The purpose of this function is to support sub-stepping. Sub-stepping is often used to maintain
    /// a fixed sized time step under a variable frame-rate.
    /// When you perform sub-stepping you will disable auto clearing of forces and instead call
    /// ClearForces after all sub-steps are complete in one pass of your game loop.
    /// @see SetAutoClearForces
    public ClearForces(): void {
      for (let body = this.bodyList; body; body = body.next) {
        body.force.SetZero();
        body.torque = 0;
      }
    }

    // #if ENABLE_PARTICLE

    public DrawParticleSystem(system: ParticleSystem): void {
      if (this.debugDraw === null) {
        return;
      }
      const particleCount = system.GetParticleCount();
      if (particleCount) {
        const radius = system.GetRadius();
        const positionBuffer = system.GetPositionBuffer();
        if (system.colorBuffer.data) {
          const colorBuffer = system.GetColorBuffer();
          this.debugDraw.DrawParticles(positionBuffer, radius, colorBuffer, particleCount);
        } else {
          this.debugDraw.DrawParticles(positionBuffer, radius, null, particleCount);
        }
      }
    }

    // #endif

    /// Call this to draw shapes and other debug draw data.
    private static DebugDraw_s_color = new Color(0, 0, 0);
    private static DebugDraw_s_vs = Vec2.MakeArray(4);
    private static DebugDraw_s_xf = new Transform();
    public DebugDraw(): void {
      if (this.debugDraw === null) {
        return;
      }

      const flags: number = this.debugDraw.GetFlags();
      const color: Color = World.DebugDraw_s_color.SetRGB(0, 0, 0);

      if (flags & DrawFlags.e_shapeBit) {
        for (let b: Body | null = this.bodyList; b; b = b.next) {
          const xf: Transform = b.xf;

          this.debugDraw.PushTransform(xf);

          for (let f: Fixture | null = b.GetFixtureList(); f; f = f.next) {
            if (b.GetType() === BodyType.dynamicBody && b.mass === 0.0) {
              // Bad body
              this.DrawShape(f, new Color(1.0, 0.0, 0.0));
            } else if (!b.IsEnabled()) {
              color.SetRGB(0.5, 0.5, 0.3);
              this.DrawShape(f, color);
            } else if (b.GetType() === BodyType.staticBody) {
              color.SetRGB(0.5, 0.9, 0.5);
              this.DrawShape(f, color);
            } else if (b.GetType() === BodyType.kinematicBody) {
              color.SetRGB(0.5, 0.5, 0.9);
              this.DrawShape(f, color);
            } else if (!b.IsAwake()) {
              color.SetRGB(0.6, 0.6, 0.6);
              this.DrawShape(f, color);
            } else {
              color.SetRGB(0.9, 0.7, 0.7);
              this.DrawShape(f, color);
            }
          }

          this.debugDraw.PopTransform(xf);
        }
      }

      // #if ENABLE_PARTICLE
      if (flags & DrawFlags.e_particleBit) {
        for (let p = this.particleSystemList; p; p = p.next) {
          this.DrawParticleSystem(p);
        }
      }
      // #endif

      if (flags & DrawFlags.e_jointBit) {
        for (let j: Joint | null = this.jointList; j; j = j.next) {
          j.Draw(this.debugDraw);
        }
      }

      if (flags & DrawFlags.e_pairBit) {
        color.SetRGB(0.3, 0.9, 0.9);
        for (let contact = this.contactManager.contactList; contact; contact = contact.next) {
          const fixtureA = contact.GetFixtureA();
          const fixtureB = contact.GetFixtureB();
          const indexA = contact.GetChildIndexA();
          const indexB = contact.GetChildIndexB();
          const cA = fixtureA.GetAABB(indexA).GetCenter();
          const cB = fixtureB.GetAABB(indexB).GetCenter();

          this.debugDraw.DrawSegment(cA, cB, color);
        }
      }

      if (flags & DrawFlags.e_aabbBit) {
        color.SetRGB(0.9, 0.3, 0.9);
        const vs: Vec2[] = World.DebugDraw_s_vs;

        for (let b: Body | null = this.bodyList; b; b = b.next) {
          if (!b.IsEnabled()) {
            continue;
          }

          for (let f: Fixture | null = b.GetFixtureList(); f; f = f.next) {
            for (let i: number = 0; i < f.proxyCount; ++i) {
              const proxy: FixtureProxy = f.proxies[i];

              const aabb: AABB = proxy.treeNode.aabb;
              vs[0].Set(aabb.lowerBound.x, aabb.lowerBound.y);
              vs[1].Set(aabb.upperBound.x, aabb.lowerBound.y);
              vs[2].Set(aabb.upperBound.x, aabb.upperBound.y);
              vs[3].Set(aabb.lowerBound.x, aabb.upperBound.y);

              this.debugDraw.DrawPolygon(vs, 4, color);
            }
          }
        }
      }

      if (flags & DrawFlags.e_centerOfMassBit) {
        for (let b: Body | null = this.bodyList; b; b = b.next) {
          const xf: Transform = World.DebugDraw_s_xf;
          xf.q.Copy(b.xf.q);
          xf.p.Copy(b.GetWorldCenter());
          this.debugDraw.DrawTransform(xf);
        }
      }

      // #if ENABLE_CONTROLLER
      // @see Controller list
      if (flags & DrawFlags.e_controllerBit) {
        for (let c = this.controllerList; c; c = c.next) {
          c.Draw(this.debugDraw);
        }
      }
      // #endif
    }

    /// Query the world for all fixtures that potentially overlap the
    /// provided AABB.
    /// @param callback a user implemented callback class.
    /// @param aabb the query box.
    public QueryAABB(callback: QueryCallback, aabb: AABB): void;
    public QueryAABB(aabb: AABB, fn: QueryCallbackFunction): void;
    public QueryAABB(...args: any[]): void {
      if (args[0] instanceof QueryCallback) {
        this._QueryAABB(args[0], args[1]);
      } else {
        this._QueryAABB(null, args[0], args[1]);
      }
    }
    private _QueryAABB(callback: QueryCallback | null, aabb: AABB, fn?: QueryCallbackFunction): void {
      this.contactManager.broadPhase.Query(aabb, (proxy: TreeNode<FixtureProxy>): boolean => {
        const fixture_proxy: FixtureProxy = proxy.userData;
        // DEBUG: Assert(fixture_proxy instanceof FixtureProxy);
        const fixture: Fixture = fixture_proxy.fixture;
        if (callback) {
          return callback.ReportFixture(fixture);
        } else if (fn) {
          return fn(fixture);
        }
        return true;
      });
      // #if ENABLE_PARTICLE
      if (callback instanceof QueryCallback) {
        for (let p = this.particleSystemList; p; p = p.next) {
          if (callback.ShouldQueryParticleSystem(p)) {
            p.QueryAABB(callback, aabb);
          }
        }
      }
      // #endif
    }

    public QueryAllAABB(aabb: AABB, out: Fixture[] = []): Fixture[] {
      this.QueryAABB(aabb, (fixture: Fixture): boolean => { out.push(fixture); return true; });
      return out;
    }

    /// Query the world for all fixtures that potentially overlap the
    /// provided point.
    /// @param callback a user implemented callback class.
    /// @param point the query point.
    public QueryPointAABB(callback: QueryCallback, point: XY): void;
    public QueryPointAABB(point: XY, fn: QueryCallbackFunction): void;
    public QueryPointAABB(...args: any[]): void {
      if (args[0] instanceof QueryCallback) {
        this._QueryPointAABB(args[0], args[1]);
      } else {
        this._QueryPointAABB(null, args[0], args[1]);
      }
    }
    private _QueryPointAABB(callback: QueryCallback | null, point: XY, fn?: QueryCallbackFunction): void {
      this.contactManager.broadPhase.QueryPoint(point, (proxy: TreeNode<FixtureProxy>): boolean => {
        const fixture_proxy: FixtureProxy = proxy.userData;
        // DEBUG: Assert(fixture_proxy instanceof FixtureProxy);
        const fixture: Fixture = fixture_proxy.fixture;
        if (callback) {
          return callback.ReportFixture(fixture);
        } else if (fn) {
          return fn(fixture);
        }
        return true;
      });
      // #if ENABLE_PARTICLE
      if (callback instanceof QueryCallback) {
        for (let p = this.particleSystemList; p; p = p.next) {
          if (callback.ShouldQueryParticleSystem(p)) {
            p.QueryPointAABB(callback, point);
          }
        }
      }
      // #endif
    }

    public QueryAllPointAABB(point: XY, out: Fixture[] = []): Fixture[] {
      this.QueryPointAABB(point, (fixture: Fixture): boolean => { out.push(fixture); return true; });
      return out;
    }

    public QueryFixtureShape(callback: QueryCallback, shape: Shape, index: number, transform: Transform): void;
    public QueryFixtureShape(shape: Shape, index: number, transform: Transform, fn: QueryCallbackFunction): void;
    public QueryFixtureShape(...args: any[]): void {
      if (args[0] instanceof QueryCallback) {
        this._QueryFixtureShape(args[0], args[1], args[2], args[3]);
      } else {
        this._QueryFixtureShape(null, args[0], args[1], args[2], args[3]);
      }
    }
    private static QueryFixtureShape_s_aabb = new AABB();
    private _QueryFixtureShape(callback: QueryCallback | null, shape: Shape, index: number, transform: Transform, fn?: QueryCallbackFunction): void {
      const aabb: AABB = World.QueryFixtureShape_s_aabb;
      shape.ComputeAABB(aabb, transform, index);
      this.contactManager.broadPhase.Query(aabb, (proxy: TreeNode<FixtureProxy>): boolean => {
        const fixture_proxy: FixtureProxy = proxy.userData;
        // DEBUG: Assert(fixture_proxy instanceof FixtureProxy);
        const fixture: Fixture = fixture_proxy.fixture;
        if (TestOverlapShape(shape, index, fixture.GetShape(), fixture_proxy.childIndex, transform, fixture.GetBody().GetTransform())) {
          if (callback) {
            return callback.ReportFixture(fixture);
          } else if (fn) {
            return fn(fixture);
          }
        }
        return true;
      });
      // #if ENABLE_PARTICLE
      if (callback instanceof QueryCallback) {
        for (let p = this.particleSystemList; p; p = p.next) {
          if (callback.ShouldQueryParticleSystem(p)) {
            p.QueryAABB(callback, aabb);
          }
        }
      }
      // #endif
    }

    public QueryAllFixtureShape(shape: Shape, index: number, transform: Transform, out: Fixture[] = []): Fixture[] {
      this.QueryFixtureShape(shape, index, transform, (fixture: Fixture): boolean => { out.push(fixture); return true; });
      return out;
    }

    public QueryFixturePoint(callback: QueryCallback, point: XY): void;
    public QueryFixturePoint(point: XY, fn: QueryCallbackFunction): void;
    public QueryFixturePoint(...args: any[]): void {
      if (args[0] instanceof QueryCallback) {
        this._QueryFixturePoint(args[0], args[1]);
      } else {
        this._QueryFixturePoint(null, args[0], args[1]);
      }
    }
    private _QueryFixturePoint(callback: QueryCallback | null, point: XY, fn?: QueryCallbackFunction): void {
      this.contactManager.broadPhase.QueryPoint(point, (proxy: TreeNode<FixtureProxy>): boolean => {
        const fixture_proxy: FixtureProxy = proxy.userData;
        // DEBUG: Assert(fixture_proxy instanceof FixtureProxy);
        const fixture: Fixture = fixture_proxy.fixture;
        if (fixture.TestPoint(point)) {
          if (callback) {
            return callback.ReportFixture(fixture);
          } else if (fn) {
            return fn(fixture);
          }
        }
        return true;
      });
      // #if ENABLE_PARTICLE
      if (callback) {
        for (let p = this.particleSystemList; p; p = p.next) {
          if (callback.ShouldQueryParticleSystem(p)) {
            p.QueryPointAABB(callback, point);
          }
        }
      }
      // #endif
    }

    public QueryAllFixturePoint(point: XY, out: Fixture[] = []): Fixture[] {
      this.QueryFixturePoint(point, (fixture: Fixture): boolean => { out.push(fixture); return true; });
      return out;
    }

    /// Ray-cast the world for all fixtures in the path of the ray. Your callback
    /// controls whether you get the closest point, any point, or n-points.
    /// The ray-cast ignores shapes that contain the starting point.
    /// @param callback a user implemented callback class.
    /// @param point1 the ray starting point
    /// @param point2 the ray ending point
    public RayCast(callback: RayCastCallback, point1: XY, point2: XY): void;
    public RayCast(point1: XY, point2: XY, fn: RayCastCallbackFunction): void;
    public RayCast(...args: any[]): void {
      if (args[0] instanceof RayCastCallback) {
        this._RayCast(args[0], args[1], args[2]);
      } else {
        this._RayCast(null, args[0], args[1], args[2]);
      }
    }
    private static RayCast_s_input = new RayCastInput();
    private static RayCast_s_output = new RayCastOutput();
    private static RayCast_s_point = new Vec2();
    private _RayCast(callback: RayCastCallback | null, point1: XY, point2: XY, fn?: RayCastCallbackFunction): void {
      const input: RayCastInput = World.RayCast_s_input;
      input.maxFraction = 1;
      input.p1.Copy(point1);
      input.p2.Copy(point2);
      this.contactManager.broadPhase.RayCast(input, (input: RayCastInput, proxy: TreeNode<FixtureProxy>): number => {
        const fixture_proxy: FixtureProxy = proxy.userData;
        // DEBUG: Assert(fixture_proxy instanceof FixtureProxy);
        const fixture: Fixture = fixture_proxy.fixture;
        const index: number = fixture_proxy.childIndex;
        const output: RayCastOutput = World.RayCast_s_output;
        const hit: boolean = fixture.RayCast(output, input, index);
        if (hit) {
          const fraction: number = output.fraction;
          const point: Vec2 = World.RayCast_s_point;
          point.Set((1 - fraction) * point1.x + fraction * point2.x, (1 - fraction) * point1.y + fraction * point2.y);
          if (callback) {
            return callback.ReportFixture(fixture, point, output.normal, fraction);
          } else if (fn) {
            return fn(fixture, point, output.normal, fraction);
          }
        }
        return input.maxFraction;
      });
      // #if ENABLE_PARTICLE
      if (callback) {
        for (let p = this.particleSystemList; p; p = p.next) {
          if (callback.ShouldQueryParticleSystem(p)) {
            p.RayCast(callback, point1, point2);
          }
        }
      }
      // #endif
    }

    public RayCastOne(point1: XY, point2: XY): Fixture | null {
      let result: Fixture | null = null;
      let min_fraction: number = 1;
      this.RayCast(point1, point2, (fixture: Fixture, point: Vec2, normal: Vec2, fraction: number): number => {
        if (fraction < min_fraction) {
          min_fraction = fraction;
          result = fixture;
        }
        return min_fraction;
      });
      return result;
    }

    public RayCastAll(point1: XY, point2: XY, out: Fixture[] = []): Fixture[] {
      this.RayCast(point1, point2, (fixture: Fixture, point: Vec2, normal: Vec2, fraction: number): number => {
        out.push(fixture);
        return 1;
      });
      return out;
    }

    /// Get the world body list. With the returned body, use Body::GetNext to get
    /// the next body in the world list. A NULL body indicates the end of the list.
    /// @return the head of the world body list.
    public GetBodyList(): Body | null {
      return this.bodyList;
    }

    /// Get the world joint list. With the returned joint, use Joint::GetNext to get
    /// the next joint in the world list. A NULL joint indicates the end of the list.
    /// @return the head of the world joint list.
    public GetJointList(): Joint | null {
      return this.jointList;
    }

    // #if ENABLE_PARTICLE
    public GetParticleSystemList(): ParticleSystem | null {
      return this.particleSystemList;
    }
    // #endif

    /// Get the world contact list. With the returned contact, use Contact::GetNext to get
    /// the next contact in the world list. A NULL contact indicates the end of the list.
    /// @return the head of the world contact list.
    /// @warning contacts are created and destroyed in the middle of a time step.
    /// Use ContactListener to avoid missing contacts.
    public GetContactList(): Contact | null {
      return this.contactManager.contactList;
    }

    /// Enable/disable sleep.
    public SetAllowSleeping(flag: boolean): void {
      if (flag === this.allowSleep) {
        return;
      }

      this.allowSleep = flag;
      if (!this.allowSleep) {
        for (let b = this.bodyList; b; b = b.next) {
          b.SetAwake(true);
        }
      }
    }

    public GetAllowSleeping(): boolean {
      return this.allowSleep;
    }

    /// Enable/disable warm starting. For testing.
    public SetWarmStarting(flag: boolean): void {
      this.warmStarting = flag;
    }

    public GetWarmStarting(): boolean {
      return this.warmStarting;
    }

    /// Enable/disable continuous physics. For testing.
    public SetContinuousPhysics(flag: boolean): void {
      this.continuousPhysics = flag;
    }

    public GetContinuousPhysics(): boolean {
      return this.continuousPhysics;
    }

    /// Enable/disable single stepped continuous physics. For testing.
    public SetSubStepping(flag: boolean): void {
      this.subStepping = flag;
    }

    public GetSubStepping(): boolean {
      return this.subStepping;
    }

    /// Get the number of broad-phase proxies.
    public GetProxyCount(): number {
      return this.contactManager.broadPhase.GetProxyCount();
    }

    /// Get the number of bodies.
    public GetBodyCount(): number {
      return this.bodyCount;
    }

    /// Get the number of joints.
    public GetJointCount(): number {
      return this.jointCount;
    }

    /// Get the number of contacts (each may have 0 or more contact points).
    public GetContactCount(): number {
      return this.contactManager.contactCount;
    }

    /// Get the height of the dynamic tree.
    public GetTreeHeight(): number {
      return this.contactManager.broadPhase.GetTreeHeight();
    }

    /// Get the balance of the dynamic tree.
    public GetTreeBalance(): number {
      return this.contactManager.broadPhase.GetTreeBalance();
    }

    /// Get the quality metric of the dynamic tree. The smaller the better.
    /// The minimum is 1.
    public GetTreeQuality(): number {
      return this.contactManager.broadPhase.GetTreeQuality();
    }

    /// Change the global gravity vector.
    public SetGravity(gravity: XY, wake: boolean = true) {
      if (!Vec2.IsEqualToV(this.gravity, gravity)) {
        this.gravity.Copy(gravity);

        if (wake) {
          for (let b: Body | null = this.bodyList; b; b = b.next) {
            b.SetAwake(true);
          }
        }
      }
    }

    /// Get the global gravity vector.
    public GetGravity(): Vec2 {
      return this.gravity;
    }

    /// Is the world locked (in the middle of a time step).
    public IsLocked(): boolean {
      return this.locked;
    }

    /// Set flag to control automatic clearing of forces after each time step.
    public SetAutoClearForces(flag: boolean): void {
      this.clearForces = flag;
    }

    /// Get the flag that controls automatic clearing of forces after each time step.
    public GetAutoClearForces(): boolean {
      return this.clearForces;
    }

    /// Shift the world origin. Useful for large worlds.
    /// The body shift formula is: position -= newOrigin
    /// @param newOrigin the new origin with respect to the old origin
    public ShiftOrigin(newOrigin: XY): void {
      if (this.IsLocked()) { throw new Error(); }

      for (let b: Body | null = this.bodyList; b; b = b.next) {
        b.xf.p.SelfSub(newOrigin);
        b.sweep.c0.SelfSub(newOrigin);
        b.sweep.c.SelfSub(newOrigin);
      }

      for (let j: Joint | null = this.jointList; j; j = j.next) {
        j.ShiftOrigin(newOrigin);
      }

      this.contactManager.broadPhase.ShiftOrigin(newOrigin);
    }

    /// Get the contact manager for testing.
    public GetContactManager(): ContactManager {
      return this.contactManager;
    }

    /// Get the current profile.
    public GetProfile(): Profile {
      return this.profile;
    }

    /// Dump the world into the log file.
    /// @warning this should be called outside of a time step.
    public Dump(log: (format: string, ...args: any[]) => void): void {
      if (this.locked) {
        return;
      }

      // OpenDump("box2d_dump.inl");

      log("const g: Vec2 = new Vec2(%.15f, %.15f);\n", this.gravity.x, this.gravity.y);
      log("this.world.SetGravity(g);\n");

      log("const bodies: Body[] = [];\n");
      log("const joints: Joint[] = [];\n");
      let i: number = 0;
      for (let b: Body | null = this.bodyList; b; b = b.next) {
        b.islandIndex = i;
        b.Dump(log);
        ++i;
      }

      i = 0;
      for (let j: Joint | null = this.jointList; j; j = j.next) {
        j.index = i;
        ++i;
      }

      // First pass on joints, skip gear joints.
      for (let j: Joint | null = this.jointList; j; j = j.next) {
        if (j.type === JointType.e_gearJoint) {
          continue;
        }

        log("{\n");
        j.Dump(log);
        log("}\n");
      }

      // Second pass on joints, only gear joints.
      for (let j: Joint | null = this.jointList; j; j = j.next) {
        if (j.type !== JointType.e_gearJoint) {
          continue;
        }

        log("{\n");
        j.Dump(log);
        log("}\n");
      }

      // CloseDump();
    }

    public DrawShape(fixture: Fixture, color: Color): void {
      if (this.debugDraw === null) {
        return;
      }
      const shape: Shape = fixture.GetShape();

      switch (shape.type) {
        case ShapeType.e_circleShape: {
          const circle: CircleShape = shape as CircleShape;
          const center: Vec2 = circle.p;
          const radius: number = circle.radius;
          const axis: Vec2 = Vec2.UNITX;
          this.debugDraw.DrawSolidCircle(center, radius, axis, color);
          break;
        }

        case ShapeType.e_edgeShape: {
          const edge: EdgeShape = shape as EdgeShape;
          const v1: Vec2 = edge.vertex1;
          const v2: Vec2 = edge.vertex2;
          this.debugDraw.DrawSegment(v1, v2, color);

          if (edge.oneSided === false) {
            this.debugDraw.DrawPoint(v1, 4.0, color);
            this.debugDraw.DrawPoint(v2, 4.0, color);
          }
          break;
        }

        case ShapeType.e_chainShape: {
          const chain: ChainShape = shape as ChainShape;
          const count: number = chain.count;
          const vertices: Vec2[] = chain.vertices;
          let v1: Vec2 = vertices[0];
          for (let i: number = 1; i < count; ++i) {
            const v2: Vec2 = vertices[i];
            this.debugDraw.DrawSegment(v1, v2, color);
            v1 = v2;
          }

          break;
        }

        case ShapeType.e_polygonShape: {
          const poly: PolygonShape = shape as PolygonShape;
          const vertexCount: number = poly.count;
          const vertices: Vec2[] = poly.vertices;
          this.debugDraw.DrawSolidPolygon(vertices, vertexCount, color);
          break;
        }
      }
    }

    public Solve(step: TimeStep): void {
      // #if ENABLE_PARTICLE
      // update previous transforms
      for (let b = this.bodyList; b; b = b.next) {
        b.xf0.Copy(b.xf);
      }
      // #endif

      // #if ENABLE_CONTROLLER
      // @see Controller list
      for (let controller = this.controllerList; controller; controller = controller.next) {
        controller.Step(step);
      }
      // #endif

      this.profile.solveInit = 0;
      this.profile.solveVelocity = 0;
      this.profile.solvePosition = 0;

      // Size the island for the worst case.
      const island: Island = this.island;
      island.Initialize(this.bodyCount,
        this.contactManager.contactCount,
        this.jointCount,
        this.contactManager.contactListener);

      // Clear all the island flags.
      for (let b: Body | null = this.bodyList; b; b = b.next) {
        b.islandFlag = false;
      }
      for (let c: Contact | null = this.contactManager.contactList; c; c = c.next) {
        c.islandFlag = false;
      }
      for (let j: Joint | null = this.jointList; j; j = j.next) {
        j.islandFlag = false;
      }

      // Build and simulate all awake islands.
      // DEBUG: const stackSize: number = this.bodyCount;
      const stack: Array<Body | null> = this.s_stack;
      for (let seed: Body | null = this.bodyList; seed; seed = seed.next) {
        if (seed.islandFlag) {
          continue;
        }

        if (!seed.IsAwake() || !seed.IsEnabled()) {
          continue;
        }

        // The seed can be dynamic or kinematic.
        if (seed.GetType() === BodyType.staticBody) {
          continue;
        }

        // Reset island and stack.
        island.Clear();
        let stackCount: number = 0;
        stack[stackCount++] = seed;
        seed.islandFlag = true;

        // Perform a depth first search (DFS) on the constraint graph.
        while (stackCount > 0) {
          // Grab the next body off the stack and add it to the island.
          const b: Body | null = stack[--stackCount];
          if (!b) { throw new Error(); }
          // DEBUG: Assert(b.IsEnabled());
          island.AddBody(b);

          // To keep islands as small as possible, we don't
          // propagate islands across static bodies.
          if (b.GetType() === BodyType.staticBody) {
            continue;
          }

          // Make sure the body is awake. (without resetting sleep timer).
          b.awakeFlag = true;

          // Search all contacts connected to this body.
          for (let ce: ContactEdge | null = b.contactList; ce; ce = ce.next) {
            const contact: Contact = ce.contact;

            // Has this contact already been added to an island?
            if (contact.islandFlag) {
              continue;
            }

            // Is this contact solid and touching?
            if (!contact.IsEnabled() || !contact.IsTouching()) {
              continue;
            }

            // Skip sensors.
            const sensorA: boolean = contact.fixtureA.isSensor;
            const sensorB: boolean = contact.fixtureB.isSensor;
            if (sensorA || sensorB) {
              continue;
            }

            island.AddContact(contact);
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
          for (let je: JointEdge | null = b.jointList; je; je = je.next) {
            if (je.joint.islandFlag) {
              continue;
            }

            const other: Body = je.other;

            // Don't simulate joints connected to disabled bodies.
            if (!other.IsEnabled()) {
              continue;
            }

            island.AddJoint(je.joint);
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
        island.Solve(profile, step, this.gravity, this.allowSleep);
        this.profile.solveInit += profile.solveInit;
        this.profile.solveVelocity += profile.solveVelocity;
        this.profile.solvePosition += profile.solvePosition;

        // Post solve cleanup.
        for (let i: number = 0; i < island.bodyCount; ++i) {
          // Allow static bodies to participate in other islands.
          const b: Body = island.bodies[i];
          if (b.GetType() === BodyType.staticBody) {
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

        if (b.GetType() === BodyType.staticBody) {
          continue;
        }

        // Update fixtures (for broad-phase).
        b.SynchronizeFixtures();
      }

      // Look for new contacts.
      this.contactManager.FindNewContacts();
      this.profile.broadphase = timer.GetMilliseconds();
    }

    private static SolveTOI_s_subStep = new TimeStep();
    private static SolveTOI_s_backup = new Sweep();
    private static SolveTOI_s_backup1 = new Sweep();
    private static SolveTOI_s_backup2 = new Sweep();
    private static SolveTOI_s_toi_input = new TOIInput();
    private static SolveTOI_s_toi_output = new TOIOutput();
    public SolveTOI(step: TimeStep): void {
      const island: Island = this.island;
      island.Initialize(2 * maxTOIContacts, maxTOIContacts, 0, this.contactManager.contactListener);

      if (this.stepComplete) {
        for (let b: Body | null = this.bodyList; b; b = b.next) {
          b.islandFlag = false;
          b.sweep.alpha0 = 0;
        }

        for (let c: Contact | null = this.contactManager.contactList; c; c = c.next) {
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
        let minContact: Contact | null = null;
        let minAlpha: number = 1;

        for (let c: Contact | null = this.contactManager.contactList; c; c = c.next) {
          // Is this contact disabled?
          if (!c.IsEnabled()) {
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
            const fA: Fixture = c.GetFixtureA();
            const fB: Fixture = c.GetFixtureB();

            // Is there a sensor?
            if (fA.IsSensor() || fB.IsSensor()) {
              continue;
            }

            const bA: Body = fA.GetBody();
            const bB: Body = fB.GetBody();

            const typeA: BodyType = bA.type;
            const typeB: BodyType = bB.type;
            // DEBUG: Assert(typeA !== BodyType.staticBody || typeB !== BodyType.staticBody);

            const activeA: boolean = bA.IsAwake() && typeA !== BodyType.staticBody;
            const activeB: boolean = bB.IsAwake() && typeB !== BodyType.staticBody;

            // Is at least one body active (awake and dynamic or kinematic)?
            if (!activeA && !activeB) {
              continue;
            }

            const collideA: boolean = bA.IsBullet() || typeA !== BodyType.dynamicBody;
            const collideB: boolean = bB.IsBullet() || typeB !== BodyType.dynamicBody;

            // Are these two non-bullet dynamic bodies?
            if (!collideA && !collideB) {
              continue;
            }

            // Compute the TOI for this contact.
            // Put the sweeps onto the same time interval.
            let alpha0: number = bA.sweep.alpha0;

            if (bA.sweep.alpha0 < bB.sweep.alpha0) {
              alpha0 = bB.sweep.alpha0;
              bA.sweep.Advance(alpha0);
            } else if (bB.sweep.alpha0 < bA.sweep.alpha0) {
              alpha0 = bA.sweep.alpha0;
              bB.sweep.Advance(alpha0);
            }

            // DEBUG: Assert(alpha0 < 1);

            const indexA: number = c.GetChildIndexA();
            const indexB: number = c.GetChildIndexB();

            // Compute the time of impact in interval [0, minTOI]
            const input: TOIInput = World.SolveTOI_s_toi_input;
            input.proxyA.SetShape(fA.GetShape(), indexA);
            input.proxyB.SetShape(fB.GetShape(), indexB);
            input.sweepA.Copy(bA.sweep);
            input.sweepB.Copy(bB.sweep);
            input.tMax = 1;

            const output: TOIOutput = World.SolveTOI_s_toi_output;
            TimeOfImpact(output, input);

            // Beta is the fraction of the remaining portion of the .
            const beta: number = output.t;
            if (output.state === TOIOutputState.e_touching) {
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
        const fA: Fixture = minContact.GetFixtureA();
        const fB: Fixture = minContact.GetFixtureB();
        const bA: Body = fA.GetBody();
        const bB: Body = fB.GetBody();

        const backup1: Sweep = World.SolveTOI_s_backup1.Copy(bA.sweep);
        const backup2: Sweep = World.SolveTOI_s_backup2.Copy(bB.sweep);

        bA.Advance(minAlpha);
        bB.Advance(minAlpha);

        // The TOI contact likely has some new contact points.
        minContact.Update(this.contactManager.contactListener);
        minContact.toiFlag = false;
        ++minContact.toiCount;

        // Is the contact solid?
        if (!minContact.IsEnabled() || !minContact.IsTouching()) {
          // Restore the sweeps.
          minContact.SetEnabled(false);
          bA.sweep.Copy(backup1);
          bB.sweep.Copy(backup2);
          bA.SynchronizeTransform();
          bB.SynchronizeTransform();
          continue;
        }

        bA.SetAwake(true);
        bB.SetAwake(true);

        // Build the island
        island.Clear();
        island.AddBody(bA);
        island.AddBody(bB);
        island.AddContact(minContact);

        bA.islandFlag = true;
        bB.islandFlag = true;
        minContact.islandFlag = true;

        // Get contacts on bodyA and bodyB.
        // const bodies: Body[] = [bA, bB];
        for (let i: number = 0; i < 2; ++i) {
          const body: Body = (i === 0) ? (bA) : (bB); // bodies[i];
          if (body.type === BodyType.dynamicBody) {
            for (let ce: ContactEdge | null = body.contactList; ce; ce = ce.next) {
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
              if (other.type === BodyType.dynamicBody &&
                !body.IsBullet() && !other.IsBullet()) {
                continue;
              }

              // Skip sensors.
              const sensorA: boolean = contact.fixtureA.isSensor;
              const sensorB: boolean = contact.fixtureB.isSensor;
              if (sensorA || sensorB) {
                continue;
              }

              // Tentatively advance the body to the TOI.
              const backup: Sweep = World.SolveTOI_s_backup.Copy(other.sweep);
              if (!other.islandFlag) {
                other.Advance(minAlpha);
              }

              // Update the contact points
              contact.Update(this.contactManager.contactListener);

              // Was the contact disabled by the user?
              if (!contact.IsEnabled()) {
                other.sweep.Copy(backup);
                other.SynchronizeTransform();
                continue;
              }

              // Are there contact points?
              if (!contact.IsTouching()) {
                other.sweep.Copy(backup);
                other.SynchronizeTransform();
                continue;
              }

              // Add the contact to the island
              contact.islandFlag = true;
              island.AddContact(contact);

              // Has the other body already been added to the island?
              if (other.islandFlag) {
                continue;
              }

              // Add the other body to the island.
              other.islandFlag = true;

              if (other.type !== BodyType.staticBody) {
                other.SetAwake(true);
              }

              island.AddBody(other);
            }
          }
        }

        const subStep: TimeStep = World.SolveTOI_s_subStep;
        subStep.dt = (1 - minAlpha) * step.dt;
        subStep.inv_dt = 1 / subStep.dt;
        subStep.dtRatio = 1;
        subStep.positionIterations = 20;
        subStep.velocityIterations = step.velocityIterations;
        // #if ENABLE_PARTICLE
        subStep.particleIterations = step.particleIterations;
        // #endif
        subStep.warmStarting = false;
        island.SolveTOI(subStep, bA.islandIndex, bB.islandIndex);

        // Reset island flags and synchronize broad-phase proxies.
        for (let i: number = 0; i < island.bodyCount; ++i) {
          const body: Body = island.bodies[i];
          body.islandFlag = false;

          if (body.type !== BodyType.dynamicBody) {
            continue;
          }

          body.SynchronizeFixtures();

          // Invalidate all contact TOIs on this displaced body.
          for (let ce: ContactEdge | null = body.contactList; ce; ce = ce.next) {
            ce.contact.toiFlag = false;
            ce.contact.islandFlag = false;
          }
        }

        // Commit fixture proxy movements to the broad-phase so that new contacts are created.
        // Also, some contacts can be destroyed.
        this.contactManager.FindNewContacts();

        if (this.subStepping) {
          this.stepComplete = false;
          break;
        }
      }
    }

    // #if ENABLE_CONTROLLER
    public AddController(controller: Controller): Controller {
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

    public RemoveController(controller: Controller): Controller {
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
