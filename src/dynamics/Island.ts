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

/*
Position Correction Notes
=========================
I tried the several algorithms for position correction of the 2D revolute joint.
I looked at these systems:
- simple pendulum (1m diameter sphere on massless 5m stick) with initial angular velocity of 100 rad/s.
- suspension bridge with 30 1m long planks of length 1m.
- multi-link chain with 30 1m long links.

Here are the algorithms:

Baumgarte - A fraction of the position error is added to the velocity error. There is no
separate position solver.

Pseudo Velocities - After the velocity solver and position integration,
the position error, Jacobian, and effective mass are recomputed. Then
the velocity constraints are solved with pseudo velocities and a fraction
of the position error is added to the pseudo velocity error. The pseudo
velocities are initialized to zero and there is no warm-starting. After
the position solver, the pseudo velocities are added to the positions.
This is also called the First Order World method or the Position LCP method.

Modified Nonlinear Gauss-Seidel (NGS) - Like Pseudo Velocities except the
position error is re-computed for each constraint and the positions are updated
after the constraint is solved. The radius vectors (aka Jacobians) are
re-computed too (otherwise the algorithm has horrible instability). The pseudo
velocity states are not needed because they are effectively zero at the beginning
of each iteration. Since we have the current position error, we allow the
iterations to terminate early if the error becomes smaller than linearSlop.

Full NGS or just NGS - Like Modified NGS except the effective mass are re-computed
each time a constraint is solved.

Here are the results:
Baumgarte - this is the cheapest algorithm but it has some stability problems,
especially with the bridge. The chain links separate easily close to the root
and they jitter as they struggle to pull together. This is one of the most common
methods in the field. The big drawback is that the position correction artificially
affects the momentum, thus leading to instabilities and false bounce. I used a
bias factor of 0.2. A larger bias factor makes the bridge less stable, a smaller
factor makes joints and contacts more spongy.

Pseudo Velocities - the is more stable than the Baumgarte method. The bridge is
stable. However, joints still separate with large angular velocities. Drag the
simple pendulum in a circle quickly and the joint will separate. The chain separates
easily and does not recover. I used a bias factor of 0.2. A larger value lead to
the bridge collapsing when a heavy cube drops on it.

Modified NGS - this algorithm is better in some ways than Baumgarte and Pseudo
Velocities, but in other ways it is worse. The bridge and chain are much more
stable, but the simple pendulum goes unstable at high angular velocities.

Full NGS - stable in all tests. The joints display good stiffness. The bridge
still sags, but this is better than infinite forces.

Recommendations
Pseudo Velocities are not really worthwhile because the bridge and chain cannot
recover from joint separation. In other cases the benefit over Baumgarte is small.

Modified NGS is not a robust method for the revolute joint due to the violent
instability seen in the simple pendulum. Perhaps it is viable with other constraint
types, especially scalar constraints where the effective mass is a scalar.

This leaves Baumgarte and Full NGS. Baumgarte has small, but manageable instabilities
and is very fast. I don't think we can escape Baumgarte, especially in highly
demanding cases where high constraint fidelity is not needed.

Full NGS is robust and easy on the eyes. I recommend this as an option for
higher fidelity simulation and certainly for suspension bridges and long chains.
Full NGS might be a good choice for ragdolls, especially motorized ragdolls where
joint separation can be problematic. The number of NGS iterations can be reduced
for better performance without harming robustness much.

Each joint in a can be handled differently in the position solver. So I recommend
a system where the user can select the algorithm on a per joint basis. I would
probably default to the slower Full NGS and let the user select the faster
Baumgarte method in performance critical scenarios.
*/

/*
Cache Performance

The Box2D solvers are dominated by cache misses. Data structures are designed
to increase the number of cache hits. Much of misses are due to random access
to body data. The constraint structures are iterated over linearly, which leads
to few cache misses.

The bodies are not accessed during iteration. Instead read only data, such as
the mass values are stored with the constraints. The mutable data are the constraint
impulses and the bodies velocities/positions. The impulses are held inside the
constraint structures. The body velocities/positions are held in compact, temporary
arrays to increase the number of cache hits. Linear and angular velocity are
stored in a single array since multiple arrays lead to multiple misses.
*/

/*
2D Rotation

R = [cos(theta) -sin(theta)]
    [sin(theta) cos(theta) ]

thetaDot = omega

Let q1 = cos(theta), q2 = sin(theta).
R = [q1 -q2]
    [q2  q1]

q1Dot = -thetaDot * q2
q2Dot = thetaDot * q1

q1_new = q1_old - dt * w * q2
q2_new = q2_old + dt * w * q1
then normalize.

This might be faster than computing sin+cos.
However, we can compute sin+cos of the same angle fast.
*/
namespace b2 {
  export class Island {
    public listener!: ContactListener;

    public readonly bodies: Body[] = [/*1024*/]; // TODO: Settings
    public readonly contacts: Contact[] = [/*1024*/]; // TODO: Settings
    public readonly joints: Joint[] = [/*1024*/]; // TODO: Settings

    public readonly positions: Position[] = Position.MakeArray(1024); // TODO: Settings
    public readonly velocities: Velocity[] = Velocity.MakeArray(1024); // TODO: Settings

    public bodyCount: number = 0;
    public jointCount: number = 0;
    public contactCount: number = 0;

    public bodyCapacity: number = 0;
    public contactCapacity: number = 0;
    public jointCapacity: number = 0;

    public Initialize(bodyCapacity: number, contactCapacity: number, jointCapacity: number, listener: ContactListener): void {
      this.bodyCapacity = bodyCapacity;
      this.contactCapacity = contactCapacity;
      this.jointCapacity = jointCapacity;
      this.bodyCount = 0;
      this.contactCount = 0;
      this.jointCount = 0;

      this.listener = listener;

      // TODO:
      // while (this.bodies.length < bodyCapacity) {
      //   this.bodies[this.bodies.length] = null;
      // }
      // TODO:
      // while (this.contacts.length < contactCapacity) {
      //   this.contacts[this.contacts.length] = null;
      // }
      // TODO:
      // while (this.joints.length < jointCapacity) {
      //   this.joints[this.joints.length] = null;
      // }

      // TODO:
      if (this.positions.length < bodyCapacity) {
        const new_length = Max(this.positions.length * 2, bodyCapacity);
        while (this.positions.length < new_length) {
          this.positions[this.positions.length] = new Position();
        }
      }
      // TODO:
      if (this.velocities.length < bodyCapacity) {
        const new_length = Max(this.velocities.length * 2, bodyCapacity);
        while (this.velocities.length < new_length) {
          this.velocities[this.velocities.length] = new Velocity();
        }
      }
    }

    public Clear(): void {
      this.bodyCount = 0;
      this.contactCount = 0;
      this.jointCount = 0;
    }

    public AddBody(body: Body): void {
      // DEBUG: Assert(this.bodyCount < this.bodyCapacity);
      body.islandIndex = this.bodyCount;
      this.bodies[this.bodyCount++] = body;
    }

    public AddContact(contact: Contact): void {
      // DEBUG: Assert(this.contactCount < this.contactCapacity);
      this.contacts[this.contactCount++] = contact;
    }

    public AddJoint(joint: Joint): void {
      // DEBUG: Assert(this.jointCount < this.jointCapacity);
      this.joints[this.jointCount++] = joint;
    }

    private static s_timer = new Timer();
    private static s_solverData = new SolverData();
    private static s_contactSolverDef = new ContactSolverDef();
    private static s_contactSolver = new ContactSolver();
    private static s_translation = new Vec2();
    public Solve(profile: Profile, step: TimeStep, gravity: Vec2, allowSleep: boolean): void {
      const timer: Timer = Island.s_timer.Reset();

      const h: number = step.dt;

      // Integrate velocities and apply damping. Initialize the body state.
      for (let i: number = 0; i < this.bodyCount; ++i) {
        const b: Body = this.bodies[i];

        // const c: Vec2 =
        this.positions[i].c.Copy(b.sweep.c);
        const a: number = b.sweep.a;
        const v: Vec2 = this.velocities[i].v.Copy(b.linearVelocity);
        let w: number = b.angularVelocity;

        // Store positions for continuous collision.
        b.sweep.c0.Copy(b.sweep.c);
        b.sweep.a0 = b.sweep.a;

        if (b.type === BodyType.dynamicBody) {
          // Integrate velocities.
          // v += h * b->invMass * (b->gravityScale * b->mass * gravity + b->force);
          v.x += h * b.invMass * (b.gravityScale * b.mass * gravity.x + b.force.x);
          v.y += h * b.invMass * (b.gravityScale * b.mass * gravity.y + b.force.y);
          w += h * b.invI * b.torque;

          // Apply damping.
          // ODE: dv/dt + c * v = 0
          // Solution: v(t) = v0 * exp(-c * t)
          // Time step: v(t + dt) = v0 * exp(-c * (t + dt)) = v0 * exp(-c * t) * exp(-c * dt) = v * exp(-c * dt)
          // v2 = exp(-c * dt) * v1
          // Pade approximation:
          // v2 = v1 * 1 / (1 + c * dt)
          v.SelfMul(1.0 / (1.0 + h * b.linearDamping));
          w *= 1.0 / (1.0 + h * b.angularDamping);
        }

        // this.positions[i].c = c;
        this.positions[i].a = a;
        // this.velocities[i].v = v;
        this.velocities[i].w = w;
      }

      timer.Reset();

      // Solver data
      const solverData: SolverData = Island.s_solverData;
      solverData.step.Copy(step);
      solverData.positions = this.positions;
      solverData.velocities = this.velocities;

      // Initialize velocity constraints.
      const contactSolverDef: ContactSolverDef = Island.s_contactSolverDef;
      contactSolverDef.step.Copy(step);
      contactSolverDef.contacts = this.contacts;
      contactSolverDef.count = this.contactCount;
      contactSolverDef.positions = this.positions;
      contactSolverDef.velocities = this.velocities;

      const contactSolver: ContactSolver = Island.s_contactSolver.Initialize(contactSolverDef);
      contactSolver.InitializeVelocityConstraints();

      if (step.warmStarting) {
        contactSolver.WarmStart();
      }

      for (let i: number = 0; i < this.jointCount; ++i) {
        this.joints[i].InitVelocityConstraints(solverData);
      }

      profile.solveInit = timer.GetMilliseconds();

      // Solve velocity constraints.
      timer.Reset();
      for (let i: number = 0; i < step.velocityIterations; ++i) {
        for (let j: number = 0; j < this.jointCount; ++j) {
          this.joints[j].SolveVelocityConstraints(solverData);
        }

        contactSolver.SolveVelocityConstraints();
      }

      // Store impulses for warm starting
      contactSolver.StoreImpulses();
      profile.solveVelocity = timer.GetMilliseconds();

      // Integrate positions.
      for (let i: number = 0; i < this.bodyCount; ++i) {
        const c: Vec2 = this.positions[i].c;
        let a: number = this.positions[i].a;
        const v: Vec2 = this.velocities[i].v;
        let w: number = this.velocities[i].w;

        // Check for large velocities
        const translation: Vec2 = Vec2.MulSV(h, v, Island.s_translation);
        if (Vec2.DotVV(translation, translation) > maxTranslationSquared) {
          const ratio: number = maxTranslation / translation.Length();
          v.SelfMul(ratio);
        }

        const rotation: number = h * w;
        if (rotation * rotation > maxRotationSquared) {
          const ratio: number = maxRotation / Abs(rotation);
          w *= ratio;
        }

        // Integrate
        c.x += h * v.x;
        c.y += h * v.y;
        a += h * w;

        // this.positions[i].c = c;
        this.positions[i].a = a;
        // this.velocities[i].v = v;
        this.velocities[i].w = w;
      }

      // Solve position constraints
      timer.Reset();
      let positionSolved: boolean = false;
      for (let i: number = 0; i < step.positionIterations; ++i) {
        const contactsOkay: boolean = contactSolver.SolvePositionConstraints();

        let jointsOkay: boolean = true;
        for (let j: number = 0; j < this.jointCount; ++j) {
          const jointOkay: boolean = this.joints[j].SolvePositionConstraints(solverData);
          jointsOkay = jointsOkay && jointOkay;
        }

        if (contactsOkay && jointsOkay) {
          // Exit early if the position errors are small.
          positionSolved = true;
          break;
        }
      }

      // Copy state buffers back to the bodies
      for (let i: number = 0; i < this.bodyCount; ++i) {
        const body: Body = this.bodies[i];
        body.sweep.c.Copy(this.positions[i].c);
        body.sweep.a = this.positions[i].a;
        body.linearVelocity.Copy(this.velocities[i].v);
        body.angularVelocity = this.velocities[i].w;
        body.SynchronizeTransform();
      }

      profile.solvePosition = timer.GetMilliseconds();

      this.Report(contactSolver.velocityConstraints);

      if (allowSleep) {
        let minSleepTime: number = maxFloat;

        const linTolSqr: number = linearSleepTolerance * linearSleepTolerance;
        const angTolSqr: number = angularSleepTolerance * angularSleepTolerance;

        for (let i: number = 0; i < this.bodyCount; ++i) {
          const b: Body = this.bodies[i];
          if (b.GetType() === BodyType.staticBody) {
            continue;
          }

          if (!b.autoSleepFlag ||
            b.angularVelocity * b.angularVelocity > angTolSqr ||
            Vec2.DotVV(b.linearVelocity, b.linearVelocity) > linTolSqr) {
            b.sleepTime = 0;
            minSleepTime = 0;
          } else {
            b.sleepTime += h;
            minSleepTime = Min(minSleepTime, b.sleepTime);
          }
        }

        if (minSleepTime >= timeToSleep && positionSolved) {
          for (let i: number = 0; i < this.bodyCount; ++i) {
            const b: Body = this.bodies[i];
            b.SetAwake(false);
          }
        }
      }
    }

    public SolveTOI(subStep: TimeStep, toiIndexA: number, toiIndexB: number): void {
      // DEBUG: Assert(toiIndexA < this.bodyCount);
      // DEBUG: Assert(toiIndexB < this.bodyCount);

      // Initialize the body state.
      for (let i: number = 0; i < this.bodyCount; ++i) {
        const b: Body = this.bodies[i];
        this.positions[i].c.Copy(b.sweep.c);
        this.positions[i].a = b.sweep.a;
        this.velocities[i].v.Copy(b.linearVelocity);
        this.velocities[i].w = b.angularVelocity;
      }

      const contactSolverDef: ContactSolverDef = Island.s_contactSolverDef;
      contactSolverDef.contacts = this.contacts;
      contactSolverDef.count = this.contactCount;
      contactSolverDef.step.Copy(subStep);
      contactSolverDef.positions = this.positions;
      contactSolverDef.velocities = this.velocities;
      const contactSolver: ContactSolver = Island.s_contactSolver.Initialize(contactSolverDef);

      // Solve position constraints.
      for (let i: number = 0; i < subStep.positionIterations; ++i) {
        const contactsOkay: boolean = contactSolver.SolveTOIPositionConstraints(toiIndexA, toiIndexB);
        if (contactsOkay) {
          break;
        }
      }

      /*
      #if 0
        // Is the new position really safe?
        for (int32 i = 0; i < this.contactCount; ++i) {
          Contact* c = this.contacts[i];
          Fixture* fA = c.GetFixtureA();
          Fixture* fB = c.GetFixtureB();

          Body* bA = fA.GetBody();
          Body* bB = fB.GetBody();

          int32 indexA = c.GetChildIndexA();
          int32 indexB = c.GetChildIndexB();

          DistanceInput input;
          input.proxyA.Set(fA.GetShape(), indexA);
          input.proxyB.Set(fB.GetShape(), indexB);
          input.transformA = bA.GetTransform();
          input.transformB = bB.GetTransform();
          input.useRadii = false;

          DistanceOutput output;
          SimplexCache cache;
          cache.count = 0;
          Distance(&output, &cache, &input);

          if (output.distance === 0 || cache.count === 3) {
            cache.count += 0;
          }
        }
      #endif
      */

      // Leap of faith to new safe state.
      this.bodies[toiIndexA].sweep.c0.Copy(this.positions[toiIndexA].c);
      this.bodies[toiIndexA].sweep.a0 = this.positions[toiIndexA].a;
      this.bodies[toiIndexB].sweep.c0.Copy(this.positions[toiIndexB].c);
      this.bodies[toiIndexB].sweep.a0 = this.positions[toiIndexB].a;

      // No warm starting is needed for TOI events because warm
      // starting impulses were applied in the discrete solver.
      contactSolver.InitializeVelocityConstraints();

      // Solve velocity constraints.
      for (let i: number = 0; i < subStep.velocityIterations; ++i) {
        contactSolver.SolveVelocityConstraints();
      }

      // Don't store the TOI contact forces for warm starting
      // because they can be quite large.

      const h: number = subStep.dt;

      // Integrate positions
      for (let i: number = 0; i < this.bodyCount; ++i) {
        const c: Vec2 = this.positions[i].c;
        let a: number = this.positions[i].a;
        const v: Vec2 = this.velocities[i].v;
        let w: number = this.velocities[i].w;

        // Check for large velocities
        const translation: Vec2 = Vec2.MulSV(h, v, Island.s_translation);
        if (Vec2.DotVV(translation, translation) > maxTranslationSquared) {
          const ratio: number = maxTranslation / translation.Length();
          v.SelfMul(ratio);
        }

        const rotation: number = h * w;
        if (rotation * rotation > maxRotationSquared) {
          const ratio: number = maxRotation / Abs(rotation);
          w *= ratio;
        }

        // Integrate
        c.SelfMulAdd(h, v);
        a += h * w;

        // this.positions[i].c = c;
        this.positions[i].a = a;
        // this.velocities[i].v = v;
        this.velocities[i].w = w;

        // Sync bodies
        const body: Body = this.bodies[i];
        body.sweep.c.Copy(c);
        body.sweep.a = a;
        body.linearVelocity.Copy(v);
        body.angularVelocity = w;
        body.SynchronizeTransform();
      }

      this.Report(contactSolver.velocityConstraints);
    }

    private static s_impulse = new ContactImpulse();
    public Report(constraints: ContactVelocityConstraint[]): void {
      if (this.listener === null) {
        return;
      }

      for (let i: number = 0; i < this.contactCount; ++i) {
        const c: Contact = this.contacts[i];

        if (!c) { continue; }

        const vc: ContactVelocityConstraint = constraints[i];

        const impulse: ContactImpulse = Island.s_impulse;
        impulse.count = vc.pointCount;
        for (let j: number = 0; j < vc.pointCount; ++j) {
          impulse.normalImpulses[j] = vc.points[j].normalImpulse;
          impulse.tangentImpulses[j] = vc.points[j].tangentImpulse;
        }

        this.listener.PostSolve(c, impulse);
      }
    }
  }
}

