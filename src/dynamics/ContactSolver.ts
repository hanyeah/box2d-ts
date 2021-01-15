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

// Solver debugging is normally disabled because the block solver sometimes has to deal with a poorly conditioned effective mass matrix.
// #define DEBUG_SOLVER 0
namespace b2 {
  export let gBlockSolve: boolean = false;

  export class VelocityConstraintPoint {
    public readonly rA: Vec2 = new Vec2();
    public readonly rB: Vec2 = new Vec2();
    public normalImpulse: number = 0;
    public tangentImpulse: number = 0;
    public normalMass: number = 0;
    public tangentMass: number = 0;
    public velocityBias: number = 0;

    public static makeArray(length: number): VelocityConstraintPoint[] {
      return MakeArray(length, (i: number) => new VelocityConstraintPoint());
    }
  }

  export class ContactVelocityConstraint {
    public readonly points: VelocityConstraintPoint[] = VelocityConstraintPoint.makeArray(maxManifoldPoints);
    public readonly normal: Vec2 = new Vec2();
    public readonly tangent: Vec2 = new Vec2();
    public readonly normalMass: Mat22 = new Mat22();
    public readonly K: Mat22 = new Mat22();
    public indexA: number = 0;
    public indexB: number = 0;
    public invMassA: number = 0;
    public invMassB: number = 0;
    public invIA: number = 0;
    public invIB: number = 0;
    public friction: number = 0;
    public restitution: number = 0;
    public threshold: number = 0;
    public tangentSpeed: number = 0;
    public pointCount: number = 0;
    public contactIndex: number = 0;

    public static makeArray(length: number): ContactVelocityConstraint[] {
      return MakeArray(length, (i: number) => new ContactVelocityConstraint());
    }
  }

  export class ContactPositionConstraint {
    public readonly localPoints: Vec2[] = Vec2.MakeArray(maxManifoldPoints);
    public readonly localNormal: Vec2 = new Vec2();
    public readonly localPoint: Vec2 = new Vec2();
    public indexA: number = 0;
    public indexB: number = 0;
    public invMassA: number = 0;
    public invMassB: number = 0;
    public readonly localCenterA: Vec2 = new Vec2();
    public readonly localCenterB: Vec2 = new Vec2();
    public invIA: number = 0;
    public invIB: number = 0;
    public type: ManifoldType = ManifoldType.Unknown;
    public radiusA: number = 0;
    public radiusB: number = 0;
    public pointCount: number = 0;

    public static makeArray(length: number): ContactPositionConstraint[] {
      return MakeArray(length, (i: number) => new ContactPositionConstraint());
    }
  }

  export class ContactSolverDef {
    public readonly step: TimeStep = new TimeStep();
    public contacts!: Contact[];
    public count: number = 0;
    public positions!: Position[];
    public velocities!: Velocity[];
  }

  export class PositionSolverManifold {
    public readonly normal: Vec2 = new Vec2();
    public readonly point: Vec2 = new Vec2();
    public separation: number = 0;

    private static Initialize_s_pointA = new Vec2();
    private static Initialize_s_pointB = new Vec2();
    private static Initialize_s_planePoint = new Vec2();
    private static Initialize_s_clipPoint = new Vec2();
    public initialize(pc: ContactPositionConstraint, xfA: Transform, xfB: Transform, index: number): void {
      const pointA: Vec2 = PositionSolverManifold.Initialize_s_pointA;
      const pointB: Vec2 = PositionSolverManifold.Initialize_s_pointB;
      const planePoint: Vec2 = PositionSolverManifold.Initialize_s_planePoint;
      const clipPoint: Vec2 = PositionSolverManifold.Initialize_s_clipPoint;

      // DEBUG: Assert(pc.pointCount > 0);

      switch (pc.type) {
        case ManifoldType.Circles: {
          // Vec2 pointA = Mul(xfA, pc->localPoint);
          Transform.mulXV(xfA, pc.localPoint, pointA);
          // Vec2 pointB = Mul(xfB, pc->localPoints[0]);
          Transform.mulXV(xfB, pc.localPoints[0], pointB);
          // normal = pointB - pointA;
          // normal.Normalize();
          Vec2.SubVV(pointB, pointA, this.normal).selfNormalize();
          // point = 0.5f * (pointA + pointB);
          Vec2.MidVV(pointA, pointB, this.point);
          // separation = Dot(pointB - pointA, normal) - pc->radius;
          this.separation = Vec2.DotVV(Vec2.SubVV(pointB, pointA, Vec2.s_t0), this.normal) - pc.radiusA - pc.radiusB;
          break;
        }

        case ManifoldType.FaceA: {
          // normal = Mul(xfA.q, pc->localNormal);
          Rot.mulRV(xfA.q, pc.localNormal, this.normal);
          // Vec2 planePoint = Mul(xfA, pc->localPoint);
          Transform.mulXV(xfA, pc.localPoint, planePoint);

          // Vec2 clipPoint = Mul(xfB, pc->localPoints[index]);
          Transform.mulXV(xfB, pc.localPoints[index], clipPoint);
          // separation = Dot(clipPoint - planePoint, normal) - pc->radius;
          this.separation = Vec2.DotVV(Vec2.SubVV(clipPoint, planePoint, Vec2.s_t0), this.normal) - pc.radiusA - pc.radiusB;
          // point = clipPoint;
          this.point.copy(clipPoint);
          break;
        }

        case ManifoldType.FaceB: {
          // normal = Mul(xfB.q, pc->localNormal);
          Rot.mulRV(xfB.q, pc.localNormal, this.normal);
          // Vec2 planePoint = Mul(xfB, pc->localPoint);
          Transform.mulXV(xfB, pc.localPoint, planePoint);

          // Vec2 clipPoint = Mul(xfA, pc->localPoints[index]);
          Transform.mulXV(xfA, pc.localPoints[index], clipPoint);
          // separation = Dot(clipPoint - planePoint, normal) - pc->radius;
          this.separation = Vec2.DotVV(Vec2.SubVV(clipPoint, planePoint, Vec2.s_t0), this.normal) - pc.radiusA - pc.radiusB;
          // point = clipPoint;
          this.point.copy(clipPoint);

          // Ensure normal points from A to B
          // normal = -normal;
          this.normal.selfNeg();
          break;
        }
      }
    }
  }

  export class ContactSolver {
    public readonly step: TimeStep = new TimeStep();
    public positions!: Position[];
    public velocities!: Velocity[];
    public readonly positionConstraints: ContactPositionConstraint[] = ContactPositionConstraint.makeArray(1024); // TODO: Settings
    public readonly velocityConstraints: ContactVelocityConstraint[] = ContactVelocityConstraint.makeArray(1024); // TODO: Settings
    public contacts!: Contact[];
    public count: number = 0;

    public initialize(def: ContactSolverDef): ContactSolver {
      this.step.copy(def.step);
      this.count = def.count;
      // TODO:
      if (this.positionConstraints.length < this.count) {
        const new_length: number = Max(this.positionConstraints.length * 2, this.count);
        while (this.positionConstraints.length < new_length) {
          this.positionConstraints[this.positionConstraints.length] = new ContactPositionConstraint();
        }
      }
      // TODO:
      if (this.velocityConstraints.length < this.count) {
        const new_length: number = Max(this.velocityConstraints.length * 2, this.count);
        while (this.velocityConstraints.length < new_length) {
          this.velocityConstraints[this.velocityConstraints.length] = new ContactVelocityConstraint();
        }
      }
      this.positions = def.positions;
      this.velocities = def.velocities;
      this.contacts = def.contacts;

      // Initialize position independent portions of the constraints.
      for (let i: number = 0; i < this.count; ++i) {
        const contact: Contact = this.contacts[i];

        const fixtureA: Fixture = contact.fixtureA;
        const fixtureB: Fixture = contact.fixtureB;
        const shapeA: Shape = fixtureA.getShape();
        const shapeB: Shape = fixtureB.getShape();
        const radiusA: number = shapeA.radius;
        const radiusB: number = shapeB.radius;
        const bodyA: Body = fixtureA.getBody();
        const bodyB: Body = fixtureB.getBody();
        const manifold: Manifold = contact.getManifold();

        const pointCount: number = manifold.pointCount;
        // DEBUG: Assert(pointCount > 0);

        const vc: ContactVelocityConstraint = this.velocityConstraints[i];
        vc.friction = contact.friction;
        vc.restitution = contact.restitution;
        vc.threshold = contact.restitutionThreshold;
        vc.tangentSpeed = contact.tangentSpeed;
        vc.indexA = bodyA.islandIndex;
        vc.indexB = bodyB.islandIndex;
        vc.invMassA = bodyA.invMass;
        vc.invMassB = bodyB.invMass;
        vc.invIA = bodyA.invI;
        vc.invIB = bodyB.invI;
        vc.contactIndex = i;
        vc.pointCount = pointCount;
        vc.K.setZero();
        vc.normalMass.setZero();

        const pc: ContactPositionConstraint = this.positionConstraints[i];
        pc.indexA = bodyA.islandIndex;
        pc.indexB = bodyB.islandIndex;
        pc.invMassA = bodyA.invMass;
        pc.invMassB = bodyB.invMass;
        pc.localCenterA.copy(bodyA.sweep.localCenter);
        pc.localCenterB.copy(bodyB.sweep.localCenter);
        pc.invIA = bodyA.invI;
        pc.invIB = bodyB.invI;
        pc.localNormal.copy(manifold.localNormal);
        pc.localPoint.copy(manifold.localPoint);
        pc.pointCount = pointCount;
        pc.radiusA = radiusA;
        pc.radiusB = radiusB;
        pc.type = manifold.type;

        for (let j: number = 0; j < pointCount; ++j) {
          const cp: ManifoldPoint = manifold.points[j];
          const vcp: VelocityConstraintPoint = vc.points[j];

          if (this.step.warmStarting) {
            vcp.normalImpulse = this.step.dtRatio * cp.normalImpulse;
            vcp.tangentImpulse = this.step.dtRatio * cp.tangentImpulse;
          } else {
            vcp.normalImpulse = 0;
            vcp.tangentImpulse = 0;
          }

          vcp.rA.setZero();
          vcp.rB.setZero();
          vcp.normalMass = 0;
          vcp.tangentMass = 0;
          vcp.velocityBias = 0;

          pc.localPoints[j].copy(cp.localPoint);
        }
      }

      return this;
    }

    private static initializeVelocityConstraints_s_xfA = new Transform();
    private static initializeVelocityConstraints_s_xfB = new Transform();
    private static initializeVelocityConstraints_s_worldManifold = new WorldManifold();
    public initializeVelocityConstraints(): void {
      const xfA: Transform = ContactSolver.initializeVelocityConstraints_s_xfA;
      const xfB: Transform = ContactSolver.initializeVelocityConstraints_s_xfB;
      const worldManifold: WorldManifold = ContactSolver.initializeVelocityConstraints_s_worldManifold;

      const k_maxConditionNumber: number = 1000;

      for (let i: number = 0; i < this.count; ++i) {
        const vc: ContactVelocityConstraint = this.velocityConstraints[i];
        const pc: ContactPositionConstraint = this.positionConstraints[i];

        const radiusA: number = pc.radiusA;
        const radiusB: number = pc.radiusB;
        const manifold: Manifold = this.contacts[vc.contactIndex].getManifold();

        const indexA: number = vc.indexA;
        const indexB: number = vc.indexB;

        const mA: number = vc.invMassA;
        const mB: number = vc.invMassB;
        const iA: number = vc.invIA;
        const iB: number = vc.invIB;
        const localCenterA: Vec2 = pc.localCenterA;
        const localCenterB: Vec2 = pc.localCenterB;

        const cA: Vec2 = this.positions[indexA].c;
        const aA: number = this.positions[indexA].a;
        const vA: Vec2 = this.velocities[indexA].v;
        const wA: number = this.velocities[indexA].w;

        const cB: Vec2 = this.positions[indexB].c;
        const aB: number = this.positions[indexB].a;
        const vB: Vec2 = this.velocities[indexB].v;
        const wB: number = this.velocities[indexB].w;

        // DEBUG: Assert(manifold.pointCount > 0);

        xfA.q.setAngle(aA);
        xfB.q.setAngle(aB);
        Vec2.SubVV(cA, Rot.mulRV(xfA.q, localCenterA, Vec2.s_t0), xfA.p);
        Vec2.SubVV(cB, Rot.mulRV(xfB.q, localCenterB, Vec2.s_t0), xfB.p);

        worldManifold.initialize(manifold, xfA, radiusA, xfB, radiusB);

        vc.normal.copy(worldManifold.normal);
        Vec2.CrossVOne(vc.normal, vc.tangent); // compute from normal

        const pointCount: number = vc.pointCount;
        for (let j: number = 0; j < pointCount; ++j) {
          const vcp: VelocityConstraintPoint = vc.points[j];

          // vcp->rA = worldManifold.points[j] - cA;
          Vec2.SubVV(worldManifold.points[j], cA, vcp.rA);
          // vcp->rB = worldManifold.points[j] - cB;
          Vec2.SubVV(worldManifold.points[j], cB, vcp.rB);

          const rnA: number = Vec2.CrossVV(vcp.rA, vc.normal);
          const rnB: number = Vec2.CrossVV(vcp.rB, vc.normal);

          const kNormal: number = mA + mB + iA * rnA * rnA + iB * rnB * rnB;

          vcp.normalMass = kNormal > 0 ? 1 / kNormal : 0;

          // Vec2 tangent = Cross(vc->normal, 1.0f);
          const tangent: Vec2 = vc.tangent; // precomputed from normal

          const rtA: number = Vec2.CrossVV(vcp.rA, tangent);
          const rtB: number = Vec2.CrossVV(vcp.rB, tangent);

          const kTangent: number = mA + mB + iA * rtA * rtA + iB * rtB * rtB;

          vcp.tangentMass = kTangent > 0 ? 1 / kTangent : 0;

          // Setup a velocity bias for restitution.
          vcp.velocityBias = 0;
          // float32 vRel = Dot(vc->normal, vB + Cross(wB, vcp->rB) - vA - Cross(wA, vcp->rA));
          const vRel: number = Vec2.DotVV(
            vc.normal,
            Vec2.SubVV(
              Vec2.AddVCrossSV(vB, wB, vcp.rB, Vec2.s_t0),
              Vec2.AddVCrossSV(vA, wA, vcp.rA, Vec2.s_t1),
              Vec2.s_t0));
          if (vRel < -vc.threshold) {
            vcp.velocityBias += (-vc.restitution * vRel);
          }
        }

        // If we have two points, then prepare the block solver.
        if (vc.pointCount === 2 && gBlockSolve) {
          const vcp1: VelocityConstraintPoint = vc.points[0];
          const vcp2: VelocityConstraintPoint = vc.points[1];

          const rn1A: number = Vec2.CrossVV(vcp1.rA, vc.normal);
          const rn1B: number = Vec2.CrossVV(vcp1.rB, vc.normal);
          const rn2A: number = Vec2.CrossVV(vcp2.rA, vc.normal);
          const rn2B: number = Vec2.CrossVV(vcp2.rB, vc.normal);

          const k11: number = mA + mB + iA * rn1A * rn1A + iB * rn1B * rn1B;
          const k22: number = mA + mB + iA * rn2A * rn2A + iB * rn2B * rn2B;
          const k12: number = mA + mB + iA * rn1A * rn2A + iB * rn1B * rn2B;

          // Ensure a reasonable condition number.
          // float32 k_maxConditionNumber = 1000.0f;
          if (k11 * k11 < k_maxConditionNumber * (k11 * k22 - k12 * k12)) {
            // K is safe to invert.
            vc.K.ex.set(k11, k12);
            vc.K.ey.set(k12, k22);
            vc.K.getInverse(vc.normalMass);
          } else {
            // The constraints are redundant, just use one.
            // TODO_ERIN use deepest?
            vc.pointCount = 1;
          }
        }
      }
    }

    private static warmStart_s_P = new Vec2();
    public warmStart(): void {
      const P: Vec2 = ContactSolver.warmStart_s_P;

      // Warm start.
      for (let i: number = 0; i < this.count; ++i) {
        const vc: ContactVelocityConstraint = this.velocityConstraints[i];

        const indexA: number = vc.indexA;
        const indexB: number = vc.indexB;
        const mA: number = vc.invMassA;
        const iA: number = vc.invIA;
        const mB: number = vc.invMassB;
        const iB: number = vc.invIB;
        const pointCount: number = vc.pointCount;

        const vA: Vec2 = this.velocities[indexA].v;
        let wA: number = this.velocities[indexA].w;
        const vB: Vec2 = this.velocities[indexB].v;
        let wB: number = this.velocities[indexB].w;

        const normal: Vec2 = vc.normal;
        // Vec2 tangent = Cross(normal, 1.0f);
        const tangent: Vec2 = vc.tangent; // precomputed from normal

        for (let j: number = 0; j < pointCount; ++j) {
          const vcp: VelocityConstraintPoint = vc.points[j];
          // Vec2 P = vcp->normalImpulse * normal + vcp->tangentImpulse * tangent;
          Vec2.AddVV(
            Vec2.MulSV(vcp.normalImpulse, normal, Vec2.s_t0),
            Vec2.MulSV(vcp.tangentImpulse, tangent, Vec2.s_t1),
            P);
          // wA -= iA * Cross(vcp->rA, P);
          wA -= iA * Vec2.CrossVV(vcp.rA, P);
          // vA -= mA * P;
          vA.selfMulSub(mA, P);
          // wB += iB * Cross(vcp->rB, P);
          wB += iB * Vec2.CrossVV(vcp.rB, P);
          // vB += mB * P;
          vB.selfMulAdd(mB, P);
        }

        // this.velocities[indexA].v = vA;
        this.velocities[indexA].w = wA;
        // this.velocities[indexB].v = vB;
        this.velocities[indexB].w = wB;
      }
    }

    private static solveVelocityConstraints_s_dv = new Vec2();
    private static solveVelocityConstraints_s_dv1 = new Vec2();
    private static solveVelocityConstraints_s_dv2 = new Vec2();
    private static solveVelocityConstraints_s_P = new Vec2();
    private static solveVelocityConstraints_s_a = new Vec2();
    private static solveVelocityConstraints_s_b = new Vec2();
    private static solveVelocityConstraints_s_x = new Vec2();
    private static solveVelocityConstraints_s_d = new Vec2();
    private static solveVelocityConstraints_s_P1 = new Vec2();
    private static solveVelocityConstraints_s_P2 = new Vec2();
    private static solveVelocityConstraints_s_P1P2 = new Vec2();
    public solveVelocityConstraints(): void {
      const dv: Vec2 = ContactSolver.solveVelocityConstraints_s_dv;
      const dv1: Vec2 = ContactSolver.solveVelocityConstraints_s_dv1;
      const dv2: Vec2 = ContactSolver.solveVelocityConstraints_s_dv2;
      const P: Vec2 = ContactSolver.solveVelocityConstraints_s_P;
      const a: Vec2 = ContactSolver.solveVelocityConstraints_s_a;
      const b: Vec2 = ContactSolver.solveVelocityConstraints_s_b;
      const x: Vec2 = ContactSolver.solveVelocityConstraints_s_x;
      const d: Vec2 = ContactSolver.solveVelocityConstraints_s_d;
      const P1: Vec2 = ContactSolver.solveVelocityConstraints_s_P1;
      const P2: Vec2 = ContactSolver.solveVelocityConstraints_s_P2;
      const P1P2: Vec2 = ContactSolver.solveVelocityConstraints_s_P1P2;

      for (let i: number = 0; i < this.count; ++i) {
        const vc: ContactVelocityConstraint = this.velocityConstraints[i];

        const indexA: number = vc.indexA;
        const indexB: number = vc.indexB;
        const mA: number = vc.invMassA;
        const iA: number = vc.invIA;
        const mB: number = vc.invMassB;
        const iB: number = vc.invIB;
        const pointCount: number = vc.pointCount;

        const vA: Vec2 = this.velocities[indexA].v;
        let wA: number = this.velocities[indexA].w;
        const vB: Vec2 = this.velocities[indexB].v;
        let wB: number = this.velocities[indexB].w;

        // Vec2 normal = vc->normal;
        const normal: Vec2 = vc.normal;
        // Vec2 tangent = Cross(normal, 1.0f);
        const tangent: Vec2 = vc.tangent; // precomputed from normal
        const friction: number = vc.friction;

        // DEBUG: Assert(pointCount === 1 || pointCount === 2);

        // Solve tangent constraints first because non-penetration is more important
        // than friction.
        for (let j: number = 0; j < pointCount; ++j) {
          const vcp: VelocityConstraintPoint = vc.points[j];

          // Relative velocity at contact
          // Vec2 dv = vB + Cross(wB, vcp->rB) - vA - Cross(wA, vcp->rA);
          Vec2.SubVV(
            Vec2.AddVCrossSV(vB, wB, vcp.rB, Vec2.s_t0),
            Vec2.AddVCrossSV(vA, wA, vcp.rA, Vec2.s_t1),
            dv);

          // Compute tangent force
          // float32 vt = Dot(dv, tangent) - vc->tangentSpeed;
          const vt: number = Vec2.DotVV(dv, tangent) - vc.tangentSpeed;
          let lambda: number = vcp.tangentMass * (-vt);

          // Clamp the accumulated force
          const maxFriction: number = friction * vcp.normalImpulse;
          const newImpulse: number = clamp(vcp.tangentImpulse + lambda, (-maxFriction), maxFriction);
          lambda = newImpulse - vcp.tangentImpulse;
          vcp.tangentImpulse = newImpulse;

          // Apply contact impulse
          // Vec2 P = lambda * tangent;
          Vec2.MulSV(lambda, tangent, P);

          // vA -= mA * P;
          vA.selfMulSub(mA, P);
          // wA -= iA * Cross(vcp->rA, P);
          wA -= iA * Vec2.CrossVV(vcp.rA, P);

          // vB += mB * P;
          vB.selfMulAdd(mB, P);
          // wB += iB * Cross(vcp->rB, P);
          wB += iB * Vec2.CrossVV(vcp.rB, P);
        }

        // Solve normal constraints
        if (vc.pointCount === 1 || gBlockSolve === false) {
          for (let j = 0; j < pointCount; ++j) {
            const vcp: VelocityConstraintPoint = vc.points[j];

            // Relative velocity at contact
            // Vec2 dv = vB + Cross(wB, vcp->rB) - vA - Cross(wA, vcp->rA);
            Vec2.SubVV(
              Vec2.AddVCrossSV(vB, wB, vcp.rB, Vec2.s_t0),
              Vec2.AddVCrossSV(vA, wA, vcp.rA, Vec2.s_t1),
              dv);

            // Compute normal impulse
            // float32 vn = Dot(dv, normal);
            const vn: number = Vec2.DotVV(dv, normal);
            let lambda: number = (-vcp.normalMass * (vn - vcp.velocityBias));

            // Clamp the accumulated impulse
            // float32 newImpulse = Max(vcp->normalImpulse + lambda, 0.0f);
            const newImpulse: number = Max(vcp.normalImpulse + lambda, 0);
            lambda = newImpulse - vcp.normalImpulse;
            vcp.normalImpulse = newImpulse;

            // Apply contact impulse
            // Vec2 P = lambda * normal;
            Vec2.MulSV(lambda, normal, P);
            // vA -= mA * P;
            vA.selfMulSub(mA, P);
            // wA -= iA * Cross(vcp->rA, P);
            wA -= iA * Vec2.CrossVV(vcp.rA, P);

            // vB += mB * P;
            vB.selfMulAdd(mB, P);
            // wB += iB * Cross(vcp->rB, P);
            wB += iB * Vec2.CrossVV(vcp.rB, P);
          }
        } else {
          // Block solver developed in collaboration with Dirk Gregorius (back in 01/07 on Box2D_Lite).
          // Build the mini LCP for this contact patch
          //
          // vn = A * x + b, vn >= 0, x >= 0 and vn_i * x_i = 0 with i = 1..2
          //
          // A = J * W * JT and J = ( -n, -r1 x n, n, r2 x n )
          // b = vn0 - velocityBias
          //
          // The system is solved using the "Total enumeration method" (s. Murty). The complementary constraint vn_i * x_i
          // implies that we must have in any solution either vn_i = 0 or x_i = 0. So for the 2D contact problem the cases
          // vn1 = 0 and vn2 = 0, x1 = 0 and x2 = 0, x1 = 0 and vn2 = 0, x2 = 0 and vn1 = 0 need to be tested. The first valid
          // solution that satisfies the problem is chosen.
          //
          // In order to account of the accumulated impulse 'a' (because of the iterative nature of the solver which only requires
          // that the accumulated impulse is clamped and not the incremental impulse) we change the impulse variable (x_i).
          //
          // Substitute:
          //
          // x = a + d
          //
          // a := old total impulse
          // x := new total impulse
          // d := incremental impulse
          //
          // For the current iteration we extend the formula for the incremental impulse
          // to compute the new total impulse:
          //
          // vn = A * d + b
          //    = A * (x - a) + b
          //    = A * x + b - A * a
          //    = A * x + b'
          // b' = b - A * a;

          const cp1: VelocityConstraintPoint = vc.points[0];
          const cp2: VelocityConstraintPoint = vc.points[1];

          // Vec2 a(cp1->normalImpulse, cp2->normalImpulse);
          a.set(cp1.normalImpulse, cp2.normalImpulse);
          // DEBUG: Assert(a.x >= 0 && a.y >= 0);

          // Relative velocity at contact
          // Vec2 dv1 = vB + Cross(wB, cp1->rB) - vA - Cross(wA, cp1->rA);
          Vec2.SubVV(
            Vec2.AddVCrossSV(vB, wB, cp1.rB, Vec2.s_t0),
            Vec2.AddVCrossSV(vA, wA, cp1.rA, Vec2.s_t1),
            dv1);
          // Vec2 dv2 = vB + Cross(wB, cp2->rB) - vA - Cross(wA, cp2->rA);
          Vec2.SubVV(
            Vec2.AddVCrossSV(vB, wB, cp2.rB, Vec2.s_t0),
            Vec2.AddVCrossSV(vA, wA, cp2.rA, Vec2.s_t1),
            dv2);

          // Compute normal velocity
          // float32 vn1 = Dot(dv1, normal);
          let vn1: number = Vec2.DotVV(dv1, normal);
          // float32 vn2 = Dot(dv2, normal);
          let vn2: number = Vec2.DotVV(dv2, normal);

          // Vec2 b;
          b.x = vn1 - cp1.velocityBias;
          b.y = vn2 - cp2.velocityBias;

          // Compute b'
          // b -= Mul(vc->K, a);
          b.selfSub(Mat22.mulMV(vc.K, a, Vec2.s_t0));

          /*
          #if DEBUG_SOLVER === 1
          const k_errorTol: number = 0.001;
          #endif
          */

          for (; ; ) {
            //
            // Case 1: vn = 0
            //
            // 0 = A * x + b'
            //
            // Solve for x:
            //
            // x = - inv(A) * b'
            //
            // Vec2 x = - Mul(vc->normalMass, b);
            Mat22.mulMV(vc.normalMass, b, x).selfNeg();

            if (x.x >= 0 && x.y >= 0) {
              // Get the incremental impulse
              // Vec2 d = x - a;
              Vec2.SubVV(x, a, d);

              // Apply incremental impulse
              // Vec2 P1 = d.x * normal;
              Vec2.MulSV(d.x, normal, P1);
              // Vec2 P2 = d.y * normal;
              Vec2.MulSV(d.y, normal, P2);
              Vec2.AddVV(P1, P2, P1P2);
              // vA -= mA * (P1 + P2);
              vA.selfMulSub(mA, P1P2);
              // wA -= iA * (Cross(cp1->rA, P1) + Cross(cp2->rA, P2));
              wA -= iA * (Vec2.CrossVV(cp1.rA, P1) + Vec2.CrossVV(cp2.rA, P2));

              // vB += mB * (P1 + P2);
              vB.selfMulAdd(mB, P1P2);
              // wB += iB * (Cross(cp1->rB, P1) + Cross(cp2->rB, P2));
              wB += iB * (Vec2.CrossVV(cp1.rB, P1) + Vec2.CrossVV(cp2.rB, P2));

              // Accumulate
              cp1.normalImpulse = x.x;
              cp2.normalImpulse = x.y;

              /*
              #if DEBUG_SOLVER === 1
              // Postconditions
              dv1 = vB + Cross(wB, cp1->rB) - vA - Cross(wA, cp1->rA);
              dv2 = vB + Cross(wB, cp2->rB) - vA - Cross(wA, cp2->rA);

              // Compute normal velocity
              vn1 = Dot(dv1, normal);
              vn2 = Dot(dv2, normal);

              Assert(Abs(vn1 - cp1->velocityBias) < k_errorTol);
              Assert(Abs(vn2 - cp2->velocityBias) < k_errorTol);
              #endif
              */
              break;
            }

            //
            // Case 2: vn1 = 0 and x2 = 0
            //
            //   0 = a11 * x1 + a12 * 0 + b1'
            // vn2 = a21 * x1 + a22 * 0 + '
            //
            x.x = (-cp1.normalMass * b.x);
            x.y = 0;
            vn1 = 0;
            vn2 = vc.K.ex.y * x.x + b.y;

            if (x.x >= 0 && vn2 >= 0) {
              // Get the incremental impulse
              // Vec2 d = x - a;
              Vec2.SubVV(x, a, d);

              // Apply incremental impulse
              // Vec2 P1 = d.x * normal;
              Vec2.MulSV(d.x, normal, P1);
              // Vec2 P2 = d.y * normal;
              Vec2.MulSV(d.y, normal, P2);
              Vec2.AddVV(P1, P2, P1P2);
              // vA -= mA * (P1 + P2);
              vA.selfMulSub(mA, P1P2);
              // wA -= iA * (Cross(cp1->rA, P1) + Cross(cp2->rA, P2));
              wA -= iA * (Vec2.CrossVV(cp1.rA, P1) + Vec2.CrossVV(cp2.rA, P2));

              // vB += mB * (P1 + P2);
              vB.selfMulAdd(mB, P1P2);
              // wB += iB * (Cross(cp1->rB, P1) + Cross(cp2->rB, P2));
              wB += iB * (Vec2.CrossVV(cp1.rB, P1) + Vec2.CrossVV(cp2.rB, P2));

              // Accumulate
              cp1.normalImpulse = x.x;
              cp2.normalImpulse = x.y;

              /*
              #if DEBUG_SOLVER === 1
              // Postconditions
              dv1 = vB + Cross(wB, cp1->rB) - vA - Cross(wA, cp1->rA);

              // Compute normal velocity
              vn1 = Dot(dv1, normal);

              Assert(Abs(vn1 - cp1->velocityBias) < k_errorTol);
              #endif
              */
              break;
            }

            //
            // Case 3: vn2 = 0 and x1 = 0
            //
            // vn1 = a11 * 0 + a12 * x2 + b1'
            //   0 = a21 * 0 + a22 * x2 + '
            //
            x.x = 0;
            x.y = (-cp2.normalMass * b.y);
            vn1 = vc.K.ey.x * x.y + b.x;
            vn2 = 0;

            if (x.y >= 0 && vn1 >= 0) {
              // Resubstitute for the incremental impulse
              // Vec2 d = x - a;
              Vec2.SubVV(x, a, d);

              // Apply incremental impulse
              // Vec2 P1 = d.x * normal;
              Vec2.MulSV(d.x, normal, P1);
              // Vec2 P2 = d.y * normal;
              Vec2.MulSV(d.y, normal, P2);
              Vec2.AddVV(P1, P2, P1P2);
              // vA -= mA * (P1 + P2);
              vA.selfMulSub(mA, P1P2);
              // wA -= iA * (Cross(cp1->rA, P1) + Cross(cp2->rA, P2));
              wA -= iA * (Vec2.CrossVV(cp1.rA, P1) + Vec2.CrossVV(cp2.rA, P2));

              // vB += mB * (P1 + P2);
              vB.selfMulAdd(mB, P1P2);
              // wB += iB * (Cross(cp1->rB, P1) + Cross(cp2->rB, P2));
              wB += iB * (Vec2.CrossVV(cp1.rB, P1) + Vec2.CrossVV(cp2.rB, P2));

              // Accumulate
              cp1.normalImpulse = x.x;
              cp2.normalImpulse = x.y;

              /*
              #if DEBUG_SOLVER === 1
              // Postconditions
              dv2 = vB + Cross(wB, cp2->rB) - vA - Cross(wA, cp2->rA);

              // Compute normal velocity
              vn2 = Dot(dv2, normal);

              Assert(Abs(vn2 - cp2->velocityBias) < k_errorTol);
              #endif
              */
              break;
            }

            //
            // Case 4: x1 = 0 and x2 = 0
            //
            // vn1 = b1
            // vn2 = ;
            x.x = 0;
            x.y = 0;
            vn1 = b.x;
            vn2 = b.y;

            if (vn1 >= 0 && vn2 >= 0) {
              // Resubstitute for the incremental impulse
              // Vec2 d = x - a;
              Vec2.SubVV(x, a, d);

              // Apply incremental impulse
              // Vec2 P1 = d.x * normal;
              Vec2.MulSV(d.x, normal, P1);
              // Vec2 P2 = d.y * normal;
              Vec2.MulSV(d.y, normal, P2);
              Vec2.AddVV(P1, P2, P1P2);
              // vA -= mA * (P1 + P2);
              vA.selfMulSub(mA, P1P2);
              // wA -= iA * (Cross(cp1->rA, P1) + Cross(cp2->rA, P2));
              wA -= iA * (Vec2.CrossVV(cp1.rA, P1) + Vec2.CrossVV(cp2.rA, P2));

              // vB += mB * (P1 + P2);
              vB.selfMulAdd(mB, P1P2);
              // wB += iB * (Cross(cp1->rB, P1) + Cross(cp2->rB, P2));
              wB += iB * (Vec2.CrossVV(cp1.rB, P1) + Vec2.CrossVV(cp2.rB, P2));

              // Accumulate
              cp1.normalImpulse = x.x;
              cp2.normalImpulse = x.y;

              break;
            }

            // No solution, give up. This is hit sometimes, but it doesn't seem to matter.
            break;
          }
        }

        // this.velocities[indexA].v = vA;
        this.velocities[indexA].w = wA;
        // this.velocities[indexB].v = vB;
        this.velocities[indexB].w = wB;
      }
    }

    public storeImpulses(): void {
      for (let i: number = 0; i < this.count; ++i) {
        const vc: ContactVelocityConstraint = this.velocityConstraints[i];
        const manifold: Manifold = this.contacts[vc.contactIndex].getManifold();

        for (let j: number = 0; j < vc.pointCount; ++j) {
          manifold.points[j].normalImpulse = vc.points[j].normalImpulse;
          manifold.points[j].tangentImpulse = vc.points[j].tangentImpulse;
        }
      }
    }

    private static solvePositionConstraints_s_xfA = new Transform();
    private static solvePositionConstraints_s_xfB = new Transform();
    private static solvePositionConstraints_s_psm = new PositionSolverManifold();
    private static solvePositionConstraints_s_rA = new Vec2();
    private static solvePositionConstraints_s_rB = new Vec2();
    private static solvePositionConstraints_s_P = new Vec2();
    public solvePositionConstraints(): boolean {
      const xfA: Transform = ContactSolver.solvePositionConstraints_s_xfA;
      const xfB: Transform = ContactSolver.solvePositionConstraints_s_xfB;
      const psm: PositionSolverManifold = ContactSolver.solvePositionConstraints_s_psm;
      const rA: Vec2 = ContactSolver.solvePositionConstraints_s_rA;
      const rB: Vec2 = ContactSolver.solvePositionConstraints_s_rB;
      const P: Vec2 = ContactSolver.solvePositionConstraints_s_P;

      let minSeparation: number = 0;

      for (let i: number = 0; i < this.count; ++i) {
        const pc: ContactPositionConstraint = this.positionConstraints[i];

        const indexA: number = pc.indexA;
        const indexB: number = pc.indexB;
        const localCenterA: Vec2 = pc.localCenterA;
        const mA: number = pc.invMassA;
        const iA: number = pc.invIA;
        const localCenterB: Vec2 = pc.localCenterB;
        const mB: number = pc.invMassB;
        const iB: number = pc.invIB;
        const pointCount: number = pc.pointCount;

        const cA: Vec2 = this.positions[indexA].c;
        let aA: number = this.positions[indexA].a;

        const cB: Vec2 = this.positions[indexB].c;
        let aB: number = this.positions[indexB].a;

        // Solve normal constraints
        for (let j: number = 0; j < pointCount; ++j) {
          xfA.q.setAngle(aA);
          xfB.q.setAngle(aB);
          Vec2.SubVV(cA, Rot.mulRV(xfA.q, localCenterA, Vec2.s_t0), xfA.p);
          Vec2.SubVV(cB, Rot.mulRV(xfB.q, localCenterB, Vec2.s_t0), xfB.p);

          psm.initialize(pc, xfA, xfB, j);
          const normal: Vec2 = psm.normal;

          const point: Vec2 = psm.point;
          const separation: number = psm.separation;

          // Vec2 rA = point - cA;
          Vec2.SubVV(point, cA, rA);
          // Vec2 rB = point - cB;
          Vec2.SubVV(point, cB, rB);

          // Track max constraint error.
          minSeparation = Min(minSeparation, separation);

          // Prevent large corrections and allow slop.
          const C: number = clamp(baumgarte * (separation + linearSlop), (-maxLinearCorrection), 0);

          // Compute the effective mass.
          // float32 rnA = Cross(rA, normal);
          const rnA: number = Vec2.CrossVV(rA, normal);
          // float32 rnB = Cross(rB, normal);
          const rnB: number = Vec2.CrossVV(rB, normal);
          // float32 K = mA + mB + iA * rnA * rnA + iB * rnB * rnB;
          const K: number = mA + mB + iA * rnA * rnA + iB * rnB * rnB;

          // Compute normal impulse
          const impulse: number = K > 0 ? - C / K : 0;

          // Vec2 P = impulse * normal;
          Vec2.MulSV(impulse, normal, P);

          // cA -= mA * P;
          cA.selfMulSub(mA, P);
          // aA -= iA * Cross(rA, P);
          aA -= iA * Vec2.CrossVV(rA, P);

          // cB += mB * P;
          cB.selfMulAdd(mB, P);
          // aB += iB * Cross(rB, P);
          aB += iB * Vec2.CrossVV(rB, P);
        }

        // this.positions[indexA].c = cA;
        this.positions[indexA].a = aA;

        // this.positions[indexB].c = cB;
        this.positions[indexB].a = aB;
      }

      // We can't expect minSpeparation >= -linearSlop because we don't
      // push the separation above -linearSlop.
      return minSeparation > (-3 * linearSlop);
    }

    private static solveTOIPositionConstraints_s_xfA = new Transform();
    private static solveTOIPositionConstraints_s_xfB = new Transform();
    private static solveTOIPositionConstraints_s_psm = new PositionSolverManifold();
    private static solveTOIPositionConstraints_s_rA = new Vec2();
    private static solveTOIPositionConstraints_s_rB = new Vec2();
    private static solveTOIPositionConstraints_s_P = new Vec2();
    public solveTOIPositionConstraints(toiIndexA: number, toiIndexB: number): boolean {
      const xfA: Transform = ContactSolver.solveTOIPositionConstraints_s_xfA;
      const xfB: Transform = ContactSolver.solveTOIPositionConstraints_s_xfB;
      const psm: PositionSolverManifold = ContactSolver.solveTOIPositionConstraints_s_psm;
      const rA: Vec2 = ContactSolver.solveTOIPositionConstraints_s_rA;
      const rB: Vec2 = ContactSolver.solveTOIPositionConstraints_s_rB;
      const P: Vec2 = ContactSolver.solveTOIPositionConstraints_s_P;

      let minSeparation: number = 0;

      for (let i: number = 0; i < this.count; ++i) {
        const pc: ContactPositionConstraint = this.positionConstraints[i];

        const indexA: number = pc.indexA;
        const indexB: number = pc.indexB;
        const localCenterA: Vec2 = pc.localCenterA;
        const localCenterB: Vec2 = pc.localCenterB;
        const pointCount: number = pc.pointCount;

        let mA: number = 0;
        let iA: number = 0;
        if (indexA === toiIndexA || indexA === toiIndexB) {
          mA = pc.invMassA;
          iA = pc.invIA;
        }

        let mB: number = 0;
        let iB: number = 0;
        if (indexB === toiIndexA || indexB === toiIndexB) {
          mB = pc.invMassB;
          iB = pc.invIB;
        }

        const cA: Vec2 = this.positions[indexA].c;
        let aA: number = this.positions[indexA].a;

        const cB: Vec2 = this.positions[indexB].c;
        let aB: number = this.positions[indexB].a;

        // Solve normal constraints
        for (let j: number = 0; j < pointCount; ++j) {
          xfA.q.setAngle(aA);
          xfB.q.setAngle(aB);
          Vec2.SubVV(cA, Rot.mulRV(xfA.q, localCenterA, Vec2.s_t0), xfA.p);
          Vec2.SubVV(cB, Rot.mulRV(xfB.q, localCenterB, Vec2.s_t0), xfB.p);

          psm.initialize(pc, xfA, xfB, j);
          const normal: Vec2 = psm.normal;

          const point: Vec2 = psm.point;
          const separation: number = psm.separation;

          // Vec2 rA = point - cA;
          Vec2.SubVV(point, cA, rA);
          // Vec2 rB = point - cB;
          Vec2.SubVV(point, cB, rB);

          // Track max constraint error.
          minSeparation = Min(minSeparation, separation);

          // Prevent large corrections and allow slop.
          const C: number = clamp(toiBaumgarte * (separation + linearSlop), (-maxLinearCorrection), 0);

          // Compute the effective mass.
          // float32 rnA = Cross(rA, normal);
          const rnA: number = Vec2.CrossVV(rA, normal);
          // float32 rnB = Cross(rB, normal);
          const rnB: number = Vec2.CrossVV(rB, normal);
          // float32 K = mA + mB + iA * rnA * rnA + iB * rnB * rnB;
          const K: number = mA + mB + iA * rnA * rnA + iB * rnB * rnB;

          // Compute normal impulse
          const impulse: number = K > 0 ? - C / K : 0;

          // Vec2 P = impulse * normal;
          Vec2.MulSV(impulse, normal, P);

          // cA -= mA * P;
          cA.selfMulSub(mA, P);
          // aA -= iA * Cross(rA, P);
          aA -= iA * Vec2.CrossVV(rA, P);

          // cB += mB * P;
          cB.selfMulAdd(mB, P);
          // aB += iB * Cross(rB, P);
          aB += iB * Vec2.CrossVV(rB, P);
        }

        // this.positions[indexA].c = cA;
        this.positions[indexA].a = aA;

        // this.positions[indexB].c = cB;
        this.positions[indexB].a = aB;
      }

      // We can't expect minSpeparation >= -linearSlop because we don't
      // push the separation above -linearSlop.
      return minSeparation >= -1.5 * linearSlop;
    }
  }

}
