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
  export interface IPrismaticJointDef extends IJointDef {
    localAnchorA?: XY;

    localAnchorB?: XY;

    localAxisA?: XY;

    referenceAngle?: number;

    enableLimit?: boolean;

    lowerTranslation?: number;

    upperTranslation?: number;

    enableMotor?: boolean;

    maxMotorForce?: number;

    motorSpeed?: number;
  }

/// Prismatic joint definition. This requires defining a line of
/// motion using an axis and an anchor point. The definition uses local
/// anchor points and a local axis so that the initial configuration
/// can violate the constraint slightly. The joint translation is zero
/// when the local anchor points coincide in world space. Using local
/// anchors and a local axis helps when saving and loading a game.
  export class PrismaticJointDef extends JointDef implements IPrismaticJointDef {
    public readonly localAnchorA: Vec2 = new Vec2();

    public readonly localAnchorB: Vec2 = new Vec2();

    public readonly localAxisA: Vec2 = new Vec2(1, 0);

    public referenceAngle: number = 0;

    public enableLimit = false;

    public lowerTranslation: number = 0;

    public upperTranslation: number = 0;

    public enableMotor = false;

    public maxMotorForce: number = 0;

    public motorSpeed: number = 0;

    constructor() {
      super(JointType.PrismaticJoint);
    }

    public initialize(bA: Body, bB: Body, anchor1: XY, anchor2: XY, axis: XY): void {
      this.bodyA = bA;
      this.bodyB = bB;
      this.localAnchorA.copy(anchor1);
      this.localAnchorB.copy(anchor2);
      this.localAxisA.copy(axis);
      this.referenceAngle = this.bodyB.getAngle() - this.bodyA.getAngle();
    }
  }

// Linear constraint (point-to-line)
// d = p2 - p1 = x2 + r2 - x1 - r1
// C = dot(perp, d)
// Cdot = dot(d, cross(w1, perp)) + dot(perp, v2 + cross(w2, r2) - v1 - cross(w1, r1))
//      = -dot(perp, v1) - dot(cross(d + r1, perp), w1) + dot(perp, v2) + dot(cross(r2, perp), v2)
// J = [-perp, -cross(d + r1, perp), perp, cross(r2,perp)]
//
// Angular constraint
// C = a2 - a1 + a_initial
// Cdot = w2 - w1
// J = [0 0 -1 0 0 1]
//
// K = J * invM * JT
//
// J = [-a -s1 a s2]
//     [0  -1  0  1]
// a = perp
// s1 = cross(d + r1, a) = cross(p2 - x1, a)
// s2 = cross(r2, a) = cross(p2 - x2, a)

// Motor/Limit linear constraint
// C = dot(ax1, d)
// Cdot = -dot(ax1, v1) - dot(cross(d + r1, ax1), w1) + dot(ax1, v2) + dot(cross(r2, ax1), v2)
// J = [-ax1 -cross(d+r1,ax1) ax1 cross(r2,ax1)]

// Predictive limit is applied even when the limit is not active.
// Prevents a constraint speed that can lead to a constraint error in one time step.
// Want C2 = C1 + h * Cdot >= 0
// Or:
// Cdot + C1/h >= 0
// I do not apply a negative constraint error because that is handled in position correction.
// So:
// Cdot + max(C1, 0)/h >= 0

// Block Solver
// We develop a block solver that includes the angular and linear constraints. This makes the limit stiffer.
//
// The Jacobian has 2 rows:
// J = [-uT -s1 uT s2] // linear
//     [0   -1   0  1] // angular
//
// u = perp
// s1 = cross(d + r1, u), s2 = cross(r2, u)
// a1 = cross(d + r1, v), a2 = cross(r2, v)

  export class PrismaticJoint extends Joint {
    public readonly localAnchorA: Vec2 = new Vec2();
    public readonly localAnchorB: Vec2 = new Vec2();
    public readonly localXAxisA: Vec2 = new Vec2();
    public readonly localYAxisA: Vec2 = new Vec2();
    public referenceAngle: number = 0;
    public readonly impulse: Vec2 = new Vec2(0, 0);
    public motorImpulse: number = 0;
    public lowerImpulse: number = 0;
    public upperImpulse: number = 0;
    public lowerTranslation: number = 0;
    public upperTranslation: number = 0;
    public maxMotorForce: number = 0;
    public motorSpeed: number = 0;
    public enableLimit: boolean = false;
    public enableMotor: boolean = false;

    // Solver temp
    public indexA: number = 0;
    public indexB: number = 0;
    public readonly localCenterA: Vec2 = new Vec2();
    public readonly localCenterB: Vec2 = new Vec2();
    public invMassA: number = 0;
    public invMassB: number = 0;
    public invIA: number = 0;
    public invIB: number = 0;
    public readonly axis: Vec2 = new Vec2(0, 0);
    public readonly perp: Vec2 = new Vec2(0, 0);
    public s1: number = 0;
    public s2: number = 0;
    public a1: number = 0;
    public a2: number = 0;
    public readonly K: Mat22 = new Mat22();
    public readonly K3: Mat33 = new Mat33();
    public readonly K2: Mat22 = new Mat22();
    public translation: number = 0;
    public axialMass: number = 0;

    public readonly qA: Rot = new Rot();
    public readonly qB: Rot = new Rot();
    public readonly lalcA: Vec2 = new Vec2();
    public readonly lalcB: Vec2 = new Vec2();
    public readonly rA: Vec2 = new Vec2();
    public readonly rB: Vec2 = new Vec2();

    constructor(def: IPrismaticJointDef) {
      super(def);

      this.localAnchorA.copy(maybe(def.localAnchorA, Vec2.ZERO));
      this.localAnchorB.copy(maybe(def.localAnchorB, Vec2.ZERO));
      this.localXAxisA.copy(maybe(def.localAxisA, new Vec2(1, 0))).selfNormalize();
      Vec2.CrossOneV(this.localXAxisA, this.localYAxisA);
      this.referenceAngle = maybe(def.referenceAngle, 0);
      this.lowerTranslation = maybe(def.lowerTranslation, 0);
      this.upperTranslation = maybe(def.upperTranslation, 0);
      // Assert(this.lowerTranslation <= this.upperTranslation);
      this.maxMotorForce = maybe(def.maxMotorForce, 0);
      this.motorSpeed = maybe(def.motorSpeed, 0);
      this.enableLimit = maybe(def.enableLimit, false);
      this.enableMotor = maybe(def.enableMotor, false);
    }

    private static initVelocityConstraints_s_d = new Vec2();
    private static initVelocityConstraints_s_P = new Vec2();
    public initVelocityConstraints(data: SolverData): void {
      this.indexA = this.bodyA.islandIndex;
      this.indexB = this.bodyB.islandIndex;
      this.localCenterA.copy(this.bodyA.sweep.localCenter);
      this.localCenterB.copy(this.bodyB.sweep.localCenter);
      this.invMassA = this.bodyA.invMass;
      this.invMassB = this.bodyB.invMass;
      this.invIA = this.bodyA.invI;
      this.invIB = this.bodyB.invI;

      const cA: Vec2 = data.positions[this.indexA].c;
      const aA: number = data.positions[this.indexA].a;
      const vA: Vec2 = data.velocities[this.indexA].v;
      let wA: number = data.velocities[this.indexA].w;

      const cB: Vec2 = data.positions[this.indexB].c;
      const aB: number = data.positions[this.indexB].a;
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;

      const qA: Rot = this.qA.setAngle(aA), qB: Rot = this.qB.setAngle(aB);

      // Compute the effective masses.
      // Vec2 rA = Mul(qA, localAnchorA - localCenterA);
      Vec2.SubVV(this.localAnchorA, this.localCenterA, this.lalcA);
      const rA: Vec2 = Rot.mulRV(qA, this.lalcA, this.rA);
      // Vec2 rB = Mul(qB, localAnchorB - localCenterB);
      Vec2.SubVV(this.localAnchorB, this.localCenterB, this.lalcB);
      const rB: Vec2 = Rot.mulRV(qB, this.lalcB, this.rB);
      // Vec2 d = (cB - cA) + rB - rA;
      const d: Vec2 = Vec2.AddVV(
        Vec2.SubVV(cB, cA, Vec2.s_t0),
        Vec2.SubVV(rB, rA, Vec2.s_t1),
        PrismaticJoint.initVelocityConstraints_s_d);

      const mA: number = this.invMassA, mB: number = this.invMassB;
      const iA: number = this.invIA, iB: number = this.invIB;

      // Compute motor Jacobian and effective mass.
      {
        // axis = Mul(qA, localXAxisA);
        Rot.mulRV(qA, this.localXAxisA, this.axis);
        // a1 = Cross(d + rA, axis);
        this.a1 = Vec2.CrossVV(Vec2.AddVV(d, rA, Vec2.s_t0), this.axis);
        // a2 = Cross(rB, axis);
        this.a2 = Vec2.CrossVV(rB, this.axis);

        this.axialMass = mA + mB + iA * this.a1 * this.a1 + iB * this.a2 * this.a2;
        if (this.axialMass > 0) {
          this.axialMass = 1 / this.axialMass;
        }
      }

      // Prismatic constraint.
      {
        // perp = Mul(qA, localYAxisA);
        Rot.mulRV(qA, this.localYAxisA, this.perp);

        // s1 = Cross(d + rA, perp);
        this.s1 = Vec2.CrossVV(Vec2.AddVV(d, rA, Vec2.s_t0), this.perp);
        // s2 = Cross(rB, perp);
        this.s2 = Vec2.CrossVV(rB, this.perp);

        // float32 k11 = mA + mB + iA * s1 * s1 + iB * s2 * s2;
        this.K.ex.x = mA + mB + iA * this.s1 * this.s1 + iB * this.s2 * this.s2;
        // float32 k12 = iA * s1 + iB * s2;
        this.K.ex.y = iA * this.s1 + iB * this.s2;
        this.K.ey.x = this.K.ex.y;
        // float32 k22 = iA + iB;
        this.K.ey.y = iA + iB;
        if (this.K.ey.y === 0) {
          // For bodies with fixed rotation.
          this.K.ey.y = 1;
        }

        // K.ex.Set(k11, k12);
        // K.ey.Set(k12, k22);
      }

      // Compute motor and limit terms.
      if (this.enableLimit) {
        this.translation = Vec2.DotVV(this.axis, d);
      } else {
        this.lowerImpulse = 0.0;
        this.upperImpulse = 0.0;
      }

      if (!this.enableMotor) {
        this.motorImpulse = 0;
      }

      if (data.step.warmStarting) {
        // Account for variable time step.
        // impulse *= data.step.dtRatio;
        this.impulse.selfMul(data.step.dtRatio);
        this.motorImpulse *= data.step.dtRatio;
        this.lowerImpulse *= data.step.dtRatio;
        this.upperImpulse *= data.step.dtRatio;

        const axialImpulse: number = this.motorImpulse + this.lowerImpulse - this.upperImpulse;
        // Vec2 P = impulse.x * perp + axialImpulse * axis;
        const P: Vec2 = Vec2.AddVV(
          Vec2.MulSV(this.impulse.x, this.perp, Vec2.s_t0),
          Vec2.MulSV(axialImpulse, this.axis, Vec2.s_t1),
          PrismaticJoint.initVelocityConstraints_s_P);
        // float LA = impulse.x * s1 + impulse.y + axialImpulse * a1;
        const LA = this.impulse.x * this.s1 + this.impulse.y + axialImpulse * this.a1;
        // float LB = impulse.x * s2 + impulse.y + axialImpulse * a2;
        const LB = this.impulse.x * this.s2 + this.impulse.y + axialImpulse * this.a2;

        // vA -= mA * P;
        vA.selfMulSub(mA, P);
        wA -= iA * LA;

        // vB += mB * P;
        vB.selfMulAdd(mB, P);
        wB += iB * LB;
      } else {
        this.impulse.setZero();
        this.motorImpulse = 0.0;
        this.lowerImpulse = 0.0;
        this.upperImpulse = 0.0;
      }

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    private static solveVelocityConstraints_s_P = new Vec2();
    // private static solveVelocityConstraints_s_f2r = new Vec2();
    // private static solveVelocityConstraints_s_f1 = new Vec3();
    // private static solveVelocityConstraints_s_df3 = new Vec3();
    private static solveVelocityConstraints_s_df = new Vec2();
    public solveVelocityConstraints(data: SolverData): void {
      const vA: Vec2 = data.velocities[this.indexA].v;
      let wA: number = data.velocities[this.indexA].w;
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;

      const mA: number = this.invMassA, mB: number = this.invMassB;
      const iA: number = this.invIA, iB: number = this.invIB;

      // Solve linear motor constraint.
      if (this.enableMotor) {
        // float32 Cdot = Dot(axis, vB - vA) + a2 * wB - a1 * wA;
        const Cdot: number = Vec2.DotVV(this.axis, Vec2.SubVV(vB, vA, Vec2.s_t0)) + this.a2 * wB - this.a1 * wA;
        let impulse = this.axialMass * (this.motorSpeed - Cdot);
        const oldImpulse = this.motorImpulse;
        const maxImpulse = data.step.dt * this.maxMotorForce;
        this.motorImpulse = clamp(this.motorImpulse + impulse, (-maxImpulse), maxImpulse);
        impulse = this.motorImpulse - oldImpulse;

        // Vec2 P = impulse * axis;
        const P: Vec2 = Vec2.MulSV(impulse, this.axis, PrismaticJoint.solveVelocityConstraints_s_P);
        const LA = impulse * this.a1;
        const LB = impulse * this.a2;

        // vA -= mA * P;
        vA.selfMulSub(mA, P);
        wA -= iA * LA;
        // vB += mB * P;
        vB.selfMulAdd(mB, P);
        wB += iB * LB;
      }

      if (this.enableLimit) {
        // Lower limit
        {
          const C: number = this.translation - this.lowerTranslation;
          const Cdot: number = Vec2.DotVV(this.axis, Vec2.SubVV(vB, vA, Vec2.s_t0)) + this.a2 * wB - this.a1 * wA;
          let impulse: number = -this.axialMass * (Cdot + Max(C, 0.0) * data.step.inv_dt);
          const oldImpulse: number = this.lowerImpulse;
          this.lowerImpulse = Max(this.lowerImpulse + impulse, 0.0);
          impulse = this.lowerImpulse - oldImpulse;

          // Vec2 P = impulse * this.axis;
          const P: Vec2 = Vec2.MulSV(impulse, this.axis, PrismaticJoint.solveVelocityConstraints_s_P);
          const LA: number = impulse * this.a1;
          const LB: number = impulse * this.a2;

          // vA -= mA * P;
          vA.selfMulSub(mA, P);
          wA -= iA * LA;
          // vB += mB * P;
          vB.selfMulAdd(mB, P);
          wB += iB * LB;
        }

        // Upper limit
        // Note: signs are flipped to keep C positive when the constraint is satisfied.
        // This also keeps the impulse positive when the limit is active.
        {
          const C: number = this.upperTranslation - this.translation;
          const Cdot: number = Vec2.DotVV(this.axis, Vec2.SubVV(vA, vB, Vec2.s_t0)) + this.a1 * wA - this.a2 * wB;
          let impulse: number = -this.axialMass * (Cdot + Max(C, 0.0) * data.step.inv_dt);
          const oldImpulse: number = this.upperImpulse;
          this.upperImpulse = Max(this.upperImpulse + impulse, 0.0);
          impulse = this.upperImpulse - oldImpulse;

          // Vec2 P = impulse * this.axis;
          const P: Vec2 = Vec2.MulSV(impulse, this.axis, PrismaticJoint.solveVelocityConstraints_s_P);
          const LA: number = impulse * this.a1;
          const LB: number = impulse * this.a2;

          // vA += mA * P;
          vA.selfMulAdd(mA, P);
          wA += iA * LA;
          // vB -= mB * P;
          vB.selfMulSub(mB, P);
          wB -= iB * LB;
        }
      }

      // Solve the prismatic constraint in block form.
      {
        // Vec2 Cdot;
        // Cdot.x = Dot(perp, vB - vA) + s2 * wB - s1 * wA;
        const Cdot_x: number = Vec2.DotVV(this.perp, Vec2.SubVV(vB, vA, Vec2.s_t0)) + this.s2 * wB - this.s1 * wA;
        // Cdot.y = wB - wA;
        const Cdot_y = wB - wA;

        // Vec2 df = K.Solve(-Cdot);
        const df = this.K.solve(-Cdot_x, -Cdot_y, PrismaticJoint.solveVelocityConstraints_s_df);
        // impulse += df;
        this.impulse.selfAdd(df);

        // Vec2 P = df.x * perp;
        const P: Vec2 = Vec2.MulSV(df.x, this.perp, PrismaticJoint.solveVelocityConstraints_s_P);
        // float32 LA = df.x * s1 + df.y;
        const LA = df.x * this.s1 + df.y;
        // float32 LB = df.x * s2 + df.y;
        const LB = df.x * this.s2 + df.y;

        // vA -= mA * P;
        vA.selfMulSub(mA, P);
        wA -= iA * LA;

        // vB += mB * P;
        vB.selfMulAdd(mB, P);
        wB += iB * LB;
      }

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    // A velocity based solver computes reaction forces(impulses) using the velocity constraint solver.Under this context,
    // the position solver is not there to resolve forces.It is only there to cope with integration error.
    //
    // Therefore, the pseudo impulses in the position solver do not have any physical meaning.Thus it is okay if they suck.
    //
    // We could take the active state from the velocity solver.However, the joint might push past the limit when the velocity
    // solver indicates the limit is inactive.
    private static solvePositionConstraints_s_d = new Vec2();
    private static solvePositionConstraints_s_impulse = new Vec3();
    private static solvePositionConstraints_s_impulse1 = new Vec2();
    private static solvePositionConstraints_s_P = new Vec2();
    public solvePositionConstraints(data: SolverData): boolean {
      const cA: Vec2 = data.positions[this.indexA].c;
      let aA: number = data.positions[this.indexA].a;
      const cB: Vec2 = data.positions[this.indexB].c;
      let aB: number = data.positions[this.indexB].a;

      const qA: Rot = this.qA.setAngle(aA), qB: Rot = this.qB.setAngle(aB);

      const mA: number = this.invMassA, mB: number = this.invMassB;
      const iA: number = this.invIA, iB: number = this.invIB;

      // Vec2 rA = Mul(qA, localAnchorA - localCenterA);
      const rA: Vec2 = Rot.mulRV(qA, this.lalcA, this.rA);
      // Vec2 rB = Mul(qB, localAnchorB - localCenterB);
      const rB: Vec2 = Rot.mulRV(qB, this.lalcB, this.rB);
      // Vec2 d = cB + rB - cA - rA;
      const d: Vec2 = Vec2.SubVV(
        Vec2.AddVV(cB, rB, Vec2.s_t0),
        Vec2.AddVV(cA, rA, Vec2.s_t1),
        PrismaticJoint.solvePositionConstraints_s_d);

      // Vec2 axis = Mul(qA, localXAxisA);
      const axis: Vec2 = Rot.mulRV(qA, this.localXAxisA, this.axis);
      // float32 a1 = Cross(d + rA, axis);
      const a1 = Vec2.CrossVV(Vec2.AddVV(d, rA, Vec2.s_t0), axis);
      // float32 a2 = Cross(rB, axis);
      const a2 = Vec2.CrossVV(rB, axis);
      // Vec2 perp = Mul(qA, localYAxisA);
      const perp: Vec2 = Rot.mulRV(qA, this.localYAxisA, this.perp);

      // float32 s1 = Cross(d + rA, perp);
      const s1 = Vec2.CrossVV(Vec2.AddVV(d, rA, Vec2.s_t0), perp);
      // float32 s2 = Cross(rB, perp);
      const s2 = Vec2.CrossVV(rB, perp);

      // Vec3 impulse;
      let impulse = PrismaticJoint.solvePositionConstraints_s_impulse;
      // Vec2 C1;
      // C1.x = Dot(perp, d);
      const C1_x: number = Vec2.DotVV(perp, d);
      // C1.y = aB - aA - referenceAngle;
      const C1_y = aB - aA - this.referenceAngle;

      let linearError = Abs(C1_x);
      const angularError = Abs(C1_y);

      let active = false;
      let C2: number = 0;
      if (this.enableLimit) {
        // float32 translation = Dot(axis, d);
        const translation: number = Vec2.DotVV(axis, d);
        if (Abs(this.upperTranslation - this.lowerTranslation) < 2 * linearSlop) {
          C2 = translation;
          linearError = Max(linearError, Abs(translation));
          active = true;
        } else if (translation <= this.lowerTranslation) {
          C2 = Min(translation - this.lowerTranslation, 0.0);
          linearError = Max(linearError, this.lowerTranslation - translation);
          active = true;
        } else if (translation >= this.upperTranslation) {
          C2 = Max(translation - this.upperTranslation, 0.0);
          linearError = Max(linearError, translation - this.upperTranslation);
          active = true;
        }
      }

      if (active) {
        // float32 k11 = mA + mB + iA * s1 * s1 + iB * s2 * s2;
        const k11 = mA + mB + iA * s1 * s1 + iB * s2 * s2;
        // float32 k12 = iA * s1 + iB * s2;
        const k12 = iA * s1 + iB * s2;
        // float32 k13 = iA * s1 * a1 + iB * s2 * a2;
        const k13 = iA * s1 * a1 + iB * s2 * a2;
        // float32 k22 = iA + iB;
        let k22 = iA + iB;
        if (k22 === 0) {
          // For fixed rotation
          k22 = 1;
        }
        // float32 k23 = iA * a1 + iB * a2;
        const k23 = iA * a1 + iB * a2;
        // float32 k33 = mA + mB + iA * a1 * a1 + iB * a2 * a2;
        const k33 = mA + mB + iA * a1 * a1 + iB * a2 * a2;

        // Mat33 K;
        const K = this.K3;
        // K.ex.Set(k11, k12, k13);
        K.ex.setXYZ(k11, k12, k13);
        // K.ey.Set(k12, k22, k23);
        K.ey.setXYZ(k12, k22, k23);
        // K.ez.Set(k13, k23, k33);
        K.ez.setXYZ(k13, k23, k33);

        // Vec3 C;
        // C.x = C1.x;
        // C.y = C1.y;
        // C.z = C2;

        // impulse = K.Solve33(-C);
        impulse = K.solve33((-C1_x), (-C1_y), (-C2), impulse);
      } else {
        // float32 k11 = mA + mB + iA * s1 * s1 + iB * s2 * s2;
        const k11 = mA + mB + iA * s1 * s1 + iB * s2 * s2;
        // float32 k12 = iA * s1 + iB * s2;
        const k12 = iA * s1 + iB * s2;
        // float32 k22 = iA + iB;
        let k22 = iA + iB;
        if (k22 === 0) {
          k22 = 1;
        }

        // Mat22 K;
        const K2 = this.K2;
        // K.ex.Set(k11, k12);
        K2.ex.set(k11, k12);
        // K.ey.Set(k12, k22);
        K2.ey.set(k12, k22);

        // Vec2 impulse1 = K.Solve(-C1);
        const impulse1 = K2.solve((-C1_x), (-C1_y), PrismaticJoint.solvePositionConstraints_s_impulse1);
        impulse.x = impulse1.x;
        impulse.y = impulse1.y;
        impulse.z = 0;
      }

      // Vec2 P = impulse.x * perp + impulse.z * axis;
      const P: Vec2 = Vec2.AddVV(
        Vec2.MulSV(impulse.x, perp, Vec2.s_t0),
        Vec2.MulSV(impulse.z, axis, Vec2.s_t1),
        PrismaticJoint.solvePositionConstraints_s_P);
      // float32 LA = impulse.x * s1 + impulse.y + impulse.z * a1;
      const LA = impulse.x * s1 + impulse.y + impulse.z * a1;
      // float32 LB = impulse.x * s2 + impulse.y + impulse.z * a2;
      const LB = impulse.x * s2 + impulse.y + impulse.z * a2;

      // cA -= mA * P;
      cA.selfMulSub(mA, P);
      aA -= iA * LA;
      // cB += mB * P;
      cB.selfMulAdd(mB, P);
      aB += iB * LB;

      // data.positions[this.indexA].c = cA;
      data.positions[this.indexA].a = aA;
      // data.positions[this.indexB].c = cB;
      data.positions[this.indexB].a = aB;

      return linearError <= linearSlop && angularError <= angularSlop;
    }

    public getAnchorA<T extends XY>(out: T): T {
      return this.bodyA.getWorldPoint(this.localAnchorA, out);
    }

    public getAnchorB<T extends XY>(out: T): T {
      return this.bodyB.getWorldPoint(this.localAnchorB, out);
    }

    public getReactionForce<T extends XY>(inv_dt: number, out: T): T {
      out.x = inv_dt * (this.impulse.x * this.perp.x + (this.motorImpulse + this.lowerImpulse - this.upperImpulse) * this.axis.x);
      out.y = inv_dt * (this.impulse.y * this.perp.y + (this.motorImpulse + this.lowerImpulse - this.upperImpulse) * this.axis.y);
      return out;
    }

    public getReactionTorque(inv_dt: number): number {
      return inv_dt * this.impulse.y;
    }
    
    private static getJointTranslation_s_pA = new Vec2();
    private static getJointTranslation_s_pB = new Vec2();
    private static getJointTranslation_s_d = new Vec2();
    private static getJointTranslation_s_axis = new Vec2();
    public getJointTranslation(): number {
      // Vec2 pA = bodyA.GetWorldPoint(localAnchorA);
      const pA = this.bodyA.getWorldPoint(this.localAnchorA, PrismaticJoint.getJointTranslation_s_pA);
      // Vec2 pB = bodyB.GetWorldPoint(localAnchorB);
      const pB = this.bodyB.getWorldPoint(this.localAnchorB, PrismaticJoint.getJointTranslation_s_pB);
      // Vec2 d = pB - pA;
      const d: Vec2 = Vec2.SubVV(pB, pA, PrismaticJoint.getJointTranslation_s_d);
      // Vec2 axis = bodyA.GetWorldVector(localXAxisA);
      const axis = this.bodyA.getWorldVector(this.localXAxisA, PrismaticJoint.getJointTranslation_s_axis);

      // float32 translation = Dot(d, axis);
      const translation: number = Vec2.DotVV(d, axis);
      return translation;
    }

    public getJointSpeed(): number {
      const bA: Body = this.bodyA;
      const bB: Body = this.bodyB;

      // Vec2 rA = Mul(bA->xf.q, localAnchorA - bA->sweep.localCenter);
      Vec2.SubVV(this.localAnchorA, bA.sweep.localCenter, this.lalcA);
      const rA: Vec2 = Rot.mulRV(bA.xf.q, this.lalcA, this.rA);
      // Vec2 rB = Mul(bB->xf.q, localAnchorB - bB->sweep.localCenter);
      Vec2.SubVV(this.localAnchorB, bB.sweep.localCenter, this.lalcB);
      const rB: Vec2 = Rot.mulRV(bB.xf.q, this.lalcB, this.rB);
      // Vec2 pA = bA->sweep.c + rA;
      const pA: Vec2 = Vec2.AddVV(bA.sweep.c, rA, Vec2.s_t0); // pA uses s_t0
      // Vec2 pB = bB->sweep.c + rB;
      const pB: Vec2 = Vec2.AddVV(bB.sweep.c, rB, Vec2.s_t1); // pB uses s_t1
      // Vec2 d = pB - pA;
      const d: Vec2 = Vec2.SubVV(pB, pA, Vec2.s_t2); // d uses s_t2
      // Vec2 axis = Mul(bA.xf.q, localXAxisA);
      const axis = bA.getWorldVector(this.localXAxisA, this.axis);

      const vA = bA.linearVelocity;
      const vB = bB.linearVelocity;
      const wA = bA.angularVelocity;
      const wB = bB.angularVelocity;

      // float32 speed = Dot(d, Cross(wA, axis)) + Dot(axis, vB + Cross(wB, rB) - vA - Cross(wA, rA));
      const speed =
        Vec2.DotVV(d, Vec2.CrossSV(wA, axis, Vec2.s_t0)) +
        Vec2.DotVV(
          axis,
          Vec2.SubVV(
            Vec2.AddVCrossSV(vB, wB, rB, Vec2.s_t0),
            Vec2.AddVCrossSV(vA, wA, rA, Vec2.s_t1),
            Vec2.s_t0));
      return speed;
    }

    public isLimitEnabled() {
      return this.enableLimit;
    }

    public setEnableLimit(flag: boolean) {
      if (flag !== this.enableLimit) {
        this.bodyA.setAwake(true);
        this.bodyB.setAwake(true);
        this.enableLimit = flag;
        this.lowerImpulse = 0.0;
        this.upperImpulse = 0.0;
      }
    }

    public getLowerLimit() {
      return this.lowerTranslation;
    }

    public getUpperLimit() {
      return this.upperTranslation;
    }

    public setLimits(lower: number, upper: number): void {
      if (lower !== this.lowerTranslation || upper !== this.upperTranslation) {
        this.bodyA.setAwake(true);
        this.bodyB.setAwake(true);
        this.lowerTranslation = lower;
        this.upperTranslation = upper;
        this.lowerImpulse = 0.0;
        this.upperImpulse = 0.0;
      }
    }

    public isMotorEnabled(): boolean {
      return this.enableMotor;
    }

    public setEnableMotor(flag: boolean): void {
      if (flag !== this.enableMotor) {
        this.bodyA.setAwake(true);
        this.bodyB.setAwake(true);
        this.enableMotor = flag;
      }
    }

    public setMotorSpeed(speed: number): void {
      if (speed !== this.motorSpeed) {
        this.bodyA.setAwake(true);
        this.bodyB.setAwake(true);
        this.motorSpeed = speed;
      }
    }

    public getMotorSpeed() {
      return this.motorSpeed;
    }

    public setMaxMotorForce(force: number): void {
      if (force !== this.maxMotorForce) {
        this.bodyA.setAwake(true);
        this.bodyB.setAwake(true);
        this.maxMotorForce = force;
      }
    }

    public getMaxMotorForce(): number { return this.maxMotorForce; }

    public getMotorForce(inv_dt: number): number {
      return inv_dt * this.motorImpulse;
    }

    public dump(log: (format: string, ...args: any[]) => void) {
      const indexA = this.bodyA.islandIndex;
      const indexB = this.bodyB.islandIndex;

      log("  const jd: PrismaticJointDef = new PrismaticJointDef();\n");
      log("  jd.bodyA = bodies[%d];\n", indexA);
      log("  jd.bodyB = bodies[%d];\n", indexB);
      log("  jd.collideConnected = %s;\n", (this.collideConnected) ? ("true") : ("false"));
      log("  jd.localAnchorA.Set(%.15f, %.15f);\n", this.localAnchorA.x, this.localAnchorA.y);
      log("  jd.localAnchorB.Set(%.15f, %.15f);\n", this.localAnchorB.x, this.localAnchorB.y);
      log("  jd.localAxisA.Set(%.15f, %.15f);\n", this.localXAxisA.x, this.localXAxisA.y);
      log("  jd.referenceAngle = %.15f;\n", this.referenceAngle);
      log("  jd.enableLimit = %s;\n", (this.enableLimit) ? ("true") : ("false"));
      log("  jd.lowerTranslation = %.15f;\n", this.lowerTranslation);
      log("  jd.upperTranslation = %.15f;\n", this.upperTranslation);
      log("  jd.enableMotor = %s;\n", (this.enableMotor) ? ("true") : ("false"));
      log("  jd.motorSpeed = %.15f;\n", this.motorSpeed);
      log("  jd.maxMotorForce = %.15f;\n", this.maxMotorForce);
      log("  joints[%d] = this.world.CreateJoint(jd);\n", this.index);
    }

    private static draw_s_pA = new Vec2();
    private static draw_s_pB = new Vec2();
    private static draw_s_axis = new Vec2();
    private static draw_s_c1 = new Color(0.7, 0.7, 0.7);
    private static draw_s_c2 = new Color(0.3, 0.9, 0.3);
    private static draw_s_c3 = new Color(0.9, 0.3, 0.3);
    private static draw_s_c4 = new Color(0.3, 0.3, 0.9);
    private static draw_s_c5 = new Color(0.4, 0.4, 0.4);
    private static draw_s_lower = new Vec2();
    private static draw_s_upper = new Vec2();
    private static draw_s_perp = new Vec2();
    public draw(draw: Draw): void {
      const xfA: Transform = this.bodyA.getTransform();
      const xfB: Transform = this.bodyB.getTransform();
      const pA = Transform.mulXV(xfA, this.localAnchorA, PrismaticJoint.draw_s_pA);
      const pB = Transform.mulXV(xfB, this.localAnchorB, PrismaticJoint.draw_s_pB);

      // Vec2 axis = Mul(xfA.q, localXAxisA);
      const axis: Vec2 = Rot.mulRV(xfA.q, this.localXAxisA, PrismaticJoint.draw_s_axis);

      const c1 = PrismaticJoint.draw_s_c1; // Color c1(0.7f, 0.7f, 0.7f);
      const c2 = PrismaticJoint.draw_s_c2; // Color c2(0.3f, 0.9f, 0.3f);
      const c3 = PrismaticJoint.draw_s_c3; // Color c3(0.9f, 0.3f, 0.3f);
      const c4 = PrismaticJoint.draw_s_c4; // Color c4(0.3f, 0.3f, 0.9f);
      const c5 = PrismaticJoint.draw_s_c5; // Color c5(0.4f, 0.4f, 0.4f);

      draw.drawSegment(pA, pB, c5);

      if (this.enableLimit) {
        // Vec2 lower = pA + lowerTranslation * axis;
        const lower = Vec2.AddVMulSV(pA, this.lowerTranslation, axis, PrismaticJoint.draw_s_lower);
        // Vec2 upper = pA + upperTranslation * axis;
        const upper = Vec2.AddVMulSV(pA, this.upperTranslation, axis, PrismaticJoint.draw_s_upper);
        // Vec2 perp = Mul(xfA.q, localYAxisA);
        const perp = Rot.mulRV(xfA.q, this.localYAxisA, PrismaticJoint.draw_s_perp);
        draw.drawSegment(lower, upper, c1);
        // draw.DrawSegment(lower - 0.5 * perp, lower + 0.5 * perp, c2);
        draw.drawSegment(Vec2.AddVMulSV(lower, -0.5, perp, Vec2.s_t0), Vec2.AddVMulSV(lower, 0.5, perp, Vec2.s_t1), c2);
        // draw.DrawSegment(upper - 0.5 * perp, upper + 0.5 * perp, c3);
        draw.drawSegment(Vec2.AddVMulSV(upper, -0.5, perp, Vec2.s_t0), Vec2.AddVMulSV(upper, 0.5, perp, Vec2.s_t1), c3);
      } else {
        // draw.DrawSegment(pA - 1.0 * axis, pA + 1.0 * axis, c1);
        draw.drawSegment(Vec2.AddVMulSV(pA, -1.0, axis, Vec2.s_t0), Vec2.AddVMulSV(pA, 1.0, axis, Vec2.s_t1), c1);
      }

      draw.drawPoint(pA, 5.0, c1);
      draw.drawPoint(pB, 5.0, c4);
    }
  }
}

