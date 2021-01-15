/*
* Copyright (c) 2006-2012 Erin Catto http://www.box2d.org
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

// Point-to-point constraint
// Cdot = v2 - v1
//      = v2 + cross(w2, r2) - v1 - cross(w1, r1)
// J = [-I -r1_skew I r2_skew ]
// Identity used:
// w k % (rx i + ry j) = w * (-ry i + rx j)
//
// r1 = offset - c1
// r2 = -c2

// Angle constraint
// Cdot = w2 - w1
// J = [0 0 -1 0 0 1]
// K = invI1 + invI2
namespace b2 {
  export interface IMotorJointDef extends IJointDef {
    linearOffset?: XY;

    angularOffset?: number;

    maxForce?: number;

    maxTorque?: number;

    correctionFactor?: number;
  }

  export class MotorJointDef extends JointDef implements IMotorJointDef {
    public readonly linearOffset: Vec2 = new Vec2(0, 0);

    public angularOffset: number = 0;

    public maxForce: number = 1;

    public maxTorque: number = 1;

    public correctionFactor: number = 0.3;

    constructor() {
      super(JointType.e_motorJoint);
    }

    public Initialize(bA: Body, bB: Body): void {
      this.bodyA = bA;
      this.bodyB = bB;
      // Vec2 xB = bodyB->GetPosition();
      // linearOffset = bodyA->GetLocalPoint(xB);
      this.bodyA.GetLocalPoint(this.bodyB.GetPosition(), this.linearOffset);

      const angleA: number = this.bodyA.GetAngle();
      const angleB: number = this.bodyB.GetAngle();
      this.angularOffset = angleB - angleA;
    }
  }

  export class MotorJoint extends Joint {
    // Solver shared
    public readonly linearOffset: Vec2 = new Vec2();
    public angularOffset: number = 0;
    public readonly linearImpulse: Vec2 = new Vec2();
    public angularImpulse: number = 0;
    public maxForce: number = 0;
    public maxTorque: number = 0;
    public correctionFactor: number = 0.3;

    // Solver temp
    public indexA: number = 0;
    public indexB: number = 0;
    public readonly rA: Vec2 = new Vec2();
    public readonly rB: Vec2 = new Vec2();
    public readonly localCenterA: Vec2 = new Vec2();
    public readonly localCenterB: Vec2 = new Vec2();
    public readonly linearError: Vec2 = new Vec2();
    public angularError: number = 0;
    public invMassA: number = 0;
    public invMassB: number = 0;
    public invIA: number = 0;
    public invIB: number = 0;
    public readonly linearMass: Mat22 = new Mat22();
    public angularMass: number = 0;

    public readonly qA: Rot = new Rot();
    public readonly qB: Rot = new Rot();
    public readonly K: Mat22 = new Mat22();

    constructor(def: IMotorJointDef) {
      super(def);

      this.linearOffset.Copy(Maybe(def.linearOffset, Vec2.ZERO));
      this.linearImpulse.SetZero();
      this.maxForce = Maybe(def.maxForce, 0);
      this.maxTorque = Maybe(def.maxTorque, 0);
      this.correctionFactor = Maybe(def.correctionFactor, 0.3);
    }

    public GetAnchorA<T extends XY>(out: T): T {
      const pos: Vec2 = this.bodyA.GetPosition();
      out.x = pos.x;
      out.y = pos.y;
      return out;
    }
    public GetAnchorB<T extends XY>(out: T): T {
      const pos: Vec2 = this.bodyB.GetPosition();
      out.x = pos.x;
      out.y = pos.y;
      return out;
    }

    public GetReactionForce<T extends XY>(inv_dt: number, out: T): T {
      // return inv_dt * linearImpulse;
      return Vec2.MulSV(inv_dt, this.linearImpulse, out);
    }

    public GetReactionTorque(inv_dt: number): number {
      return inv_dt * this.angularImpulse;
    }

    public SetLinearOffset(linearOffset: Vec2): void {
      if (!Vec2.IsEqualToV(linearOffset, this.linearOffset)) {
        this.bodyA.SetAwake(true);
        this.bodyB.SetAwake(true);
        this.linearOffset.Copy(linearOffset);
      }
    }
    public GetLinearOffset() {
      return this.linearOffset;
    }

    public SetAngularOffset(angularOffset: number): void {
      if (angularOffset !== this.angularOffset) {
        this.bodyA.SetAwake(true);
        this.bodyB.SetAwake(true);
        this.angularOffset = angularOffset;
      }
    }
    public GetAngularOffset() {
      return this.angularOffset;
    }

    public SetMaxForce(force: number): void {
      // DEBUG: Assert(IsValid(force) && force >= 0);
      this.maxForce = force;
    }

    public GetMaxForce() {
      return this.maxForce;
    }

    public SetMaxTorque(torque: number): void {
      // DEBUG: Assert(IsValid(torque) && torque >= 0);
      this.maxTorque = torque;
    }

    public GetMaxTorque() {
      return this.maxTorque;
    }

    public InitVelocityConstraints(data: SolverData): void {
      this.indexA = this.bodyA.islandIndex;
      this.indexB = this.bodyB.islandIndex;
      this.localCenterA.Copy(this.bodyA.sweep.localCenter);
      this.localCenterB.Copy(this.bodyB.sweep.localCenter);
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

      const qA: Rot = this.qA.SetAngle(aA), qB: Rot = this.qB.SetAngle(aB);

      // Compute the effective mass matrix.
      // this.rA = Mul(qA, linearOffset - this.localCenterA);
      const rA: Vec2 = Rot.MulRV(qA, Vec2.SubVV(this.linearOffset, this.localCenterA, Vec2.s_t0), this.rA);
      // this.rB = Mul(qB, -this.localCenterB);
      const rB: Vec2 = Rot.MulRV(qB, Vec2.NegV(this.localCenterB, Vec2.s_t0), this.rB);

      // J = [-I -r1_skew I r2_skew]
      // r_skew = [-ry; rx]

      // Matlab
      // K = [ mA+r1y^2*iA+mB+r2y^2*iB,  -r1y*iA*r1x-r2y*iB*r2x,          -r1y*iA-r2y*iB]
      //     [  -r1y*iA*r1x-r2y*iB*r2x, mA+r1x^2*iA+mB+r2x^2*iB,           r1x*iA+r2x*iB]
      //     [          -r1y*iA-r2y*iB,           r1x*iA+r2x*iB,                   iA+iB]

      const mA: number = this.invMassA, mB: number = this.invMassB;
      const iA: number = this.invIA, iB: number = this.invIB;

      // Upper 2 by 2 of K for point to point
      const K: Mat22 = this.K;
      K.ex.x = mA + mB + iA * rA.y * rA.y + iB * rB.y * rB.y;
      K.ex.y = -iA * rA.x * rA.y - iB * rB.x * rB.y;
      K.ey.x = K.ex.y;
      K.ey.y = mA + mB + iA * rA.x * rA.x + iB * rB.x * rB.x;

      // this.linearMass = K.GetInverse();
      K.GetInverse(this.linearMass);

      this.angularMass = iA + iB;
      if (this.angularMass > 0) {
        this.angularMass = 1 / this.angularMass;
      }

      // this.linearError = cB + rB - cA - rA;
      Vec2.SubVV(
        Vec2.AddVV(cB, rB, Vec2.s_t0),
        Vec2.AddVV(cA, rA, Vec2.s_t1),
        this.linearError);
      this.angularError = aB - aA - this.angularOffset;

      if (data.step.warmStarting) {
        // Scale impulses to support a variable time step.
        // this.linearImpulse *= data.step.dtRatio;
        this.linearImpulse.SelfMul(data.step.dtRatio);
        this.angularImpulse *= data.step.dtRatio;

        // Vec2 P(this.linearImpulse.x, this.linearImpulse.y);
        const P: Vec2 = this.linearImpulse;
        // vA -= mA * P;
        vA.SelfMulSub(mA, P);
        wA -= iA * (Vec2.CrossVV(rA, P) + this.angularImpulse);
        // vB += mB * P;
        vB.SelfMulAdd(mB, P);
        wB += iB * (Vec2.CrossVV(rB, P) + this.angularImpulse);
      } else {
        this.linearImpulse.SetZero();
        this.angularImpulse = 0;
      }

      // data.velocities[this.indexA].v = vA; // vA is a reference
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB; // vB is a reference
      data.velocities[this.indexB].w = wB;
    }

    private static SolveVelocityConstraints_s_Cdot_v2 = new Vec2();
    private static SolveVelocityConstraints_s_impulse_v2 = new Vec2();
    private static SolveVelocityConstraints_s_oldImpulse_v2 = new Vec2();
    public SolveVelocityConstraints(data: SolverData): void {
      const vA: Vec2 = data.velocities[this.indexA].v;
      let wA: number = data.velocities[this.indexA].w;
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;

      const mA: number = this.invMassA, mB: number = this.invMassB;
      const iA: number = this.invIA, iB: number = this.invIB;

      const h: number = data.step.dt;
      const inv_h: number = data.step.inv_dt;

      // Solve angular friction
      {
        const Cdot: number = wB - wA + inv_h * this.correctionFactor * this.angularError;
        let impulse: number = -this.angularMass * Cdot;

        const oldImpulse: number = this.angularImpulse;
        const maxImpulse: number = h * this.maxTorque;
        this.angularImpulse = Clamp(this.angularImpulse + impulse, -maxImpulse, maxImpulse);
        impulse = this.angularImpulse - oldImpulse;

        wA -= iA * impulse;
        wB += iB * impulse;
      }

      // Solve linear friction
      {
        const rA = this.rA;
        const rB = this.rB;

        // Vec2 Cdot = vB + Vec2.CrossSV(wB, rB) - vA - Vec2.CrossSV(wA, rA) + inv_h * this.correctionFactor * this.linearError;
        const Cdot_v2 =
          Vec2.AddVV(
            Vec2.SubVV(
              Vec2.AddVV(vB, Vec2.CrossSV(wB, rB, Vec2.s_t0), Vec2.s_t0),
              Vec2.AddVV(vA, Vec2.CrossSV(wA, rA, Vec2.s_t1), Vec2.s_t1), Vec2.s_t2),
            Vec2.MulSV(inv_h * this.correctionFactor, this.linearError, Vec2.s_t3),
            MotorJoint.SolveVelocityConstraints_s_Cdot_v2);

        // Vec2 impulse = -Mul(this.linearMass, Cdot);
        const impulse_v2: Vec2 = Mat22.MulMV(this.linearMass, Cdot_v2, MotorJoint.SolveVelocityConstraints_s_impulse_v2).SelfNeg();
        // Vec2 oldImpulse = this.linearImpulse;
        const oldImpulse_v2 = MotorJoint.SolveVelocityConstraints_s_oldImpulse_v2.Copy(this.linearImpulse);
        // this.linearImpulse += impulse;
        this.linearImpulse.SelfAdd(impulse_v2);

        const maxImpulse: number = h * this.maxForce;

        if (this.linearImpulse.LengthSquared() > maxImpulse * maxImpulse) {
          this.linearImpulse.Normalize();
          // this.linearImpulse *= maxImpulse;
          this.linearImpulse.SelfMul(maxImpulse);
        }

        // impulse = this.linearImpulse - oldImpulse;
        Vec2.SubVV(this.linearImpulse, oldImpulse_v2, impulse_v2);

        // vA -= mA * impulse;
        vA.SelfMulSub(mA, impulse_v2);
        // wA -= iA * Vec2.CrossVV(rA, impulse);
        wA -= iA * Vec2.CrossVV(rA, impulse_v2);

        // vB += mB * impulse;
        vB.SelfMulAdd(mB, impulse_v2);
        // wB += iB * Vec2.CrossVV(rB, impulse);
        wB += iB * Vec2.CrossVV(rB, impulse_v2);
      }

      // data.velocities[this.indexA].v = vA; // vA is a reference
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB; // vB is a reference
      data.velocities[this.indexB].w = wB;
    }

    public SolvePositionConstraints(data: SolverData): boolean {
      return true;
    }

    public Dump(log: (format: string, ...args: any[]) => void) {
      const indexA = this.bodyA.islandIndex;
      const indexB = this.bodyB.islandIndex;

      log("  const jd: MotorJointDef = new MotorJointDef();\n");

      log("  jd.bodyA = bodies[%d];\n", indexA);
      log("  jd.bodyB = bodies[%d];\n", indexB);
      log("  jd.collideConnected = %s;\n", (this.collideConnected) ? ("true") : ("false"));

      log("  jd.linearOffset.Set(%.15f, %.15f);\n", this.linearOffset.x, this.linearOffset.y);
      log("  jd.angularOffset = %.15f;\n", this.angularOffset);
      log("  jd.maxForce = %.15f;\n", this.maxForce);
      log("  jd.maxTorque = %.15f;\n", this.maxTorque);
      log("  jd.correctionFactor = %.15f;\n", this.correctionFactor);
      log("  joints[%d] = this.world.CreateJoint(jd);\n", this.index);
    }
  }

}
