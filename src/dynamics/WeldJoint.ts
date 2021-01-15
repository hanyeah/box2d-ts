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
  export interface IWeldJointDef extends IJointDef {
    localAnchorA?: XY;

    localAnchorB?: XY;

    referenceAngle?: number;

    stiffness?: number;

    damping?: number;
  }

/// Weld joint definition. You need to specify local anchor points
/// where they are attached and the relative body angle. The position
/// of the anchor points is important for computing the reaction torque.
  export class WeldJointDef extends JointDef implements IWeldJointDef {
    public readonly localAnchorA: Vec2 = new Vec2();

    public readonly localAnchorB: Vec2 = new Vec2();

    public referenceAngle: number = 0;

    public stiffness: number = 0;

    public damping: number = 0;

    constructor() {
      super(JointType.WeldJoint);
    }

    public initialize(bA: Body, bB: Body, anchor: Vec2): void {
      this.bodyA = bA;
      this.bodyB = bB;
      this.bodyA.getLocalPoint(anchor, this.localAnchorA);
      this.bodyB.getLocalPoint(anchor, this.localAnchorB);
      this.referenceAngle = this.bodyB.getAngle() - this.bodyA.getAngle();
    }
  }

  export class WeldJoint extends Joint {
    public stiffness: number = 0;
    public damping: number = 0;
    public bias: number = 0;

    // Solver shared
    public readonly localAnchorA: Vec2 = new Vec2();
    public readonly localAnchorB: Vec2 = new Vec2();
    public referenceAngle: number = 0;
    public gamma: number = 0;
    public readonly impulse: Vec3 = new Vec3(0, 0, 0);

    // Solver temp
    public indexA: number = 0;
    public indexB: number = 0;
    public readonly rA: Vec2 = new Vec2();
    public readonly rB: Vec2 = new Vec2();
    public readonly localCenterA: Vec2 = new Vec2();
    public readonly localCenterB: Vec2 = new Vec2();
    public invMassA: number = 0;
    public invMassB: number = 0;
    public invIA: number = 0;
    public invIB: number = 0;
    public readonly mass: Mat33 = new Mat33();

    public readonly qA: Rot = new Rot();
    public readonly qB: Rot = new Rot();
    public readonly lalcA: Vec2 = new Vec2();
    public readonly lalcB: Vec2 = new Vec2();
    public readonly K: Mat33 = new Mat33();

    constructor(def: IWeldJointDef) {
      super(def);

      this.stiffness = maybe(def.stiffness, 0);
      this.damping = maybe(def.damping, 0);

      this.localAnchorA.copy(maybe(def.localAnchorA, Vec2.ZERO));
      this.localAnchorB.copy(maybe(def.localAnchorB, Vec2.ZERO));
      this.referenceAngle = maybe(def.referenceAngle, 0);
      this.impulse.setZero();
    }

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

      const aA: number = data.positions[this.indexA].a;
      const vA: Vec2 = data.velocities[this.indexA].v;
      let wA: number = data.velocities[this.indexA].w;

      const aB: number = data.positions[this.indexB].a;
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;

      const qA: Rot = this.qA.setAngle(aA), qB: Rot = this.qB.setAngle(aB);

      // rA = Mul(qA, localAnchorA - localCenterA);
      Vec2.SubVV(this.localAnchorA, this.localCenterA, this.lalcA);
      Rot.mulRV(qA, this.lalcA, this.rA);
      // rB = Mul(qB, localAnchorB - localCenterB);
      Vec2.SubVV(this.localAnchorB, this.localCenterB, this.lalcB);
      Rot.mulRV(qB, this.lalcB, this.rB);

      // J = [-I -r1_skew I r2_skew]
      //     [ 0       -1 0       1]
      // r_skew = [-ry; rx]

      // Matlab
      // K = [ mA+r1y^2*iA+mB+r2y^2*iB,  -r1y*iA*r1x-r2y*iB*r2x,          -r1y*iA-r2y*iB]
      //     [  -r1y*iA*r1x-r2y*iB*r2x, mA+r1x^2*iA+mB+r2x^2*iB,           r1x*iA+r2x*iB]
      //     [          -r1y*iA-r2y*iB,           r1x*iA+r2x*iB,                   iA+iB]

      const mA: number = this.invMassA, mB: number = this.invMassB;
      const iA: number = this.invIA, iB: number = this.invIB;

      const K: Mat33 = this.K;
      K.ex.x = mA + mB + this.rA.y * this.rA.y * iA + this.rB.y * this.rB.y * iB;
      K.ey.x = -this.rA.y * this.rA.x * iA - this.rB.y * this.rB.x * iB;
      K.ez.x = -this.rA.y * iA - this.rB.y * iB;
      K.ex.y = K.ey.x;
      K.ey.y = mA + mB + this.rA.x * this.rA.x * iA + this.rB.x * this.rB.x * iB;
      K.ez.y = this.rA.x * iA + this.rB.x * iB;
      K.ex.z = K.ez.x;
      K.ey.z = K.ez.y;
      K.ez.z = iA + iB;

      if (this.stiffness > 0) {
        K.getInverse22(this.mass);

        let invM: number = iA + iB;

        const C: number = aB - aA - this.referenceAngle;

        // Damping coefficient
        const d: number = this.damping;

        // Spring stiffness
        const k: number = this.stiffness;

        // magic formulas
        const h: number = data.step.dt;
        this.gamma = h * (d + h * k);
        this.gamma = this.gamma !== 0 ? 1 / this.gamma : 0;
        this.bias = C * h * k * this.gamma;

        invM += this.gamma;
        this.mass.ez.z = invM !== 0 ? 1 / invM : 0;
      } else {
        K.getSymInverse33(this.mass);
        this.gamma = 0;
        this.bias = 0;
      }

      if (data.step.warmStarting) {
        // Scale impulses to support a variable time step.
        this.impulse.selfMul(data.step.dtRatio);

        // Vec2 P(impulse.x, impulse.y);
        const P: Vec2 = WeldJoint.initVelocityConstraints_s_P.set(this.impulse.x, this.impulse.y);

        // vA -= mA * P;
        vA.selfMulSub(mA, P);
        wA -= iA * (Vec2.CrossVV(this.rA, P) + this.impulse.z);

        // vB += mB * P;
        vB.selfMulAdd(mB, P);
        wB += iB * (Vec2.CrossVV(this.rB, P) + this.impulse.z);
      } else {
        this.impulse.setZero();
      }

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    private static solveVelocityConstraints_s_Cdot1 = new Vec2();
    private static solveVelocityConstraints_s_impulse1 = new Vec2();
    private static solveVelocityConstraints_s_impulse = new Vec3();
    private static solveVelocityConstraints_s_P = new Vec2();
    public solveVelocityConstraints(data: SolverData): void {
      const vA: Vec2 = data.velocities[this.indexA].v;
      let wA: number = data.velocities[this.indexA].w;
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;

      const mA: number = this.invMassA, mB: number = this.invMassB;
      const iA: number = this.invIA, iB: number = this.invIB;

      if (this.stiffness > 0) {
        const Cdot2: number = wB - wA;

        const impulse2: number = -this.mass.ez.z * (Cdot2 + this.bias + this.gamma * this.impulse.z);
        this.impulse.z += impulse2;

        wA -= iA * impulse2;
        wB += iB * impulse2;

        // Vec2 Cdot1 = vB + Vec2.CrossSV(wB, this.rB) - vA - Vec2.CrossSV(wA, this.rA);
        const Cdot1: Vec2 = Vec2.SubVV(
          Vec2.AddVCrossSV(vB, wB, this.rB, Vec2.s_t0),
          Vec2.AddVCrossSV(vA, wA, this.rA, Vec2.s_t1),
          WeldJoint.solveVelocityConstraints_s_Cdot1);

        // Vec2 impulse1 = -Mul22(mass, Cdot1);
        const impulse1: Vec2 = Mat33.mulM33XY(this.mass, Cdot1.x, Cdot1.y, WeldJoint.solveVelocityConstraints_s_impulse1).selfNeg();
        this.impulse.x += impulse1.x;
        this.impulse.y += impulse1.y;

        // Vec2 P = impulse1;
        const P: Vec2 = impulse1;

        // vA -= mA * P;
        vA.selfMulSub(mA, P);
        // wA -= iA * Cross(rA, P);
        wA -= iA * Vec2.CrossVV(this.rA, P);

        // vB += mB * P;
        vB.selfMulAdd(mB, P);
        // wB += iB * Cross(rB, P);
        wB += iB * Vec2.CrossVV(this.rB, P);
      } else {
        // Vec2 Cdot1 = vB + Cross(wB, this.rB) - vA - Cross(wA, this.rA);
        const Cdot1: Vec2 = Vec2.SubVV(
          Vec2.AddVCrossSV(vB, wB, this.rB, Vec2.s_t0),
          Vec2.AddVCrossSV(vA, wA, this.rA, Vec2.s_t1),
          WeldJoint.solveVelocityConstraints_s_Cdot1);
        const Cdot2: number = wB - wA;
        // Vec3 const Cdot(Cdot1.x, Cdot1.y, Cdot2);

        // Vec3 impulse = -Mul(mass, Cdot);
        const impulse: Vec3 = Mat33.mulM33XYZ(this.mass, Cdot1.x, Cdot1.y, Cdot2, WeldJoint.solveVelocityConstraints_s_impulse).selfNeg();
        this.impulse.selfAdd(impulse);

        // Vec2 P(impulse.x, impulse.y);
        const P: Vec2 = WeldJoint.solveVelocityConstraints_s_P.set(impulse.x, impulse.y);

        // vA -= mA * P;
        vA.selfMulSub(mA, P);
        wA -= iA * (Vec2.CrossVV(this.rA, P) + impulse.z);

        // vB += mB * P;
        vB.selfMulAdd(mB, P);
        wB += iB * (Vec2.CrossVV(this.rB, P) + impulse.z);
      }

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    private static solvePositionConstraints_s_C1 = new Vec2();
    private static solvePositionConstraints_s_P = new Vec2();
    private static solvePositionConstraints_s_impulse = new Vec3();
    public solvePositionConstraints(data: SolverData): boolean {
      const cA: Vec2 = data.positions[this.indexA].c;
      let aA: number = data.positions[this.indexA].a;
      const cB: Vec2 = data.positions[this.indexB].c;
      let aB: number = data.positions[this.indexB].a;

      const qA: Rot = this.qA.setAngle(aA), qB: Rot = this.qB.setAngle(aB);

      const mA: number = this.invMassA, mB: number = this.invMassB;
      const iA: number = this.invIA, iB: number = this.invIB;

      // Vec2 rA = Mul(qA, localAnchorA - localCenterA);
      Vec2.SubVV(this.localAnchorA, this.localCenterA, this.lalcA);
      const rA: Vec2 = Rot.mulRV(qA, this.lalcA, this.rA);
      // Vec2 rB = Mul(qB, localAnchorB - localCenterB);
      Vec2.SubVV(this.localAnchorB, this.localCenterB, this.lalcB);
      const rB: Vec2 = Rot.mulRV(qB, this.lalcB, this.rB);

      let positionError: number, angularError: number;

      const K: Mat33 = this.K;
      K.ex.x = mA + mB + rA.y * rA.y * iA + rB.y * rB.y * iB;
      K.ey.x = -rA.y * rA.x * iA - rB.y * rB.x * iB;
      K.ez.x = -rA.y * iA - rB.y * iB;
      K.ex.y = K.ey.x;
      K.ey.y = mA + mB + rA.x * rA.x * iA + rB.x * rB.x * iB;
      K.ez.y = rA.x * iA + rB.x * iB;
      K.ex.z = K.ez.x;
      K.ey.z = K.ez.y;
      K.ez.z = iA + iB;

      if (this.stiffness > 0) {
        // Vec2 C1 =  cB + rB - cA - rA;
        const C1 =
          Vec2.SubVV(
            Vec2.AddVV(cB, rB, Vec2.s_t0),
            Vec2.AddVV(cA, rA, Vec2.s_t1),
            WeldJoint.solvePositionConstraints_s_C1);
        positionError = C1.length();
        angularError = 0;

        // Vec2 P = -K.Solve22(C1);
        const P: Vec2 = K.solve22(C1.x, C1.y, WeldJoint.solvePositionConstraints_s_P).selfNeg();

        // cA -= mA * P;
        cA.selfMulSub(mA, P);
        aA -= iA * Vec2.CrossVV(rA, P);

        // cB += mB * P;
        cB.selfMulAdd(mB, P);
        aB += iB * Vec2.CrossVV(rB, P);
      } else {
        // Vec2 C1 =  cB + rB - cA - rA;
        const C1 =
          Vec2.SubVV(
            Vec2.AddVV(cB, rB, Vec2.s_t0),
            Vec2.AddVV(cA, rA, Vec2.s_t1),
            WeldJoint.solvePositionConstraints_s_C1);
        const C2: number = aB - aA - this.referenceAngle;

        positionError = C1.length();
        angularError = Abs(C2);

        // Vec3 C(C1.x, C1.y, C2);

        // Vec3 impulse = -K.Solve33(C);
        const impulse: Vec3 = K.solve33(C1.x, C1.y, C2, WeldJoint.solvePositionConstraints_s_impulse).selfNeg();

        // Vec2 P(impulse.x, impulse.y);
        const P: Vec2 = WeldJoint.solvePositionConstraints_s_P.set(impulse.x, impulse.y);

        // cA -= mA * P;
        cA.selfMulSub(mA, P);
        aA -= iA * (Vec2.CrossVV(this.rA, P) + impulse.z);

        // cB += mB * P;
        cB.selfMulAdd(mB, P);
        aB += iB * (Vec2.CrossVV(this.rB, P) + impulse.z);
      }

      // data.positions[this.indexA].c = cA;
      data.positions[this.indexA].a = aA;
      // data.positions[this.indexB].c = cB;
      data.positions[this.indexB].a = aB;

      return positionError <= linearSlop && angularError <= angularSlop;
    }

    public getAnchorA<T extends XY>(out: T): T {
      return this.bodyA.getWorldPoint(this.localAnchorA, out);
    }

    public getAnchorB<T extends XY>(out: T): T {
      return this.bodyB.getWorldPoint(this.localAnchorB, out);
    }

    public getReactionForce<T extends XY>(inv_dt: number, out: T): T {
      // Vec2 P(this.impulse.x, this.impulse.y);
      // return inv_dt * P;
      out.x = inv_dt * this.impulse.x;
      out.y = inv_dt * this.impulse.y;
      return out;
    }

    public getReactionTorque(inv_dt: number): number {
      return inv_dt * this.impulse.z;
    }

    public setStiffness(stiffness: number): void { this.stiffness = stiffness; }
    public getStiffness(): number { return this.stiffness; }

    public setDamping(damping: number) { this.damping = damping; }
    public getDamping() { return this.damping; }

    public dump(log: (format: string, ...args: any[]) => void) {
      const indexA = this.bodyA.islandIndex;
      const indexB = this.bodyB.islandIndex;

      log("  const jd: WeldJointDef = new WeldJointDef();\n");
      log("  jd.bodyA = bodies[%d];\n", indexA);
      log("  jd.bodyB = bodies[%d];\n", indexB);
      log("  jd.collideConnected = %s;\n", (this.collideConnected) ? ("true") : ("false"));
      log("  jd.localAnchorA.Set(%.15f, %.15f);\n", this.localAnchorA.x, this.localAnchorA.y);
      log("  jd.localAnchorB.Set(%.15f, %.15f);\n", this.localAnchorB.x, this.localAnchorB.y);
      log("  jd.referenceAngle = %.15f;\n", this.referenceAngle);
      log("  jd.stiffness = %.15f;\n", this.stiffness);
      log("  jd.damping = %.15f;\n", this.damping);
      log("  joints[%d] = this.world.CreateJoint(jd);\n", this.index);
    }
  }

}
