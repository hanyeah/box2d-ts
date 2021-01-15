/*
* Copyright (c) 2006-2007 Erin Catto http://www.box2d.org
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
  export interface IFrictionJointDef extends IJointDef {
    localAnchorA?: XY;

    localAnchorB?: XY;

    maxForce?: number;

    maxTorque?: number;
  }

/// Friction joint definition.
  export class FrictionJointDef extends JointDef implements IFrictionJointDef {
    public readonly localAnchorA: Vec2 = new Vec2();

    public readonly localAnchorB: Vec2 = new Vec2();

    public maxForce: number = 0;

    public maxTorque: number = 0;

    constructor() {
      super(JointType.e_frictionJoint);
    }

    public Initialize(bA: Body, bB: Body, anchor: Vec2): void {
      this.bodyA = bA;
      this.bodyB = bB;
      this.bodyA.GetLocalPoint(anchor, this.localAnchorA);
      this.bodyB.GetLocalPoint(anchor, this.localAnchorB);
    }
  }

  export class FrictionJoint extends Joint {
    public readonly localAnchorA: Vec2 = new Vec2();
    public readonly localAnchorB: Vec2 = new Vec2();

    // Solver shared
    public readonly linearImpulse: Vec2 = new Vec2();
    public angularImpulse: number = 0;
    public maxForce: number = 0;
    public maxTorque: number = 0;

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
    public readonly linearMass: Mat22 = new Mat22();
    public angularMass: number = 0;

    public readonly qA: Rot = new Rot();
    public readonly qB: Rot = new Rot();
    public readonly lalcA: Vec2 = new Vec2();
    public readonly lalcB: Vec2 = new Vec2();
    public readonly K: Mat22 = new Mat22();

    constructor(def: IFrictionJointDef) {
      super(def);

      this.localAnchorA.Copy(Maybe(def.localAnchorA, Vec2.ZERO));
      this.localAnchorB.Copy(Maybe(def.localAnchorB, Vec2.ZERO));

      this.linearImpulse.SetZero();
      this.maxForce = Maybe(def.maxForce, 0);
      this.maxTorque = Maybe(def.maxTorque, 0);

      this.linearMass.SetZero();
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

      // const cA: Vec2 = data.positions[this.indexA].c;
      const aA: number = data.positions[this.indexA].a;
      const vA: Vec2 = data.velocities[this.indexA].v;
      let wA: number = data.velocities[this.indexA].w;

      // const cB: Vec2 = data.positions[this.indexB].c;
      const aB: number = data.positions[this.indexB].a;
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;

      // const qA: Rot = new Rot(aA), qB: Rot = new Rot(aB);
      const qA: Rot = this.qA.SetAngle(aA), qB: Rot = this.qB.SetAngle(aB);

      // Compute the effective mass matrix.
      // rA = Mul(qA, localAnchorA - localCenterA);
      Vec2.SubVV(this.localAnchorA, this.localCenterA, this.lalcA);
      const rA: Vec2 = Rot.MulRV(qA, this.lalcA, this.rA);
      // rB = Mul(qB, localAnchorB - localCenterB);
      Vec2.SubVV(this.localAnchorB, this.localCenterB, this.lalcB);
      const rB: Vec2 = Rot.MulRV(qB, this.lalcB, this.rB);

      // J = [-I -r1_skew I r2_skew]
      //     [ 0       -1 0       1]
      // r_skew = [-ry; rx]

      // Matlab
      // K = [ mA+r1y^2*iA+mB+r2y^2*iB,  -r1y*iA*r1x-r2y*iB*r2x,          -r1y*iA-r2y*iB]
      //     [  -r1y*iA*r1x-r2y*iB*r2x, mA+r1x^2*iA+mB+r2x^2*iB,           r1x*iA+r2x*iB]
      //     [          -r1y*iA-r2y*iB,           r1x*iA+r2x*iB,                   iA+iB]

      const mA: number = this.invMassA, mB: number = this.invMassB;
      const iA: number = this.invIA, iB: number = this.invIB;

      const K: Mat22 = this.K; // new Mat22();
      K.ex.x = mA + mB + iA * rA.y * rA.y + iB * rB.y * rB.y;
      K.ex.y = -iA * rA.x * rA.y - iB * rB.x * rB.y;
      K.ey.x = K.ex.y;
      K.ey.y = mA + mB + iA * rA.x * rA.x + iB * rB.x * rB.x;

      K.GetInverse(this.linearMass);

      this.angularMass = iA + iB;
      if (this.angularMass > 0) {
        this.angularMass = 1 / this.angularMass;
      }

      if (data.step.warmStarting) {
        // Scale impulses to support a variable time step.
        // linearImpulse *= data.step.dtRatio;
        this.linearImpulse.SelfMul(data.step.dtRatio);
        this.angularImpulse *= data.step.dtRatio;

        // const P: Vec2(linearImpulse.x, linearImpulse.y);
        const P: Vec2 = this.linearImpulse;

        // vA -= mA * P;
        vA.SelfMulSub(mA, P);
        // wA -= iA * (Cross(rA, P) + angularImpulse);
        wA -= iA * (Vec2.CrossVV(this.rA, P) + this.angularImpulse);
        // vB += mB * P;
        vB.SelfMulAdd(mB, P);
        // wB += iB * (Cross(rB, P) + angularImpulse);
        wB += iB * (Vec2.CrossVV(this.rB, P) + this.angularImpulse);
      } else {
        this.linearImpulse.SetZero();
        this.angularImpulse = 0;
      }

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    private static SolveVelocityConstraints_s_Cdot_v2 = new Vec2();
    private static SolveVelocityConstraints_s_impulseV = new Vec2();
    private static SolveVelocityConstraints_s_oldImpulseV = new Vec2();
    public SolveVelocityConstraints(data: SolverData): void {
      const vA: Vec2 = data.velocities[this.indexA].v;
      let wA: number = data.velocities[this.indexA].w;
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;

      const mA: number = this.invMassA, mB: number = this.invMassB;
      const iA: number = this.invIA, iB: number = this.invIB;

      const h: number = data.step.dt;

      // Solve angular friction
      {
        const Cdot: number = wB - wA;
        let impulse: number = (-this.angularMass * Cdot);

        const oldImpulse: number = this.angularImpulse;
        const maxImpulse: number = h * this.maxTorque;
        this.angularImpulse = Clamp(this.angularImpulse + impulse, (-maxImpulse), maxImpulse);
        impulse = this.angularImpulse - oldImpulse;

        wA -= iA * impulse;
        wB += iB * impulse;
      }

      // Solve linear friction
      {
        // Vec2 Cdot = vB + Cross(wB, rB) - vA - Cross(wA, rA);
        const Cdot_v2: Vec2 = Vec2.SubVV(
          Vec2.AddVCrossSV(vB, wB, this.rB, Vec2.s_t0),
          Vec2.AddVCrossSV(vA, wA, this.rA, Vec2.s_t1),
          FrictionJoint.SolveVelocityConstraints_s_Cdot_v2);

        // Vec2 impulse = -Mul(linearMass, Cdot);
        const impulseV: Vec2 = Mat22.MulMV(this.linearMass, Cdot_v2, FrictionJoint.SolveVelocityConstraints_s_impulseV).SelfNeg();
        // Vec2 oldImpulse = linearImpulse;
        const oldImpulseV = FrictionJoint.SolveVelocityConstraints_s_oldImpulseV.Copy(this.linearImpulse);
        // linearImpulse += impulse;
        this.linearImpulse.SelfAdd(impulseV);

        const maxImpulse: number = h * this.maxForce;

        if (this.linearImpulse.LengthSquared() > maxImpulse * maxImpulse) {
          this.linearImpulse.Normalize();
          this.linearImpulse.SelfMul(maxImpulse);
        }

        // impulse = linearImpulse - oldImpulse;
        Vec2.SubVV(this.linearImpulse, oldImpulseV, impulseV);

        // vA -= mA * impulse;
        vA.SelfMulSub(mA, impulseV);
        // wA -= iA * Cross(rA, impulse);
        wA -= iA * Vec2.CrossVV(this.rA, impulseV);

        // vB += mB * impulse;
        vB.SelfMulAdd(mB, impulseV);
        // wB += iB * Cross(rB, impulse);
        wB += iB * Vec2.CrossVV(this.rB, impulseV);
      }

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    public SolvePositionConstraints(data: SolverData): boolean {
      return true;
    }

    public GetAnchorA<T extends XY>(out: T): T {
      return this.bodyA.GetWorldPoint(this.localAnchorA, out);
    }

    public GetAnchorB<T extends XY>(out: T): T {
      return this.bodyB.GetWorldPoint(this.localAnchorB, out);
    }

    public GetReactionForce<T extends XY>(inv_dt: number, out: T): T {
      out.x = inv_dt * this.linearImpulse.x;
      out.y = inv_dt * this.linearImpulse.y;
      return out;
    }

    public GetReactionTorque(inv_dt: number): number {
      return inv_dt * this.angularImpulse;
    }

    public GetLocalAnchorA(): Vec2 { return this.localAnchorA; }

    public GetLocalAnchorB(): Vec2 { return this.localAnchorB; }

    public SetMaxForce(force: number): void {
      this.maxForce = force;
    }

    public GetMaxForce(): number {
      return this.maxForce;
    }

    public SetMaxTorque(torque: number): void {
      this.maxTorque = torque;
    }

    public GetMaxTorque(): number {
      return this.maxTorque;
    }

    public Dump(log: (format: string, ...args: any[]) => void): void {
      const indexA: number = this.bodyA.islandIndex;
      const indexB: number = this.bodyB.islandIndex;

      log("  const jd: FrictionJointDef = new FrictionJointDef();\n");
      log("  jd.bodyA = bodies[%d];\n", indexA);
      log("  jd.bodyB = bodies[%d];\n", indexB);
      log("  jd.collideConnected = %s;\n", (this.collideConnected) ? ("true") : ("false"));
      log("  jd.localAnchorA.Set(%.15f, %.15f);\n", this.localAnchorA.x, this.localAnchorA.y);
      log("  jd.localAnchorB.Set(%.15f, %.15f);\n", this.localAnchorB.x, this.localAnchorB.y);
      log("  jd.maxForce = %.15f;\n", this.maxForce);
      log("  jd.maxTorque = %.15f;\n", this.maxTorque);
      log("  joints[%d] = this.world.CreateJoint(jd);\n", this.index);
    }
  }

}
