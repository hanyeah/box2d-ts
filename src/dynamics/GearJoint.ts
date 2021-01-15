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
  export interface IGearJointDef extends IJointDef {
    joint1: RevoluteJoint | PrismaticJoint;

    joint2: RevoluteJoint | PrismaticJoint;

    ratio?: number;
  }

/// Gear joint definition. This definition requires two existing
/// revolute or prismatic joints (any combination will work).
  export class GearJointDef extends JointDef implements IGearJointDef {
    public joint1!: RevoluteJoint | PrismaticJoint;

    public joint2!: RevoluteJoint | PrismaticJoint;

    public ratio: number = 1;

    constructor() {
      super(JointType.e_gearJoint);
    }
  }

  export class GearJoint extends Joint {
    public joint1: RevoluteJoint | PrismaticJoint;
    public joint2: RevoluteJoint | PrismaticJoint;

    public typeA: JointType = JointType.e_unknownJoint;
    public typeB: JointType = JointType.e_unknownJoint;

    // Body A is connected to body C
    // Body B is connected to body D
    public bodyC: Body;
    public bodyD: Body;

    // Solver shared
    public readonly localAnchorA: Vec2 = new Vec2();
    public readonly localAnchorB: Vec2 = new Vec2();
    public readonly localAnchorC: Vec2 = new Vec2();
    public readonly localAnchorD: Vec2 = new Vec2();

    public readonly localAxisC: Vec2 = new Vec2();
    public readonly localAxisD: Vec2 = new Vec2();

    public referenceAngleA: number = 0;
    public referenceAngleB: number = 0;

    public constant: number = 0;
    public ratio: number = 0;

    public impulse: number = 0;

    // Solver temp
    public indexA: number = 0;
    public indexB: number = 0;
    public indexC: number = 0;
    public indexD: number = 0;
    public readonly lcA: Vec2 = new Vec2();
    public readonly lcB: Vec2 = new Vec2();
    public readonly lcC: Vec2 = new Vec2();
    public readonly lcD: Vec2 = new Vec2();
    public mA: number = 0;
    public mB: number = 0;
    public mC: number = 0;
    public mD: number = 0;
    public iA: number = 0;
    public iB: number = 0;
    public iC: number = 0;
    public iD: number = 0;
    public readonly JvAC: Vec2 = new Vec2();
    public readonly JvBD: Vec2 = new Vec2();
    public JwA: number = 0;
    public JwB: number = 0;
    public JwC: number = 0;
    public JwD: number = 0;
    public mass: number = 0;

    public readonly qA: Rot = new Rot();
    public readonly qB: Rot = new Rot();
    public readonly qC: Rot = new Rot();
    public readonly qD: Rot = new Rot();
    public readonly lalcA: Vec2 = new Vec2();
    public readonly lalcB: Vec2 = new Vec2();
    public readonly lalcC: Vec2 = new Vec2();
    public readonly lalcD: Vec2 = new Vec2();

    constructor(def: IGearJointDef) {
      super(def);

      this.joint1 = def.joint1;
      this.joint2 = def.joint2;

      this.typeA = this.joint1.GetType();
      this.typeB = this.joint2.GetType();

      // DEBUG: Assert(this.typeA === JointType.e_revoluteJoint || this.typeA === JointType.e_prismaticJoint);
      // DEBUG: Assert(this.typeB === JointType.e_revoluteJoint || this.typeB === JointType.e_prismaticJoint);

      let coordinateA: number, coordinateB: number;

      // TODO_ERIN there might be some problem with the joint edges in Joint.

      this.bodyC = this.joint1.GetBodyA();
      this.bodyA = this.joint1.GetBodyB();

      // Body B on joint1 must be dynamic
      // Assert(this.bodyA.type === dynamicBody);

      // Get geometry of joint1
      const xfA: Transform = this.bodyA.xf;
      const aA: number = this.bodyA.sweep.a;
      const xfC: Transform = this.bodyC.xf;
      const aC: number = this.bodyC.sweep.a;

      if (this.typeA === JointType.e_revoluteJoint) {
        const revolute: RevoluteJoint = def.joint1 as RevoluteJoint;
        this.localAnchorC.Copy(revolute.localAnchorA);
        this.localAnchorA.Copy(revolute.localAnchorB);
        this.referenceAngleA = revolute.referenceAngle;
        this.localAxisC.SetZero();

        coordinateA = aA - aC - this.referenceAngleA;
      } else {
        const prismatic: PrismaticJoint = def.joint1 as PrismaticJoint;
        this.localAnchorC.Copy(prismatic.localAnchorA);
        this.localAnchorA.Copy(prismatic.localAnchorB);
        this.referenceAngleA = prismatic.referenceAngle;
        this.localAxisC.Copy(prismatic.localXAxisA);

        // Vec2 pC = localAnchorC;
        const pC = this.localAnchorC;
        // Vec2 pA = MulT(xfC.q, Mul(xfA.q, localAnchorA) + (xfA.p - xfC.p));
        const pA: Vec2 = Rot.MulTRV(
          xfC.q,
          Vec2.AddVV(
            Rot.MulRV(xfA.q, this.localAnchorA, Vec2.s_t0),
            Vec2.SubVV(xfA.p, xfC.p, Vec2.s_t1),
            Vec2.s_t0),
          Vec2.s_t0); // pA uses s_t0
        // coordinateA = Dot(pA - pC, localAxisC);
        coordinateA = Vec2.DotVV(Vec2.SubVV(pA, pC, Vec2.s_t0), this.localAxisC);
      }

      this.bodyD = this.joint2.GetBodyA();
      this.bodyB = this.joint2.GetBodyB();

      // Body B on joint2 must be dynamic
      // Assert(this.bodyB.type === dynamicBody);

      // Get geometry of joint2
      const xfB: Transform = this.bodyB.xf;
      const aB: number = this.bodyB.sweep.a;
      const xfD: Transform = this.bodyD.xf;
      const aD: number = this.bodyD.sweep.a;

      if (this.typeB === JointType.e_revoluteJoint) {
        const revolute: RevoluteJoint = def.joint2 as RevoluteJoint;
        this.localAnchorD.Copy(revolute.localAnchorA);
        this.localAnchorB.Copy(revolute.localAnchorB);
        this.referenceAngleB = revolute.referenceAngle;
        this.localAxisD.SetZero();

        coordinateB = aB - aD - this.referenceAngleB;
      } else {
        const prismatic: PrismaticJoint = def.joint2 as PrismaticJoint;
        this.localAnchorD.Copy(prismatic.localAnchorA);
        this.localAnchorB.Copy(prismatic.localAnchorB);
        this.referenceAngleB = prismatic.referenceAngle;
        this.localAxisD.Copy(prismatic.localXAxisA);

        // Vec2 pD = localAnchorD;
        const pD = this.localAnchorD;
        // Vec2 pB = MulT(xfD.q, Mul(xfB.q, localAnchorB) + (xfB.p - xfD.p));
        const pB: Vec2 = Rot.MulTRV(
          xfD.q,
          Vec2.AddVV(
            Rot.MulRV(xfB.q, this.localAnchorB, Vec2.s_t0),
            Vec2.SubVV(xfB.p, xfD.p, Vec2.s_t1),
            Vec2.s_t0),
          Vec2.s_t0); // pB uses s_t0
        // coordinateB = Dot(pB - pD, localAxisD);
        coordinateB = Vec2.DotVV(Vec2.SubVV(pB, pD, Vec2.s_t0), this.localAxisD);
      }

      this.ratio = Maybe(def.ratio, 1);

      this.constant = coordinateA + this.ratio * coordinateB;

      this.impulse = 0;
    }

    private static InitVelocityConstraints_s_u = new Vec2();
    private static InitVelocityConstraints_s_rA = new Vec2();
    private static InitVelocityConstraints_s_rB = new Vec2();
    private static InitVelocityConstraints_s_rC = new Vec2();
    private static InitVelocityConstraints_s_rD = new Vec2();
    public InitVelocityConstraints(data: SolverData): void {
      this.indexA = this.bodyA.islandIndex;
      this.indexB = this.bodyB.islandIndex;
      this.indexC = this.bodyC.islandIndex;
      this.indexD = this.bodyD.islandIndex;
      this.lcA.Copy(this.bodyA.sweep.localCenter);
      this.lcB.Copy(this.bodyB.sweep.localCenter);
      this.lcC.Copy(this.bodyC.sweep.localCenter);
      this.lcD.Copy(this.bodyD.sweep.localCenter);
      this.mA = this.bodyA.invMass;
      this.mB = this.bodyB.invMass;
      this.mC = this.bodyC.invMass;
      this.mD = this.bodyD.invMass;
      this.iA = this.bodyA.invI;
      this.iB = this.bodyB.invI;
      this.iC = this.bodyC.invI;
      this.iD = this.bodyD.invI;

      const aA: number = data.positions[this.indexA].a;
      const vA: Vec2 = data.velocities[this.indexA].v;
      let wA: number = data.velocities[this.indexA].w;

      const aB: number = data.positions[this.indexB].a;
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;

      const aC: number = data.positions[this.indexC].a;
      const vC: Vec2 = data.velocities[this.indexC].v;
      let wC: number = data.velocities[this.indexC].w;

      const aD: number = data.positions[this.indexD].a;
      const vD: Vec2 = data.velocities[this.indexD].v;
      let wD: number = data.velocities[this.indexD].w;

      // Rot qA(aA), qB(aB), qC(aC), qD(aD);
      const qA: Rot = this.qA.SetAngle(aA),
        qB: Rot = this.qB.SetAngle(aB),
        qC: Rot = this.qC.SetAngle(aC),
        qD: Rot = this.qD.SetAngle(aD);

      this.mass = 0;

      if (this.typeA === JointType.e_revoluteJoint) {
        this.JvAC.SetZero();
        this.JwA = 1;
        this.JwC = 1;
        this.mass += this.iA + this.iC;
      } else {
        // Vec2 u = Mul(qC, localAxisC);
        const u: Vec2 = Rot.MulRV(qC, this.localAxisC, GearJoint.InitVelocityConstraints_s_u);
        // Vec2 rC = Mul(qC, localAnchorC - lcC);
        Vec2.SubVV(this.localAnchorC, this.lcC, this.lalcC);
        const rC: Vec2 = Rot.MulRV(qC, this.lalcC, GearJoint.InitVelocityConstraints_s_rC);
        // Vec2 rA = Mul(qA, localAnchorA - lcA);
        Vec2.SubVV(this.localAnchorA, this.lcA, this.lalcA);
        const rA: Vec2 = Rot.MulRV(qA, this.lalcA, GearJoint.InitVelocityConstraints_s_rA);
        // JvAC = u;
        this.JvAC.Copy(u);
        // JwC = Cross(rC, u);
        this.JwC = Vec2.CrossVV(rC, u);
        // JwA = Cross(rA, u);
        this.JwA = Vec2.CrossVV(rA, u);
        this.mass += this.mC + this.mA + this.iC * this.JwC * this.JwC + this.iA * this.JwA * this.JwA;
      }

      if (this.typeB === JointType.e_revoluteJoint) {
        this.JvBD.SetZero();
        this.JwB = this.ratio;
        this.JwD = this.ratio;
        this.mass += this.ratio * this.ratio * (this.iB + this.iD);
      } else {
        // Vec2 u = Mul(qD, localAxisD);
        const u: Vec2 = Rot.MulRV(qD, this.localAxisD, GearJoint.InitVelocityConstraints_s_u);
        // Vec2 rD = Mul(qD, localAnchorD - lcD);
        Vec2.SubVV(this.localAnchorD, this.lcD, this.lalcD);
        const rD: Vec2 = Rot.MulRV(qD, this.lalcD, GearJoint.InitVelocityConstraints_s_rD);
        // Vec2 rB = Mul(qB, localAnchorB - lcB);
        Vec2.SubVV(this.localAnchorB, this.lcB, this.lalcB);
        const rB: Vec2 = Rot.MulRV(qB, this.lalcB, GearJoint.InitVelocityConstraints_s_rB);
        // JvBD = ratio * u;
        Vec2.MulSV(this.ratio, u, this.JvBD);
        // JwD = ratio * Cross(rD, u);
        this.JwD = this.ratio * Vec2.CrossVV(rD, u);
        // JwB = ratio * Cross(rB, u);
        this.JwB = this.ratio * Vec2.CrossVV(rB, u);
        this.mass += this.ratio * this.ratio * (this.mD + this.mB) + this.iD * this.JwD * this.JwD + this.iB * this.JwB * this.JwB;
      }

      // Compute effective mass.
      this.mass = this.mass > 0 ? 1 / this.mass : 0;

      if (data.step.warmStarting) {
        // vA += (mA * impulse) * JvAC;
        vA.SelfMulAdd(this.mA * this.impulse, this.JvAC);
        wA += this.iA * this.impulse * this.JwA;
        // vB += (mB * impulse) * JvBD;
        vB.SelfMulAdd(this.mB * this.impulse, this.JvBD);
        wB += this.iB * this.impulse * this.JwB;
        // vC -= (mC * impulse) * JvAC;
        vC.SelfMulSub(this.mC * this.impulse, this.JvAC);
        wC -= this.iC * this.impulse * this.JwC;
        // vD -= (mD * impulse) * JvBD;
        vD.SelfMulSub(this.mD * this.impulse, this.JvBD);
        wD -= this.iD * this.impulse * this.JwD;
      } else {
        this.impulse = 0;
      }

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
      // data.velocities[this.indexC].v = vC;
      data.velocities[this.indexC].w = wC;
      // data.velocities[this.indexD].v = vD;
      data.velocities[this.indexD].w = wD;
    }

    public SolveVelocityConstraints(data: SolverData): void {
      const vA: Vec2 = data.velocities[this.indexA].v;
      let wA: number = data.velocities[this.indexA].w;
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;
      const vC: Vec2 = data.velocities[this.indexC].v;
      let wC: number = data.velocities[this.indexC].w;
      const vD: Vec2 = data.velocities[this.indexD].v;
      let wD: number = data.velocities[this.indexD].w;

      // float32 Cdot = Dot(JvAC, vA - vC) + Dot(JvBD, vB - vD);
      let Cdot =
        Vec2.DotVV(this.JvAC, Vec2.SubVV(vA, vC, Vec2.s_t0)) +
        Vec2.DotVV(this.JvBD, Vec2.SubVV(vB, vD, Vec2.s_t0));
      Cdot += (this.JwA * wA - this.JwC * wC) + (this.JwB * wB - this.JwD * wD);

      const impulse: number = -this.mass * Cdot;
      this.impulse += impulse;

      // vA += (mA * impulse) * JvAC;
      vA.SelfMulAdd((this.mA * impulse), this.JvAC);
      wA += this.iA * impulse * this.JwA;
      // vB += (mB * impulse) * JvBD;
      vB.SelfMulAdd((this.mB * impulse), this.JvBD);
      wB += this.iB * impulse * this.JwB;
      // vC -= (mC * impulse) * JvAC;
      vC.SelfMulSub((this.mC * impulse), this.JvAC);
      wC -= this.iC * impulse * this.JwC;
      // vD -= (mD * impulse) * JvBD;
      vD.SelfMulSub((this.mD * impulse), this.JvBD);
      wD -= this.iD * impulse * this.JwD;

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
      // data.velocities[this.indexC].v = vC;
      data.velocities[this.indexC].w = wC;
      // data.velocities[this.indexD].v = vD;
      data.velocities[this.indexD].w = wD;
    }

    private static SolvePositionConstraints_s_u = new Vec2();
    private static SolvePositionConstraints_s_rA = new Vec2();
    private static SolvePositionConstraints_s_rB = new Vec2();
    private static SolvePositionConstraints_s_rC = new Vec2();
    private static SolvePositionConstraints_s_rD = new Vec2();
    public SolvePositionConstraints(data: SolverData): boolean {
      const cA: Vec2 = data.positions[this.indexA].c;
      let aA: number = data.positions[this.indexA].a;
      const cB: Vec2 = data.positions[this.indexB].c;
      let aB: number = data.positions[this.indexB].a;
      const cC: Vec2 = data.positions[this.indexC].c;
      let aC: number = data.positions[this.indexC].a;
      const cD: Vec2 = data.positions[this.indexD].c;
      let aD: number = data.positions[this.indexD].a;

      // Rot qA(aA), qB(aB), qC(aC), qD(aD);
      const qA: Rot = this.qA.SetAngle(aA),
        qB: Rot = this.qB.SetAngle(aB),
        qC: Rot = this.qC.SetAngle(aC),
        qD: Rot = this.qD.SetAngle(aD);

      const linearError: number = 0;

      let coordinateA: number, coordinateB: number;

      const JvAC: Vec2 = this.JvAC, JvBD: Vec2 = this.JvBD;
      let JwA: number, JwB: number, JwC: number, JwD: number;
      let mass: number = 0;

      if (this.typeA === JointType.e_revoluteJoint) {
        JvAC.SetZero();
        JwA = 1;
        JwC = 1;
        mass += this.iA + this.iC;

        coordinateA = aA - aC - this.referenceAngleA;
      } else {
        // Vec2 u = Mul(qC, localAxisC);
        const u: Vec2 = Rot.MulRV(qC, this.localAxisC, GearJoint.SolvePositionConstraints_s_u);
        // Vec2 rC = Mul(qC, localAnchorC - lcC);
        const rC: Vec2 = Rot.MulRV(qC, this.lalcC, GearJoint.SolvePositionConstraints_s_rC);
        // Vec2 rA = Mul(qA, localAnchorA - lcA);
        const rA: Vec2 = Rot.MulRV(qA, this.lalcA, GearJoint.SolvePositionConstraints_s_rA);
        // JvAC = u;
        JvAC.Copy(u);
        // JwC = Cross(rC, u);
        JwC = Vec2.CrossVV(rC, u);
        // JwA = Cross(rA, u);
        JwA = Vec2.CrossVV(rA, u);
        mass += this.mC + this.mA + this.iC * JwC * JwC + this.iA * JwA * JwA;

        // Vec2 pC = localAnchorC - lcC;
        const pC = this.lalcC;
        // Vec2 pA = MulT(qC, rA + (cA - cC));
        const pA: Vec2 = Rot.MulTRV(
          qC,
          Vec2.AddVV(
            rA,
            Vec2.SubVV(cA, cC, Vec2.s_t0),
            Vec2.s_t0),
          Vec2.s_t0); // pA uses s_t0
        // coordinateA = Dot(pA - pC, localAxisC);
        coordinateA = Vec2.DotVV(Vec2.SubVV(pA, pC, Vec2.s_t0), this.localAxisC);
      }

      if (this.typeB === JointType.e_revoluteJoint) {
        JvBD.SetZero();
        JwB = this.ratio;
        JwD = this.ratio;
        mass += this.ratio * this.ratio * (this.iB + this.iD);

        coordinateB = aB - aD - this.referenceAngleB;
      } else {
        // Vec2 u = Mul(qD, localAxisD);
        const u: Vec2 = Rot.MulRV(qD, this.localAxisD, GearJoint.SolvePositionConstraints_s_u);
        // Vec2 rD = Mul(qD, localAnchorD - lcD);
        const rD: Vec2 = Rot.MulRV(qD, this.lalcD, GearJoint.SolvePositionConstraints_s_rD);
        // Vec2 rB = Mul(qB, localAnchorB - lcB);
        const rB: Vec2 = Rot.MulRV(qB, this.lalcB, GearJoint.SolvePositionConstraints_s_rB);
        // JvBD = ratio * u;
        Vec2.MulSV(this.ratio, u, JvBD);
        // JwD = ratio * Cross(rD, u);
        JwD = this.ratio * Vec2.CrossVV(rD, u);
        // JwB = ratio * Cross(rB, u);
        JwB = this.ratio * Vec2.CrossVV(rB, u);
        mass += this.ratio * this.ratio * (this.mD + this.mB) + this.iD * JwD * JwD + this.iB * JwB * JwB;

        // Vec2 pD = localAnchorD - lcD;
        const pD = this.lalcD;
        // Vec2 pB = MulT(qD, rB + (cB - cD));
        const pB: Vec2 = Rot.MulTRV(
          qD,
          Vec2.AddVV(
            rB,
            Vec2.SubVV(cB, cD, Vec2.s_t0),
            Vec2.s_t0),
          Vec2.s_t0); // pB uses s_t0
        // coordinateB = Dot(pB - pD, localAxisD);
        coordinateB = Vec2.DotVV(Vec2.SubVV(pB, pD, Vec2.s_t0), this.localAxisD);
      }

      const C: number = (coordinateA + this.ratio * coordinateB) - this.constant;

      let impulse: number = 0;
      if (mass > 0) {
        impulse = -C / mass;
      }

      // cA += mA * impulse * JvAC;
      cA.SelfMulAdd(this.mA * impulse, JvAC);
      aA += this.iA * impulse * JwA;
      // cB += mB * impulse * JvBD;
      cB.SelfMulAdd(this.mB * impulse, JvBD);
      aB += this.iB * impulse * JwB;
      // cC -= mC * impulse * JvAC;
      cC.SelfMulSub(this.mC * impulse, JvAC);
      aC -= this.iC * impulse * JwC;
      // cD -= mD * impulse * JvBD;
      cD.SelfMulSub(this.mD * impulse, JvBD);
      aD -= this.iD * impulse * JwD;

      // data.positions[this.indexA].c = cA;
      data.positions[this.indexA].a = aA;
      // data.positions[this.indexB].c = cB;
      data.positions[this.indexB].a = aB;
      // data.positions[this.indexC].c = cC;
      data.positions[this.indexC].a = aC;
      // data.positions[this.indexD].c = cD;
      data.positions[this.indexD].a = aD;

      // TODO_ERIN not implemented
      return linearError < linearSlop;
    }

    public GetAnchorA<T extends XY>(out: T): T {
      return this.bodyA.GetWorldPoint(this.localAnchorA, out);
    }

    public GetAnchorB<T extends XY>(out: T): T {
      return this.bodyB.GetWorldPoint(this.localAnchorB, out);
    }

    public GetReactionForce<T extends XY>(inv_dt: number, out: T): T {
      // Vec2 P = impulse * JvAC;
      // return inv_dt * P;
      return Vec2.MulSV(inv_dt * this.impulse, this.JvAC, out);
    }

    public GetReactionTorque(inv_dt: number): number {
      // float32 L = impulse * JwA;
      // return inv_dt * L;
      return inv_dt * this.impulse * this.JwA;
    }

    public GetJoint1() { return this.joint1; }

    public GetJoint2() { return this.joint2; }

    public GetRatio() {
      return this.ratio;
    }

    public SetRatio(ratio: number): void {
      // DEBUG: Assert(IsValid(ratio));
      this.ratio = ratio;
    }

    public Dump(log: (format: string, ...args: any[]) => void) {
      const indexA = this.bodyA.islandIndex;
      const indexB = this.bodyB.islandIndex;

      const index1 = this.joint1.index;
      const index2 = this.joint2.index;

      log("  const jd: GearJointDef = new GearJointDef();\n");
      log("  jd.bodyA = bodies[%d];\n", indexA);
      log("  jd.bodyB = bodies[%d];\n", indexB);
      log("  jd.collideConnected = %s;\n", (this.collideConnected) ? ("true") : ("false"));
      log("  jd.joint1 = joints[%d];\n", index1);
      log("  jd.joint2 = joints[%d];\n", index2);
      log("  jd.ratio = %.15f;\n", this.ratio);
      log("  joints[%d] = this.world.CreateJoint(jd);\n", this.index);
    }
  }

}
