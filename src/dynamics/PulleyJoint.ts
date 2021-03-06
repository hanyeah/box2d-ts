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
  export const minPulleyLength: number = 2;

  export interface IPulleyJointDef extends IJointDef {
    groundAnchorA?: XY;

    groundAnchorB?: XY;

    localAnchorA?: XY;

    localAnchorB?: XY;

    lengthA?: number;

    lengthB?: number;

    ratio?: number;
  }

/// Pulley joint definition. This requires two ground anchors,
/// two dynamic body anchor points, and a pulley ratio.
  export class PulleyJointDef extends JointDef implements IPulleyJointDef {
    public readonly groundAnchorA: Vec2 = new Vec2(-1, 1);

    public readonly groundAnchorB: Vec2 = new Vec2(1, 1);

    public readonly localAnchorA: Vec2 = new Vec2(-1, 0);

    public readonly localAnchorB: Vec2 = new Vec2(1, 0);

    public lengthA: number = 0;

    public lengthB: number = 0;

    public ratio: number = 1;

    constructor() {
      super(JointType.PulleyJoint);
      this.collideConnected = true;
    }

    public initialize(bA: Body, bB: Body, groundA: XY, groundB: XY, anchorA: XY, anchorB: XY, r: number): void {
      this.bodyA = bA;
      this.bodyB = bB;
      this.groundAnchorA.copy(groundA);
      this.groundAnchorB.copy(groundB);
      this.localAnchorA.copy(anchorA);
      this.localAnchorB.copy(anchorB);
      this.lengthA = Vec2.DistanceVV(bA.getWorldPoint(anchorA, {x: 0, y: 0}), groundA);
      this.lengthB = Vec2.DistanceVV(bB.getWorldPoint(anchorB, {x: 0, y: 0}), groundB);
      this.ratio = r;
      // DEBUG: Assert(this.ratio > epsilon);
    }
  }

  export class PulleyJoint extends Joint {
    public readonly groundAnchorA: Vec2 = new Vec2();
    public readonly groundAnchorB: Vec2 = new Vec2();

    public lengthA: number = 0;
    public lengthB: number = 0;

    // Solver shared
    public readonly localAnchorA: Vec2 = new Vec2();
    public readonly localAnchorB: Vec2 = new Vec2();

    public constant: number = 0;
    public ratio: number = 0;
    public impulse: number = 0;

    // Solver temp
    public indexA: number = 0;
    public indexB: number = 0;
    public readonly uA: Vec2 = new Vec2();
    public readonly uB: Vec2 = new Vec2();
    public readonly rA: Vec2 = new Vec2();
    public readonly rB: Vec2 = new Vec2();
    public readonly localCenterA: Vec2 = new Vec2();
    public readonly localCenterB: Vec2 = new Vec2();

    public invMassA: number = 0;
    public invMassB: number = 0;
    public invIA: number = 0;
    public invIB: number = 0;
    public mass: number = 0;

    public readonly qA: Rot = new Rot();
    public readonly qB: Rot = new Rot();
    public readonly lalcA: Vec2 = new Vec2();
    public readonly lalcB: Vec2 = new Vec2();

    constructor(def: IPulleyJointDef) {
      super(def);

      this.groundAnchorA.copy(maybe(def.groundAnchorA, new Vec2(-1, 1)));
      this.groundAnchorB.copy(maybe(def.groundAnchorB, new Vec2(1, 0)));
      this.localAnchorA.copy(maybe(def.localAnchorA, new Vec2(-1, 0)));
      this.localAnchorB.copy(maybe(def.localAnchorB, new Vec2(1, 0)));

      this.lengthA = maybe(def.lengthA, 0);
      this.lengthB = maybe(def.lengthB, 0);

      // DEBUG: Assert(Maybe(def.ratio, 1) !== 0);
      this.ratio = maybe(def.ratio, 1);

      this.constant = maybe(def.lengthA, 0) + this.ratio * maybe(def.lengthB, 0);

      this.impulse = 0;
    }

    private static initVelocityConstraints_s_PA = new Vec2();
    private static initVelocityConstraints_s_PB = new Vec2();
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

      // Rot qA(aA), qB(aB);
      const qA: Rot = this.qA.setAngle(aA), qB: Rot = this.qB.setAngle(aB);

      // rA = Mul(qA, localAnchorA - localCenterA);
      Vec2.SubVV(this.localAnchorA, this.localCenterA, this.lalcA);
      Rot.mulRV(qA, this.lalcA, this.rA);
      // rB = Mul(qB, localAnchorB - localCenterB);
      Vec2.SubVV(this.localAnchorB, this.localCenterB, this.lalcB);
      Rot.mulRV(qB, this.lalcB, this.rB);

      // Get the pulley axes.
      // uA = cA + rA - groundAnchorA;
      this.uA.copy(cA).selfAdd(this.rA).selfSub(this.groundAnchorA);
      // uB = cB + rB - groundAnchorB;
      this.uB.copy(cB).selfAdd(this.rB).selfSub(this.groundAnchorB);

      const lengthA: number = this.uA.length();
      const lengthB: number = this.uB.length();

      if (lengthA > 10 * linearSlop) {
        this.uA.selfMul(1 / lengthA);
      } else {
        this.uA.setZero();
      }

      if (lengthB > 10 * linearSlop) {
        this.uB.selfMul(1 / lengthB);
      } else {
        this.uB.setZero();
      }

      // Compute effective mass.
      const ruA: number = Vec2.CrossVV(this.rA, this.uA);
      const ruB: number = Vec2.CrossVV(this.rB, this.uB);

      const mA: number = this.invMassA + this.invIA * ruA * ruA;
      const mB: number = this.invMassB + this.invIB * ruB * ruB;

      this.mass = mA + this.ratio * this.ratio * mB;

      if (this.mass > 0) {
        this.mass = 1 / this.mass;
      }

      if (data.step.warmStarting) {
        // Scale impulses to support variable time steps.
        this.impulse *= data.step.dtRatio;

        // Warm starting.
        // Vec2 PA = -(impulse) * uA;
        const PA: Vec2 = Vec2.MulSV(-(this.impulse), this.uA, PulleyJoint.initVelocityConstraints_s_PA);
        // Vec2 PB = (-ratio * impulse) * uB;
        const PB: Vec2 = Vec2.MulSV((-this.ratio * this.impulse), this.uB, PulleyJoint.initVelocityConstraints_s_PB);

        // vA += invMassA * PA;
        vA.selfMulAdd(this.invMassA, PA);
        wA += this.invIA * Vec2.CrossVV(this.rA, PA);
        // vB += invMassB * PB;
        vB.selfMulAdd(this.invMassB, PB);
        wB += this.invIB * Vec2.CrossVV(this.rB, PB);
      } else {
        this.impulse = 0;
      }

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    private static solveVelocityConstraints_s_vpA = new Vec2();
    private static solveVelocityConstraints_s_vpB = new Vec2();
    private static solveVelocityConstraints_s_PA = new Vec2();
    private static solveVelocityConstraints_s_PB = new Vec2();
    public solveVelocityConstraints(data: SolverData): void {
      const vA: Vec2 = data.velocities[this.indexA].v;
      let wA: number = data.velocities[this.indexA].w;
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;

      // Vec2 vpA = vA + Cross(wA, rA);
      const vpA: Vec2 = Vec2.AddVCrossSV(vA, wA, this.rA, PulleyJoint.solveVelocityConstraints_s_vpA);
      // Vec2 vpB = vB + Cross(wB, rB);
      const vpB: Vec2 = Vec2.AddVCrossSV(vB, wB, this.rB, PulleyJoint.solveVelocityConstraints_s_vpB);

      const Cdot: number = -Vec2.DotVV(this.uA, vpA) - this.ratio * Vec2.DotVV(this.uB, vpB);
      const impulse: number = -this.mass * Cdot;
      this.impulse += impulse;

      // Vec2 PA = -impulse * uA;
      const PA: Vec2 = Vec2.MulSV(-impulse, this.uA, PulleyJoint.solveVelocityConstraints_s_PA);
      // Vec2 PB = -ratio * impulse * uB;
      const PB: Vec2 = Vec2.MulSV(-this.ratio * impulse, this.uB, PulleyJoint.solveVelocityConstraints_s_PB);
      // vA += invMassA * PA;
      vA.selfMulAdd(this.invMassA, PA);
      wA += this.invIA * Vec2.CrossVV(this.rA, PA);
      // vB += invMassB * PB;
      vB.selfMulAdd(this.invMassB, PB);
      wB += this.invIB * Vec2.CrossVV(this.rB, PB);

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    private static solvePositionConstraints_s_PA = new Vec2();
    private static solvePositionConstraints_s_PB = new Vec2();
    public solvePositionConstraints(data: SolverData): boolean {
      const cA: Vec2 = data.positions[this.indexA].c;
      let aA: number = data.positions[this.indexA].a;
      const cB: Vec2 = data.positions[this.indexB].c;
      let aB: number = data.positions[this.indexB].a;

      // Rot qA(aA), qB(aB);
      const qA: Rot = this.qA.setAngle(aA), qB: Rot = this.qB.setAngle(aB);

      // Vec2 rA = Mul(qA, localAnchorA - localCenterA);
      Vec2.SubVV(this.localAnchorA, this.localCenterA, this.lalcA);
      const rA: Vec2 = Rot.mulRV(qA, this.lalcA, this.rA);
      // Vec2 rB = Mul(qB, localAnchorB - localCenterB);
      Vec2.SubVV(this.localAnchorB, this.localCenterB, this.lalcB);
      const rB: Vec2 = Rot.mulRV(qB, this.lalcB, this.rB);

      // Get the pulley axes.
      // Vec2 uA = cA + rA - groundAnchorA;
      const uA = this.uA.copy(cA).selfAdd(rA).selfSub(this.groundAnchorA);
      // Vec2 uB = cB + rB - groundAnchorB;
      const uB = this.uB.copy(cB).selfAdd(rB).selfSub(this.groundAnchorB);

      const lengthA: number = uA.length();
      const lengthB: number = uB.length();

      if (lengthA > 10 * linearSlop) {
        uA.selfMul(1 / lengthA);
      } else {
        uA.setZero();
      }

      if (lengthB > 10 * linearSlop) {
        uB.selfMul(1 / lengthB);
      } else {
        uB.setZero();
      }

      // Compute effective mass.
      const ruA: number = Vec2.CrossVV(rA, uA);
      const ruB: number = Vec2.CrossVV(rB, uB);

      const mA: number = this.invMassA + this.invIA * ruA * ruA;
      const mB: number = this.invMassB + this.invIB * ruB * ruB;

      let mass: number = mA + this.ratio * this.ratio * mB;

      if (mass > 0) {
        mass = 1 / mass;
      }

      const C: number = this.constant - lengthA - this.ratio * lengthB;
      const linearError: number = Abs(C);

      const impulse: number = -mass * C;

      // Vec2 PA = -impulse * uA;
      const PA: Vec2 = Vec2.MulSV(-impulse, uA, PulleyJoint.solvePositionConstraints_s_PA);
      // Vec2 PB = -ratio * impulse * uB;
      const PB: Vec2 = Vec2.MulSV(-this.ratio * impulse, uB, PulleyJoint.solvePositionConstraints_s_PB);

      // cA += invMassA * PA;
      cA.selfMulAdd(this.invMassA, PA);
      aA += this.invIA * Vec2.CrossVV(rA, PA);
      // cB += invMassB * PB;
      cB.selfMulAdd(this.invMassB, PB);
      aB += this.invIB * Vec2.CrossVV(rB, PB);

      // data.positions[this.indexA].c = cA;
      data.positions[this.indexA].a = aA;
      // data.positions[this.indexB].c = cB;
      data.positions[this.indexB].a = aB;

      return linearError < linearSlop;
    }

    public getAnchorA<T extends XY>(out: T): T {
      return this.bodyA.getWorldPoint(this.localAnchorA, out);
    }

    public getAnchorB<T extends XY>(out: T): T {
      return this.bodyB.getWorldPoint(this.localAnchorB, out);
    }

    public getReactionForce<T extends XY>(inv_dt: number, out: T): T {
      // Vec2 P = impulse * uB;
      // return inv_dt * P;
      out.x = inv_dt * this.impulse * this.uB.x;
      out.y = inv_dt * this.impulse * this.uB.y;
      return out;
    }

    public getReactionTorque(inv_dt: number): number {
      return 0;
    }

    public getGroundAnchorA() {
      return this.groundAnchorA;
    }

    public getGroundAnchorB() {
      return this.groundAnchorB;
    }

    private static getCurrentLengthA_s_p = new Vec2();
    public getCurrentLengthA() {
      // Vec2 p = bodyA->GetWorldPoint(localAnchorA);
      // Vec2 s = groundAnchorA;
      // Vec2 d = p - s;
      // return d.Length();
      const p = this.bodyA.getWorldPoint(this.localAnchorA, PulleyJoint.getCurrentLengthA_s_p);
      const s = this.groundAnchorA;
      return Vec2.DistanceVV(p, s);
    }

    private static getCurrentLengthB_s_p = new Vec2();
    public getCurrentLengthB() {
      // Vec2 p = bodyB->GetWorldPoint(localAnchorB);
      // Vec2 s = groundAnchorB;
      // Vec2 d = p - s;
      // return d.Length();
      const p = this.bodyB.getWorldPoint(this.localAnchorB, PulleyJoint.getCurrentLengthB_s_p);
      const s = this.groundAnchorB;
      return Vec2.DistanceVV(p, s);
    }

    public dump(log: (format: string, ...args: any[]) => void) {
      const indexA = this.bodyA.islandIndex;
      const indexB = this.bodyB.islandIndex;

      log("  const jd: PulleyJointDef = new PulleyJointDef();\n");
      log("  jd.bodyA = bodies[%d];\n", indexA);
      log("  jd.bodyB = bodies[%d];\n", indexB);
      log("  jd.collideConnected = %s;\n", (this.collideConnected) ? ("true") : ("false"));
      log("  jd.groundAnchorA.Set(%.15f, %.15f);\n", this.groundAnchorA.x, this.groundAnchorA.y);
      log("  jd.groundAnchorB.Set(%.15f, %.15f);\n", this.groundAnchorB.x, this.groundAnchorB.y);
      log("  jd.localAnchorA.Set(%.15f, %.15f);\n", this.localAnchorA.x, this.localAnchorA.y);
      log("  jd.localAnchorB.Set(%.15f, %.15f);\n", this.localAnchorB.x, this.localAnchorB.y);
      log("  jd.lengthA = %.15f;\n", this.lengthA);
      log("  jd.lengthB = %.15f;\n", this.lengthB);
      log("  jd.ratio = %.15f;\n", this.ratio);
      log("  joints[%d] = this.world.CreateJoint(jd);\n", this.index);
    }

    public shiftOrigin(newOrigin: Vec2) {
      this.groundAnchorA.selfSub(newOrigin);
      this.groundAnchorB.selfSub(newOrigin);
    }
  }

}
