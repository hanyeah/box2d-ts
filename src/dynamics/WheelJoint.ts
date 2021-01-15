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
  export interface IWheelJointDef extends IJointDef {
    /// The local anchor point relative to bodyA's origin.
    localAnchorA?: XY;

    /// The local anchor point relative to bodyB's origin.
    localAnchorB?: XY;

    /// The local translation axis in bodyA.
    localAxisA?: XY;

    /// Enable/disable the joint limit.
    enableLimit?: boolean;

    /// The lower translation limit, usually in meters.
    lowerTranslation?: number;

    /// The upper translation limit, usually in meters.
    upperTranslation?: number;

    /// Enable/disable the joint motor.
    enableMotor?: boolean;

    /// The maximum motor torque, usually in N-m.
    maxMotorTorque?: number;

    /// The desired motor speed in radians per second.
    motorSpeed?: number;

    /// Suspension stiffness. Typically in units N/m.
    stiffness?: number;

    /// Suspension damping. Typically in units of N*s/m.
    damping?: number;
  }

/// Wheel joint definition. This requires defining a line of
/// motion using an axis and an anchor point. The definition uses local
/// anchor points and a local axis so that the initial configuration
/// can violate the constraint slightly. The joint translation is zero
/// when the local anchor points coincide in world space. Using local
/// anchors and a local axis helps when saving and loading a game.
  export class WheelJointDef extends JointDef implements IWheelJointDef {
    public readonly localAnchorA: Vec2 = new Vec2(0, 0);

    public readonly localAnchorB: Vec2 = new Vec2(0, 0);

    public readonly localAxisA: Vec2 = new Vec2(1, 0);

    public enableLimit: boolean = false;

    public lowerTranslation: number = 0;

    public upperTranslation: number = 0;

    public enableMotor = false;

    public maxMotorTorque: number = 0;

    public motorSpeed: number = 0;

    public stiffness: number = 0;

    public damping: number = 0;

    constructor() {
      super(JointType.WheelJoint);
    }

    public initialize(bA: Body, bB: Body, anchor: Vec2, axis: Vec2): void {
      this.bodyA = bA;
      this.bodyB = bB;
      this.bodyA.getLocalPoint(anchor, this.localAnchorA);
      this.bodyB.getLocalPoint(anchor, this.localAnchorB);
      this.bodyA.getLocalVector(axis, this.localAxisA);
    }
  }

  export class WheelJoint extends Joint {
    public readonly localAnchorA: Vec2 = new Vec2();
    public readonly localAnchorB: Vec2 = new Vec2();
    public readonly localXAxisA: Vec2 = new Vec2();
    public readonly localYAxisA: Vec2 = new Vec2();

    public impulse: number = 0;
    public motorImpulse: number = 0;
    public springImpulse: number = 0;

    public lowerImpulse: number = 0;
    public upperImpulse: number = 0;
    public translation: number = 0;
    public lowerTranslation: number = 0;
    public upperTranslation: number = 0;

    public maxMotorTorque: number = 0;
    public motorSpeed: number = 0;

    public enableLimit = false;
    public enableMotor = false;

    public stiffness: number = 0;
    public damping: number = 0;

    // Solver temp
    public indexA: number = 0;
    public indexB: number = 0;
    public readonly localCenterA: Vec2 = new Vec2();
    public readonly localCenterB: Vec2 = new Vec2();
    public invMassA: number = 0;
    public invMassB: number = 0;
    public invIA: number = 0;
    public invIB: number = 0;

    public readonly ax: Vec2 = new Vec2();
    public readonly ay: Vec2 = new Vec2();
    public sAx: number = 0;
    public sBx: number = 0;
    public sAy: number = 0;
    public sBy: number = 0;

    public mass: number = 0;
    public motorMass: number = 0;
    public axialMass: number = 0;
    public springMass: number = 0;

    public bias: number = 0;
    public gamma: number = 0;

    public readonly qA: Rot = new Rot();
    public readonly qB: Rot = new Rot();
    public readonly lalcA: Vec2 = new Vec2();
    public readonly lalcB: Vec2 = new Vec2();
    public readonly rA: Vec2 = new Vec2();
    public readonly rB: Vec2 = new Vec2();

    constructor(def: IWheelJointDef) {
      super(def);

      this.localAnchorA.copy(maybe(def.localAnchorA, Vec2.ZERO));
      this.localAnchorB.copy(maybe(def.localAnchorB, Vec2.ZERO));
      this.localXAxisA.copy(maybe(def.localAxisA, Vec2.UNITX));
      Vec2.CrossOneV(this.localXAxisA, this.localYAxisA);

      this.lowerTranslation = maybe(def.lowerTranslation, 0);
      this.upperTranslation = maybe(def.upperTranslation, 0);
      this.enableLimit = maybe(def.enableLimit, false);

      this.maxMotorTorque = maybe(def.maxMotorTorque, 0);
      this.motorSpeed = maybe(def.motorSpeed, 0);
      this.enableMotor = maybe(def.enableMotor, false);

      this.ax.setZero();
      this.ay.setZero();

      this.stiffness = maybe(def.stiffness, 0);
      this.damping = maybe(def.damping, 0);
    }

    public getMotorSpeed(): number {
      return this.motorSpeed;
    }

    public getMaxMotorTorque(): number {
      return this.maxMotorTorque;
    }

    public setSpringFrequencyHz(hz: number): void {
      this.stiffness = hz;
    }

    public getSpringFrequencyHz(): number {
      return this.stiffness;
    }

    public setSpringDampingRatio(ratio: number): void {
      this.damping = ratio;
    }

    public getSpringDampingRatio(): number {
      return this.damping;
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

      const mA: number = this.invMassA, mB: number = this.invMassB;
      const iA: number = this.invIA, iB: number = this.invIB;

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
      // Vec2 d = cB + rB - cA - rA;
      const d: Vec2 = Vec2.SubVV(
        Vec2.AddVV(cB, rB, Vec2.s_t0),
        Vec2.AddVV(cA, rA, Vec2.s_t1),
        WheelJoint.initVelocityConstraints_s_d);

      // Point to line constraint
      {
        // ay = Mul(qA, localYAxisA);
        Rot.mulRV(qA, this.localYAxisA, this.ay);
        // sAy = Cross(d + rA, ay);
        this.sAy = Vec2.CrossVV(Vec2.AddVV(d, rA, Vec2.s_t0), this.ay);
        // sBy = Cross(rB, ay);
        this.sBy = Vec2.CrossVV(rB, this.ay);

        this.mass = mA + mB + iA * this.sAy * this.sAy + iB * this.sBy * this.sBy;

        if (this.mass > 0) {
          this.mass = 1 / this.mass;
        }
      }

      // Spring constraint
      Rot.mulRV(qA, this.localXAxisA, this.ax); // ax = Mul(qA, localXAxisA);
      this.sAx = Vec2.CrossVV(Vec2.AddVV(d, rA, Vec2.s_t0), this.ax);
      this.sBx = Vec2.CrossVV(rB, this.ax);

      const invMass: number = mA + mB + iA * this.sAx * this.sAx + iB * this.sBx * this.sBx;
      if (invMass > 0.0) {
        this.axialMass = 1.0 / invMass;
      } else {
        this.axialMass = 0.0;
      }

      this.springMass = 0;
      this.bias = 0;
      this.gamma = 0;

      if (this.stiffness > 0.0 && invMass > 0.0) {
        this.springMass = 1.0 / invMass;

        const C: number = Vec2.DotVV(d, this.ax);

        // magic formulas
        const h: number = data.step.dt;
        this.gamma = h * (this.damping + h * this.stiffness);
        if (this.gamma > 0.0) {
          this.gamma = 1.0 / this.gamma;
        }

        this.bias = C * h * this.stiffness * this.gamma;

        this.springMass = invMass + this.gamma;
        if (this.springMass > 0.0) {
          this.springMass = 1.0 / this.springMass;
        }
      } else {
        this.springImpulse = 0.0;
      }

      if (this.enableLimit) {
        this.translation = Vec2.DotVV(this.ax, d);
      } else {
        this.lowerImpulse = 0.0;
        this.upperImpulse = 0.0;
      }

      if (this.enableMotor) {
        this.motorMass = iA + iB;
        if (this.motorMass > 0) {
          this.motorMass = 1 / this.motorMass;
        }
      } else {
        this.motorMass = 0;
        this.motorImpulse = 0;
      }

      if (data.step.warmStarting) {
        // Account for variable time step.
        this.impulse *= data.step.dtRatio;
        this.springImpulse *= data.step.dtRatio;
        this.motorImpulse *= data.step.dtRatio;

        const axialImpulse: number = this.springImpulse + this.lowerImpulse - this.upperImpulse;
        // Vec2 P = impulse * ay + springImpulse * ax;
        const P: Vec2 = Vec2.AddVV(
          Vec2.MulSV(this.impulse, this.ay, Vec2.s_t0),
          Vec2.MulSV(axialImpulse, this.ax, Vec2.s_t1),
          WheelJoint.initVelocityConstraints_s_P);
        // float32 LA = impulse * sAy + springImpulse * sAx + motorImpulse;
        const LA: number = this.impulse * this.sAy + axialImpulse * this.sAx + this.motorImpulse;
        // float32 LB = impulse * sBy + springImpulse * sBx + motorImpulse;
        const LB: number = this.impulse * this.sBy + axialImpulse * this.sBx + this.motorImpulse;

        // vA -= invMassA * P;
        vA.selfMulSub(this.invMassA, P);
        wA -= this.invIA * LA;

        // vB += invMassB * P;
        vB.selfMulAdd(this.invMassB, P);
        wB += this.invIB * LB;
      } else {
        this.impulse = 0;
        this.springImpulse = 0;
        this.motorImpulse = 0;
        this.lowerImpulse = 0;
        this.upperImpulse = 0;
      }

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    private static solveVelocityConstraints_s_P = new Vec2();
    public solveVelocityConstraints(data: SolverData): void {
      const mA: number = this.invMassA, mB: number = this.invMassB;
      const iA: number = this.invIA, iB: number = this.invIB;

      const vA: Vec2 = data.velocities[this.indexA].v;
      let wA: number = data.velocities[this.indexA].w;
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;

      // Solve spring constraint
      {
        const Cdot: number = Vec2.DotVV(this.ax, Vec2.SubVV(vB, vA, Vec2.s_t0)) + this.sBx * wB - this.sAx * wA;
        const impulse: number = -this.springMass * (Cdot + this.bias + this.gamma * this.springImpulse);
        this.springImpulse += impulse;

        // Vec2 P = impulse * ax;
        const P: Vec2 = Vec2.MulSV(impulse, this.ax, WheelJoint.solveVelocityConstraints_s_P);
        const LA: number = impulse * this.sAx;
        const LB: number = impulse * this.sBx;

        // vA -= mA * P;
        vA.selfMulSub(mA, P);
        wA -= iA * LA;

        // vB += mB * P;
        vB.selfMulAdd(mB, P);
        wB += iB * LB;
      }

      // Solve rotational motor constraint
      {
        const Cdot: number = wB - wA - this.motorSpeed;
        let impulse: number = -this.motorMass * Cdot;

        const oldImpulse: number = this.motorImpulse;
        const maxImpulse: number = data.step.dt * this.maxMotorTorque;
        this.motorImpulse = clamp(this.motorImpulse + impulse, -maxImpulse, maxImpulse);
        impulse = this.motorImpulse - oldImpulse;

        wA -= iA * impulse;
        wB += iB * impulse;
      }

      if (this.enableLimit) {
        // Lower limit
        {
          const C: number = this.translation - this.lowerTranslation;
          const Cdot: number = Vec2.DotVV(this.ax, Vec2.SubVV(vB, vA, Vec2.s_t0)) + this.sBx * wB - this.sAx * wA;
          let impulse: number = -this.axialMass * (Cdot + Max(C, 0.0) * data.step.inv_dt);
          const oldImpulse: number = this.lowerImpulse;
          this.lowerImpulse = Max(this.lowerImpulse + impulse, 0.0);
          impulse = this.lowerImpulse - oldImpulse;

          // Vec2 P = impulse * this.ax;
          const P: Vec2 = Vec2.MulSV(impulse, this.ax, WheelJoint.solveVelocityConstraints_s_P);
          const LA: number = impulse * this.sAx;
          const LB: number = impulse * this.sBx;

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
          const Cdot: number = Vec2.DotVV(this.ax, Vec2.SubVV(vA, vB, Vec2.s_t0)) + this.sAx * wA - this.sBx * wB;
          let impulse: number = -this.axialMass * (Cdot + Max(C, 0.0) * data.step.inv_dt);
          const oldImpulse: number = this.upperImpulse;
          this.upperImpulse = Max(this.upperImpulse + impulse, 0.0);
          impulse = this.upperImpulse - oldImpulse;

          // Vec2 P = impulse * this.ax;
          const P: Vec2 = Vec2.MulSV(impulse, this.ax, WheelJoint.solveVelocityConstraints_s_P);
          const LA: number = impulse * this.sAx;
          const LB: number = impulse * this.sBx;

          // vA += mA * P;
          vA.selfMulAdd(mA, P);
          wA += iA * LA;
          // vB -= mB * P;
          vB.selfMulSub(mB, P);
          wB -= iB * LB;
        }
      }

      // Solve point to line constraint
      {
        const Cdot: number = Vec2.DotVV(this.ay, Vec2.SubVV(vB, vA, Vec2.s_t0)) + this.sBy * wB - this.sAy * wA;
        const impulse: number = -this.mass * Cdot;
        this.impulse += impulse;

        // Vec2 P = impulse * ay;
        const P: Vec2 = Vec2.MulSV(impulse, this.ay, WheelJoint.solveVelocityConstraints_s_P);
        const LA: number = impulse * this.sAy;
        const LB: number = impulse * this.sBy;

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

    private static solvePositionConstraints_s_d = new Vec2();
    private static solvePositionConstraints_s_P = new Vec2();
    public solvePositionConstraints(data: SolverData): boolean {
      const cA: Vec2 = data.positions[this.indexA].c;
      let aA: number = data.positions[this.indexA].a;
      const cB: Vec2 = data.positions[this.indexB].c;
      let aB: number = data.positions[this.indexB].a;

      // const qA: Rot = this.qA.SetAngle(aA), qB: Rot = this.qB.SetAngle(aB);

      // // Vec2 rA = Mul(qA, localAnchorA - localCenterA);
      // Vec2.SubVV(this.localAnchorA, this.localCenterA, this.lalcA);
      // const rA: Vec2 = Rot.MulRV(qA, this.lalcA, this.rA);
      // // Vec2 rB = Mul(qB, localAnchorB - localCenterB);
      // Vec2.SubVV(this.localAnchorB, this.localCenterB, this.lalcB);
      // const rB: Vec2 = Rot.MulRV(qB, this.lalcB, this.rB);
      // // Vec2 d = (cB - cA) + rB - rA;
      // const d: Vec2 = Vec2.AddVV(
      //   Vec2.SubVV(cB, cA, Vec2.s_t0),
      //   Vec2.SubVV(rB, rA, Vec2.s_t1),
      //   WheelJoint.solvePositionConstraints_s_d);

      // // Vec2 ay = Mul(qA, localYAxisA);
      // const ay: Vec2 = Rot.MulRV(qA, this.localYAxisA, this.ay);

      // // float32 sAy = Cross(d + rA, ay);
      // const sAy = Vec2.CrossVV(Vec2.AddVV(d, rA, Vec2.s_t0), ay);
      // // float32 sBy = Cross(rB, ay);
      // const sBy = Vec2.CrossVV(rB, ay);

      // // float32 C = Dot(d, ay);
      // const C: number = Vec2.DotVV(d, this.ay);

      // const k: number = this.invMassA + this.invMassB + this.invIA * this.sAy * this.sAy + this.invIB * this.sBy * this.sBy;

      // let impulse: number;
      // if (k !== 0) {
      //   impulse = - C / k;
      // } else {
      //   impulse = 0;
      // }

      // // Vec2 P = impulse * ay;
      // const P: Vec2 = Vec2.MulSV(impulse, ay, WheelJoint.solvePositionConstraints_s_P);
      // const LA: number = impulse * sAy;
      // const LB: number = impulse * sBy;

      // // cA -= invMassA * P;
      // cA.SelfMulSub(this.invMassA, P);
      // aA -= this.invIA * LA;
      // // cB += invMassB * P;
      // cB.SelfMulAdd(this.invMassB, P);
      // aB += this.invIB * LB;

      let linearError: number = 0.0;

      if (this.enableLimit) {
        // Rot qA(aA), qB(aB);
        const qA: Rot = this.qA.setAngle(aA), qB: Rot = this.qB.setAngle(aB);

        // Vec2 rA = Mul(qA, this.localAnchorA - this.localCenterA);
        // Vec2 rB = Mul(qB, this.localAnchorB - this.localCenterB);
        // Vec2 d = (cB - cA) + rB - rA;

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
          WheelJoint.solvePositionConstraints_s_d);

        // Vec2 ax = Mul(qA, this.localXAxisA);
        const ax: Vec2 = Rot.mulRV(qA, this.localXAxisA, this.ax);
        // float sAx = Cross(d + rA, this.ax);
        const sAx = Vec2.CrossVV(Vec2.AddVV(d, rA, Vec2.s_t0), this.ax);
        // float sBx = Cross(rB, this.ax);
        const sBx = Vec2.CrossVV(rB, this.ax);

        let C: number = 0.0;
        const translation: number = Vec2.DotVV(ax, d);
        if (Abs(this.upperTranslation - this.lowerTranslation) < 2.0 * linearSlop) {
          C = translation;
        } else if (translation <= this.lowerTranslation) {
          C = Min(translation - this.lowerTranslation, 0.0);
        } else if (translation >= this.upperTranslation) {
          C = Max(translation - this.upperTranslation, 0.0);
        }

        if (C !== 0.0) {

          const invMass: number = this.invMassA + this.invMassB + this.invIA * sAx * sAx + this.invIB * sBx * sBx;
          let impulse: number = 0.0;
          if (invMass !== 0.0) {
            impulse = -C / invMass;
          }

          const P: Vec2 = Vec2.MulSV(impulse, ax, WheelJoint.solvePositionConstraints_s_P);
          const LA: number = impulse * sAx;
          const LB: number = impulse * sBx;

          // cA -= invMassA * P;
          cA.selfMulSub(this.invMassA, P);
          aA -= this.invIA * LA;
          // cB += invMassB * P;
          cB.selfMulAdd(this.invMassB, P);
          // aB += invIB * LB;
          aB += this.invIB * LB;

          linearError = Abs(C);
        }
      }

      // Solve perpendicular constraint
      {
        // Rot qA(aA), qB(aB);
        const qA: Rot = this.qA.setAngle(aA), qB: Rot = this.qB.setAngle(aB);

        // Vec2 rA = Mul(qA, localAnchorA - localCenterA);
        // Vec2 rB = Mul(qB, localAnchorB - localCenterB);
        // Vec2 d = (cB - cA) + rB - rA;

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
          WheelJoint.solvePositionConstraints_s_d);

        // Vec2 ay = Mul(qA, localYAxisA);
        const ay: Vec2 = Rot.mulRV(qA, this.localYAxisA, this.ay);

        // float sAy = Cross(d + rA, ay);
        const sAy = Vec2.CrossVV(Vec2.AddVV(d, rA, Vec2.s_t0), ay);
        // float sBy = Cross(rB, ay);
        const sBy = Vec2.CrossVV(rB, ay);

        // float C = Dot(d, ay);
        const C: number = Vec2.DotVV(d, ay);

        const invMass: number = this.invMassA + this.invMassB + this.invIA * this.sAy * this.sAy + this.invIB * this.sBy * this.sBy;

        let impulse: number = 0.0;
        if (invMass !== 0.0) {
          impulse = - C / invMass;
        }

        // Vec2 P = impulse * ay;
        // const LA: number = impulse * sAy;
        // const LB: number = impulse * sBy;
        const P: Vec2 = Vec2.MulSV(impulse, ay, WheelJoint.solvePositionConstraints_s_P);
        const LA: number = impulse * sAy;
        const LB: number = impulse * sBy;

        // cA -= invMassA * P;
        cA.selfMulSub(this.invMassA, P);
        aA -= this.invIA * LA;
        // cB += invMassB * P;
        cB.selfMulAdd(this.invMassB, P);
        aB += this.invIB * LB;

        linearError = Max(linearError, Abs(C));
      }

      // data.positions[this.indexA].c = cA;
      data.positions[this.indexA].a = aA;
      // data.positions[this.indexB].c = cB;
      data.positions[this.indexB].a = aB;

      return linearError <= linearSlop;
    }

    public getDefinition(def: WheelJointDef): WheelJointDef {
      // DEBUG: Assert(false); // TODO
      return def;
    }

    public getAnchorA<T extends XY>(out: T): T {
      return this.bodyA.getWorldPoint(this.localAnchorA, out);
    }

    public getAnchorB<T extends XY>(out: T): T {
      return this.bodyB.getWorldPoint(this.localAnchorB, out);
    }

    public getReactionForce<T extends XY>(inv_dt: number, out: T): T {
      out.x = inv_dt * (this.impulse * this.ay.x + (this.springImpulse + this.lowerImpulse - this.upperImpulse) * this.ax.x);
      out.y = inv_dt * (this.impulse * this.ay.y + (this.springImpulse + this.lowerImpulse - this.upperImpulse) * this.ax.y);
      return out;
    }

    public getReactionTorque(inv_dt: number): number {
      return inv_dt * this.motorImpulse;
    }

    public getJointTranslation(): number {
      return this.getPrismaticJointTranslation();
    }

    public getJointLinearSpeed(): number {
      return this.getPrismaticJointSpeed();
    }

    public getJointAngle(): number {
      return this.getRevoluteJointAngle();
    }

    public getJointAngularSpeed(): number {
      return this.getRevoluteJointSpeed();
    }

    public getPrismaticJointTranslation(): number {
      const bA: Body = this.bodyA;
      const bB: Body = this.bodyB;

      const pA: Vec2 = bA.getWorldPoint(this.localAnchorA, new Vec2());
      const pB: Vec2 = bB.getWorldPoint(this.localAnchorB, new Vec2());
      const d: Vec2 = Vec2.SubVV(pB, pA, new Vec2());
      const axis: Vec2 = bA.getWorldVector(this.localXAxisA, new Vec2());

      const translation: number = Vec2.DotVV(d, axis);
      return translation;
    }

    public getPrismaticJointSpeed(): number {
      const bA: Body = this.bodyA;
      const bB: Body = this.bodyB;

      // Vec2 rA = Mul(bA.xf.q, localAnchorA - bA.sweep.localCenter);
      Vec2.SubVV(this.localAnchorA, bA.sweep.localCenter, this.lalcA);
      const rA = Rot.mulRV(bA.xf.q, this.lalcA, this.rA);
      // Vec2 rB = Mul(bB.xf.q, localAnchorB - bB.sweep.localCenter);
      Vec2.SubVV(this.localAnchorB, bB.sweep.localCenter, this.lalcB);
      const rB = Rot.mulRV(bB.xf.q, this.lalcB, this.rB);
      // Vec2 pA = bA.sweep.c + rA;
      const pA = Vec2.AddVV(bA.sweep.c, rA, Vec2.s_t0); // pA uses s_t0
      // Vec2 pB = bB.sweep.c + rB;
      const pB = Vec2.AddVV(bB.sweep.c, rB, Vec2.s_t1); // pB uses s_t1
      // Vec2 d = pB - pA;
      const d = Vec2.SubVV(pB, pA, Vec2.s_t2); // d uses s_t2
      // Vec2 axis = Mul(bA.xf.q, localXAxisA);
      const axis = bA.getWorldVector(this.localXAxisA, new Vec2());

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

    public getRevoluteJointAngle(): number {
      // Body* bA = this.bodyA;
      // Body* bB = this.bodyB;
      // return bB.this.sweep.a - bA.this.sweep.a;
      return this.bodyB.sweep.a - this.bodyA.sweep.a;
    }

    public getRevoluteJointSpeed(): number {
      const wA: number = this.bodyA.angularVelocity;
      const wB: number = this.bodyB.angularVelocity;
      return wB - wA;
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

    public setMaxMotorTorque(force: number): void {
      if (force !== this.maxMotorTorque) {
        this.bodyA.setAwake(true);
        this.bodyB.setAwake(true);
        this.maxMotorTorque = force;
      }
    }

    public getMotorTorque(inv_dt: number): number {
      return inv_dt * this.motorImpulse;
    }

    /// Is the joint limit enabled?
    public isLimitEnabled(): boolean {
      return this.enableLimit;
    }

    /// Enable/disable the joint translation limit.
    public setEnableLimit(flag: boolean): void {
      if (flag !== this.enableLimit) {
        this.bodyA.setAwake(true);
        this.bodyB.setAwake(true);
        this.enableLimit = flag;
        this.lowerImpulse = 0.0;
        this.upperImpulse = 0.0;
      }
    }

    /// Get the lower joint translation limit, usually in meters.
    public getLowerLimit(): number {
      return this.lowerTranslation;
    }

    /// Get the upper joint translation limit, usually in meters.
    public getUpperLimit(): number {
      return this.upperTranslation;
    }

    /// Set the joint translation limits, usually in meters.
    public setLimits(lower: number, upper: number): void {
      // Assert(lower <= upper);
      if (lower !== this.lowerTranslation || upper !== this.upperTranslation) {
        this.bodyA.setAwake(true);
        this.bodyB.setAwake(true);
        this.lowerTranslation = lower;
        this.upperTranslation = upper;
        this.lowerImpulse = 0.0;
        this.upperImpulse = 0.0;
      }
    }

    public dump(log: (format: string, ...args: any[]) => void): void {
      const indexA = this.bodyA.islandIndex;
      const indexB = this.bodyB.islandIndex;

      log("  const jd: WheelJointDef = new WheelJointDef();\n");
      log("  jd.bodyA = bodies[%d];\n", indexA);
      log("  jd.bodyB = bodies[%d];\n", indexB);
      log("  jd.collideConnected = %s;\n", (this.collideConnected) ? ("true") : ("false"));
      log("  jd.localAnchorA.Set(%.15f, %.15f);\n", this.localAnchorA.x, this.localAnchorA.y);
      log("  jd.localAnchorB.Set(%.15f, %.15f);\n", this.localAnchorB.x, this.localAnchorB.y);
      log("  jd.localAxisA.Set(%.15f, %.15f);\n", this.localXAxisA.x, this.localXAxisA.y);
      log("  jd.enableMotor = %s;\n", (this.enableMotor) ? ("true") : ("false"));
      log("  jd.motorSpeed = %.15f;\n", this.motorSpeed);
      log("  jd.maxMotorTorque = %.15f;\n", this.maxMotorTorque);
      log("  jd.stiffness = %.15f;\n", this.stiffness);
      log("  jd.damping = %.15f;\n", this.damping);
      log("  joints[%d] = this.world.CreateJoint(jd);\n", this.index);
    }

    ///
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
      const pA = Transform.mulXV(xfA, this.localAnchorA, WheelJoint.draw_s_pA);
      const pB = Transform.mulXV(xfB, this.localAnchorB, WheelJoint.draw_s_pB);

      // Vec2 axis = Mul(xfA.q, localXAxisA);
      const axis: Vec2 = Rot.mulRV(xfA.q, this.localXAxisA, WheelJoint.draw_s_axis);

      const c1 = WheelJoint.draw_s_c1; // Color c1(0.7f, 0.7f, 0.7f);
      const c2 = WheelJoint.draw_s_c2; // Color c2(0.3f, 0.9f, 0.3f);
      const c3 = WheelJoint.draw_s_c3; // Color c3(0.9f, 0.3f, 0.3f);
      const c4 = WheelJoint.draw_s_c4; // Color c4(0.3f, 0.3f, 0.9f);
      const c5 = WheelJoint.draw_s_c5; // Color c5(0.4f, 0.4f, 0.4f);

      draw.drawSegment(pA, pB, c5);

      if (this.enableLimit) {
        // Vec2 lower = pA + lowerTranslation * axis;
        const lower = Vec2.AddVMulSV(pA, this.lowerTranslation, axis, WheelJoint.draw_s_lower);
        // Vec2 upper = pA + upperTranslation * axis;
        const upper = Vec2.AddVMulSV(pA, this.upperTranslation, axis, WheelJoint.draw_s_upper);
        // Vec2 perp = Mul(xfA.q, localYAxisA);
        const perp = Rot.mulRV(xfA.q, this.localYAxisA, WheelJoint.draw_s_perp);
        // draw.DrawSegment(lower, upper, c1);
        draw.drawSegment(lower, upper, c1);
        // draw.DrawSegment(lower - 0.5f * perp, lower + 0.5f * perp, c2);
        draw.drawSegment(Vec2.AddVMulSV(lower, -0.5, perp, Vec2.s_t0), Vec2.AddVMulSV(lower, 0.5, perp, Vec2.s_t1), c2);
        // draw.DrawSegment(upper - 0.5f * perp, upper + 0.5f * perp, c3);
        draw.drawSegment(Vec2.AddVMulSV(upper, -0.5, perp, Vec2.s_t0), Vec2.AddVMulSV(upper, 0.5, perp, Vec2.s_t1), c3);
      } else {
        // draw.DrawSegment(pA - 1.0f * axis, pA + 1.0f * axis, c1);
        draw.drawSegment(Vec2.AddVMulSV(pA, -1.0, axis, Vec2.s_t0), Vec2.AddVMulSV(pA, 1.0, axis, Vec2.s_t1), c1);
      }

      draw.drawPoint(pA, 5.0, c1);
      draw.drawPoint(pB, 5.0, c4);
    }
  }

}
