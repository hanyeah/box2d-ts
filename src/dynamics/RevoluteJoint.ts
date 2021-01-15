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
  export interface IRevoluteJointDef extends IJointDef {
    localAnchorA?: XY;

    localAnchorB?: XY;

    referenceAngle?: number;

    enableLimit?: boolean;

    lowerAngle?: number;

    upperAngle?: number;

    enableMotor?: boolean;

    motorSpeed?: number;

    maxMotorTorque?: number;
  }

/// Revolute joint definition. This requires defining an anchor point where the
/// bodies are joined. The definition uses local anchor points so that the
/// initial configuration can violate the constraint slightly. You also need to
/// specify the initial relative angle for joint limits. This helps when saving
/// and loading a game.
/// The local anchor points are measured from the body's origin
/// rather than the center of mass because:
/// 1. you might not know where the center of mass will be.
/// 2. if you add/remove shapes from a body and recompute the mass,
///    the joints will be broken.
  export class RevoluteJointDef extends JointDef implements IRevoluteJointDef {
    public readonly localAnchorA: Vec2 = new Vec2(0, 0);

    public readonly localAnchorB: Vec2 = new Vec2(0, 0);

    public referenceAngle: number = 0;

    public enableLimit = false;

    public lowerAngle: number = 0;

    public upperAngle: number = 0;

    public enableMotor = false;

    public motorSpeed: number = 0;

    public maxMotorTorque: number = 0;

    constructor() {
      super(JointType.RevoluteJoint);
    }

    public initialize(bA: Body, bB: Body, anchor: XY): void {
      this.bodyA = bA;
      this.bodyB = bB;
      this.bodyA.getLocalPoint(anchor, this.localAnchorA);
      this.bodyB.getLocalPoint(anchor, this.localAnchorB);
      this.referenceAngle = this.bodyB.getAngle() - this.bodyA.getAngle();
    }
  }

  export class RevoluteJoint extends Joint {
    // Solver shared
    public readonly localAnchorA: Vec2 = new Vec2();
    public readonly localAnchorB: Vec2 = new Vec2();
    public readonly impulse: Vec2 = new Vec2();
    public motorImpulse: number = 0;
    public lowerImpulse: number = 0;
    public upperImpulse: number = 0;
    public enableMotor: boolean = false;
    public maxMotorTorque: number = 0;
    public motorSpeed: number = 0;
    public enableLimit: boolean = false;
    public referenceAngle: number = 0;
    public lowerAngle: number = 0;
    public upperAngle: number = 0;

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
    public readonly K: Mat22 = new Mat22();
    public angle: number = 0;
    public axialMass: number = 0;

    public readonly qA: Rot = new Rot();
    public readonly qB: Rot = new Rot();
    public readonly lalcA: Vec2 = new Vec2();
    public readonly lalcB: Vec2 = new Vec2();

    constructor(def: IRevoluteJointDef) {
      super(def);

      this.localAnchorA.copy(maybe(def.localAnchorA, Vec2.ZERO));
      this.localAnchorB.copy(maybe(def.localAnchorB, Vec2.ZERO));
      this.referenceAngle = maybe(def.referenceAngle, 0);

      this.impulse.setZero();
      this.motorImpulse = 0;

      this.lowerAngle = maybe(def.lowerAngle, 0);
      this.upperAngle = maybe(def.upperAngle, 0);
      this.maxMotorTorque = maybe(def.maxMotorTorque, 0);
      this.motorSpeed = maybe(def.motorSpeed, 0);
      this.enableLimit = maybe(def.enableLimit, false);
      this.enableMotor = maybe(def.enableMotor, false);
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

      // Rot qA(aA), qB(aB);
      const qA: Rot = this.qA.setAngle(aA), qB: Rot = this.qB.setAngle(aB);

      // rA = Mul(qA, localAnchorA - localCenterA);
      Vec2.SubVV(this.localAnchorA, this.localCenterA, this.lalcA);
      Rot.mulRV(qA, this.lalcA, this.rA);
      // rB = Mul(qB, localAnchorB - localCenterB);
      Vec2.SubVV(this.localAnchorB, this.localCenterB, this.lalcB);
      Rot.mulRV(qB, this.lalcB, this.rB);

      // J = [-I -r1_skew I r2_skew]
      // r_skew = [-ry; rx]

      // Matlab
      // K = [ mA+r1y^2*iA+mB+r2y^2*iB,  -r1y*iA*r1x-r2y*iB*r2x]
      //     [  -r1y*iA*r1x-r2y*iB*r2x, mA+r1x^2*iA+mB+r2x^2*iB]

      const mA: number = this.invMassA, mB: number = this.invMassB;
      const iA: number = this.invIA, iB: number = this.invIB;

      this.K.ex.x = mA + mB + this.rA.y * this.rA.y * iA + this.rB.y * this.rB.y * iB;
      this.K.ey.x = -this.rA.y * this.rA.x * iA - this.rB.y * this.rB.x * iB;
      this.K.ex.y = this.K.ey.x;
      this.K.ey.y = mA + mB + this.rA.x * this.rA.x * iA + this.rB.x * this.rB.x * iB;

      this.axialMass = iA + iB;
      let fixedRotation: boolean;
      if (this.axialMass > 0.0) {
        this.axialMass = 1.0 / this.axialMass;
        fixedRotation = false;
      } else {
        fixedRotation = true;
      }

      this.angle = aB - aA - this.referenceAngle;
      if (this.enableLimit === false || fixedRotation) {
        this.lowerImpulse = 0.0;
        this.upperImpulse = 0.0;
      }

      if (this.enableMotor === false || fixedRotation) {
        this.motorImpulse = 0.0;
      }

      if (data.step.warmStarting) {
        // Scale impulses to support a variable time step.
        this.impulse.selfMul(data.step.dtRatio);
        this.motorImpulse *= data.step.dtRatio;
        this.lowerImpulse *= data.step.dtRatio;
        this.upperImpulse *= data.step.dtRatio;

        const axialImpulse: number = this.motorImpulse + this.lowerImpulse - this.upperImpulse;
        // Vec2 P(impulse.x, impulse.y);
        const P: Vec2 = RevoluteJoint.initVelocityConstraints_s_P.set(this.impulse.x, this.impulse.y);

        // vA -= mA * P;
        vA.selfMulSub(mA, P);
        wA -= iA * (Vec2.CrossVV(this.rA, P) + axialImpulse);

        // vB += mB * P;
        vB.selfMulAdd(mB, P);
        wB += iB * (Vec2.CrossVV(this.rB, P) + axialImpulse);
      } else {
        this.impulse.setZero();
        this.motorImpulse = 0;
        this.lowerImpulse = 0;
        this.upperImpulse = 0;
      }

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    // private static solveVelocityConstraints_s_P: Vec2 = new Vec2();
    private static solveVelocityConstraints_s_Cdot_v2: Vec2 = new Vec2();
    // private static solveVelocityConstraints_s_Cdot1: Vec2 = new Vec2();
    // private static solveVelocityConstraints_s_impulse_v3: Vec3 = new Vec3();
    // private static solveVelocityConstraints_s_reduced_v2: Vec2 = new Vec2();
    private static solveVelocityConstraints_s_impulse_v2: Vec2 = new Vec2();
    public solveVelocityConstraints(data: SolverData): void {
      const vA: Vec2 = data.velocities[this.indexA].v;
      let wA: number = data.velocities[this.indexA].w;
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;

      const mA: number = this.invMassA, mB: number = this.invMassB;
      const iA: number = this.invIA, iB: number = this.invIB;

      const fixedRotation: boolean = (iA + iB === 0);

      // Solve motor constraint.
      if (this.enableMotor && !fixedRotation) {
        const Cdot: number = wB - wA - this.motorSpeed;
        let impulse: number = -this.axialMass * Cdot;
        const oldImpulse: number = this.motorImpulse;
        const maxImpulse: number = data.step.dt * this.maxMotorTorque;
        this.motorImpulse = clamp(this.motorImpulse + impulse, -maxImpulse, maxImpulse);
        impulse = this.motorImpulse - oldImpulse;

        wA -= iA * impulse;
        wB += iB * impulse;
      }

      // Solve limit constraint.
      if (this.enableLimit && !fixedRotation) {
        // Lower limit
        {
          const C: number = this.angle - this.lowerAngle;
          const Cdot: number = wB - wA;
          let impulse: number = -this.axialMass * (Cdot + Max(C, 0.0) * data.step.inv_dt);
          const oldImpulse: number = this.lowerImpulse;
          this.lowerImpulse = Max(this.lowerImpulse + impulse, 0.0);
          impulse = this.lowerImpulse - oldImpulse;

          wA -= iA * impulse;
          wB += iB * impulse;
        }

        // Upper limit
        // Note: signs are flipped to keep C positive when the constraint is satisfied.
        // This also keeps the impulse positive when the limit is active.
        {
          const C: number = this.upperAngle - this.angle;
          const Cdot: number = wA - wB;
          let impulse: number = -this.axialMass * (Cdot + Max(C, 0.0) * data.step.inv_dt);
          const oldImpulse: number = this.upperImpulse;
          this.upperImpulse = Max(this.upperImpulse + impulse, 0.0);
          impulse = this.upperImpulse - oldImpulse;

          wA += iA * impulse;
          wB -= iB * impulse;
        }
      }

      // Solve point-to-point constraint
      {
        // Vec2 Cdot = vB + Cross(wB, rB) - vA - Cross(wA, rA);
        const Cdot_v2: Vec2 = Vec2.SubVV(
          Vec2.AddVCrossSV(vB, wB, this.rB, Vec2.s_t0),
          Vec2.AddVCrossSV(vA, wA, this.rA, Vec2.s_t1),
          RevoluteJoint.solveVelocityConstraints_s_Cdot_v2);
        // Vec2 impulse = K.Solve(-Cdot);
        const impulse_v2: Vec2 = this.K.solve(-Cdot_v2.x, -Cdot_v2.y, RevoluteJoint.solveVelocityConstraints_s_impulse_v2);

        this.impulse.x += impulse_v2.x;
        this.impulse.y += impulse_v2.y;

        // vA -= mA * impulse;
        vA.selfMulSub(mA, impulse_v2);
        wA -= iA * Vec2.CrossVV(this.rA, impulse_v2);

        // vB += mB * impulse;
        vB.selfMulAdd(mB, impulse_v2);
        wB += iB * Vec2.CrossVV(this.rB, impulse_v2);
      }

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    private static solvePositionConstraints_s_C_v2 = new Vec2();
    private static solvePositionConstraints_s_impulse = new Vec2();
    public solvePositionConstraints(data: SolverData): boolean {
      const cA: Vec2 = data.positions[this.indexA].c;
      let aA: number = data.positions[this.indexA].a;
      const cB: Vec2 = data.positions[this.indexB].c;
      let aB: number = data.positions[this.indexB].a;

      // Rot qA(aA), qB(aB);
      const qA: Rot = this.qA.setAngle(aA), qB: Rot = this.qB.setAngle(aB);

      let angularError: number = 0;
      let positionError: number = 0;

      const fixedRotation: boolean = (this.invIA + this.invIB === 0);

      // Solve angular limit constraint.
      if (this.enableLimit && !fixedRotation) {
        const angle: number = aB - aA - this.referenceAngle;
        let C: number = 0.0;

        if (Abs(this.upperAngle - this.lowerAngle) < 2.0 * angularSlop) {
          // Prevent large angular corrections
          C = clamp(angle - this.lowerAngle, -maxAngularCorrection, maxAngularCorrection);
        } else if (angle <= this.lowerAngle) {
          // Prevent large angular corrections and allow some slop.
          C = clamp(angle - this.lowerAngle + angularSlop, -maxAngularCorrection, 0.0);
        } else if (angle >= this.upperAngle) {
          // Prevent large angular corrections and allow some slop.
          C = clamp(angle - this.upperAngle - angularSlop, 0.0, maxAngularCorrection);
        }

        const limitImpulse: number = -this.axialMass * C;
        aA -= this.invIA * limitImpulse;
        aB += this.invIB * limitImpulse;
        angularError = Abs(C);
      }

      // Solve point-to-point constraint.
      {
        qA.setAngle(aA);
        qB.setAngle(aB);
        // Vec2 rA = Mul(qA, localAnchorA - localCenterA);
        Vec2.SubVV(this.localAnchorA, this.localCenterA, this.lalcA);
        const rA: Vec2 = Rot.mulRV(qA, this.lalcA, this.rA);
        // Vec2 rB = Mul(qB, localAnchorB - localCenterB);
        Vec2.SubVV(this.localAnchorB, this.localCenterB, this.lalcB);
        const rB: Vec2 = Rot.mulRV(qB, this.lalcB, this.rB);

        // Vec2 C = cB + rB - cA - rA;
        const C_v2 =
          Vec2.SubVV(
            Vec2.AddVV(cB, rB, Vec2.s_t0),
            Vec2.AddVV(cA, rA, Vec2.s_t1),
            RevoluteJoint.solvePositionConstraints_s_C_v2);
        // positionError = C.Length();
        positionError = C_v2.length();

        const mA: number = this.invMassA, mB: number = this.invMassB;
        const iA: number = this.invIA, iB: number = this.invIB;

        const K: Mat22 = this.K;
        K.ex.x = mA + mB + iA * rA.y * rA.y + iB * rB.y * rB.y;
        K.ex.y = -iA * rA.x * rA.y - iB * rB.x * rB.y;
        K.ey.x = K.ex.y;
        K.ey.y = mA + mB + iA * rA.x * rA.x + iB * rB.x * rB.x;

        // Vec2 impulse = -K.Solve(C);
        const impulse: Vec2 = K.solve(C_v2.x, C_v2.y, RevoluteJoint.solvePositionConstraints_s_impulse).selfNeg();

        // cA -= mA * impulse;
        cA.selfMulSub(mA, impulse);
        aA -= iA * Vec2.CrossVV(rA, impulse);

        // cB += mB * impulse;
        cB.selfMulAdd(mB, impulse);
        aB += iB * Vec2.CrossVV(rB, impulse);
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
      return inv_dt * (this.lowerImpulse - this.upperImpulse);
    }

    public getJointAngle(): number {
      // Body* bA = this.bodyA;
      // Body* bB = this.bodyB;
      // return bB.this.sweep.a - bA.this.sweep.a - this.referenceAngle;
      return this.bodyB.sweep.a - this.bodyA.sweep.a - this.referenceAngle;
    }

    public getJointSpeed(): number {
      // Body* bA = this.bodyA;
      // Body* bB = this.bodyB;
      // return bB.this.angularVelocity - bA.this.angularVelocity;
      return this.bodyB.angularVelocity - this.bodyA.angularVelocity;
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

    public getMotorTorque(inv_dt: number): number {
      return inv_dt * this.motorImpulse;
    }

    public getMotorSpeed(): number {
      return this.motorSpeed;
    }

    public setMaxMotorTorque(torque: number): void {
      if (torque !== this.maxMotorTorque) {
        this.bodyA.setAwake(true);
        this.bodyB.setAwake(true);
        this.maxMotorTorque = torque;
      }
    }

    public getMaxMotorTorque(): number { return this.maxMotorTorque; }

    public isLimitEnabled(): boolean {
      return this.enableLimit;
    }

    public setEnableLimit(flag: boolean): void {
      if (flag !== this.enableLimit) {
        this.bodyA.setAwake(true);
        this.bodyB.setAwake(true);
        this.enableLimit = flag;
        this.lowerImpulse = 0.0;
        this.upperImpulse = 0.0;
      }
    }

    public getLowerLimit(): number {
      return this.lowerAngle;
    }

    public getUpperLimit(): number {
      return this.upperAngle;
    }

    public setLimits(lower: number, upper: number): void {

      if (lower !== this.lowerAngle || upper !== this.upperAngle) {
        this.bodyA.setAwake(true);
        this.bodyB.setAwake(true);
        this.lowerImpulse = 0.0;
        this.upperImpulse = 0.0;
        this.lowerAngle = lower;
        this.upperAngle = upper;
      }
    }

    public setMotorSpeed(speed: number): void {
      if (speed !== this.motorSpeed) {
        this.bodyA.setAwake(true);
        this.bodyB.setAwake(true);
        this.motorSpeed = speed;
      }
    }

    public dump(log: (format: string, ...args: any[]) => void) {
      const indexA = this.bodyA.islandIndex;
      const indexB = this.bodyB.islandIndex;

      log("  const jd: RevoluteJointDef = new RevoluteJointDef();\n");
      log("  jd.bodyA = bodies[%d];\n", indexA);
      log("  jd.bodyB = bodies[%d];\n", indexB);
      log("  jd.collideConnected = %s;\n", (this.collideConnected) ? ("true") : ("false"));
      log("  jd.localAnchorA.Set(%.15f, %.15f);\n", this.localAnchorA.x, this.localAnchorA.y);
      log("  jd.localAnchorB.Set(%.15f, %.15f);\n", this.localAnchorB.x, this.localAnchorB.y);
      log("  jd.referenceAngle = %.15f;\n", this.referenceAngle);
      log("  jd.enableLimit = %s;\n", (this.enableLimit) ? ("true") : ("false"));
      log("  jd.lowerAngle = %.15f;\n", this.lowerAngle);
      log("  jd.upperAngle = %.15f;\n", this.upperAngle);
      log("  jd.enableMotor = %s;\n", (this.enableMotor) ? ("true") : ("false"));
      log("  jd.motorSpeed = %.15f;\n", this.motorSpeed);
      log("  jd.maxMotorTorque = %.15f;\n", this.maxMotorTorque);
      log("  joints[%d] = this.world.CreateJoint(jd);\n", this.index);
    }

    private static draw_s_pA = new Vec2();
    private static draw_s_pB = new Vec2();
    private static draw_s_c1 = new Color(0.7, 0.7, 0.7);
    private static draw_s_c2 = new Color(0.3, 0.9, 0.3);
    private static draw_s_c3 = new Color(0.9, 0.3, 0.3);
    private static draw_s_c4 = new Color(0.3, 0.3, 0.9);
    private static draw_s_c5 = new Color(0.4, 0.4, 0.4);
    private static draw_s_color_ = new Color(0.5, 0.8, 0.8);
    private static draw_s_r = new Vec2();
    private static draw_s_rlo = new Vec2();
    private static draw_s_rhi = new Vec2();
    public draw(draw: Draw): void {
      const xfA: Transform = this.bodyA.getTransform();
      const xfB: Transform = this.bodyB.getTransform();
      const pA = Transform.mulXV(xfA, this.localAnchorA, RevoluteJoint.draw_s_pA);
      const pB = Transform.mulXV(xfB, this.localAnchorB, RevoluteJoint.draw_s_pB);

      const c1 = RevoluteJoint.draw_s_c1; // Color c1(0.7f, 0.7f, 0.7f);
      const c2 = RevoluteJoint.draw_s_c2; // Color c2(0.3f, 0.9f, 0.3f);
      const c3 = RevoluteJoint.draw_s_c3; // Color c3(0.9f, 0.3f, 0.3f);
      const c4 = RevoluteJoint.draw_s_c4; // Color c4(0.3f, 0.3f, 0.9f);
      const c5 = RevoluteJoint.draw_s_c5; // Color c5(0.4f, 0.4f, 0.4f);

      draw.drawPoint(pA, 5.0, c4);
      draw.drawPoint(pB, 5.0, c5);

      const aA: number = this.bodyA.getAngle();
      const aB: number = this.bodyB.getAngle();
      const angle: number = aB - aA - this.referenceAngle;

      const L: number = 0.5;

      // Vec2 r = L * Vec2(Math.cos(angle), Math.sin(angle));
      const r = RevoluteJoint.draw_s_r.set(L * Math.cos(angle), L * Math.sin(angle));
      // draw.DrawSegment(pB, pB + r, c1);
      draw.drawSegment(pB, Vec2.AddVV(pB, r, Vec2.s_t0), c1);
      draw.drawCircle(pB, L, c1);

      if (this.enableLimit) {
        // Vec2 rlo = L * Vec2(Math.cos(lowerAngle), Math.sin(lowerAngle));
        const rlo = RevoluteJoint.draw_s_rlo.set(L * Math.cos(this.lowerAngle), L * Math.sin(this.lowerAngle));
        // Vec2 rhi = L * Vec2(Math.cos(upperAngle), Math.sin(upperAngle));
        const rhi = RevoluteJoint.draw_s_rhi.set(L * Math.cos(this.upperAngle), L * Math.sin(this.upperAngle));

        // draw.DrawSegment(pB, pB + rlo, c2);
        draw.drawSegment(pB, Vec2.AddVV(pB, rlo, Vec2.s_t0), c2);
        // draw.DrawSegment(pB, pB + rhi, c3);
        draw.drawSegment(pB, Vec2.AddVV(pB, rhi, Vec2.s_t0), c3);
      }

      const color = RevoluteJoint.draw_s_color_; // Color color(0.5f, 0.8f, 0.8f);
      draw.drawSegment(xfA.p, pA, color);
      draw.drawSegment(pA, pB, color);
      draw.drawSegment(xfB.p, pB, color);
    }
  }

}
