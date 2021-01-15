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
  export interface IDistanceJointDef extends IJointDef {
    localAnchorA?: XY;
    localAnchorB?: XY;
    length?: number;
    minLength?: number;
    maxLength?: number;
    stiffness?: number;
    damping?: number;
  }

/// Distance joint definition. This requires defining an
/// anchor point on both bodies and the non-zero length of the
/// distance joint. The definition uses local anchor points
/// so that the initial configuration can violate the constraint
/// slightly. This helps when saving and loading a game.
/// @warning Do not use a zero or short length.
  export class DistanceJointDef extends JointDef implements IDistanceJointDef {
    public readonly localAnchorA: Vec2 = new Vec2();
    public readonly localAnchorB: Vec2 = new Vec2();
    public length: number = 1;
    public minLength: number = 0;
    public maxLength: number = maxFloat; // FLT_MAX;
    public stiffness: number = 0;
    public damping: number = 0;

    constructor() {
      super(JointType.e_distanceJoint);
    }

    public Initialize(b1: Body, b2: Body, anchor1: XY, anchor2: XY): void {
      this.bodyA = b1;
      this.bodyB = b2;
      this.bodyA.GetLocalPoint(anchor1, this.localAnchorA);
      this.bodyB.GetLocalPoint(anchor2, this.localAnchorB);
      this.length = Max(Vec2.DistanceVV(anchor1, anchor2), linearSlop);
      this.minLength = this.length;
      this.maxLength = this.length;
    }
  }

  export class DistanceJoint extends Joint {
    public stiffness: number = 0;
    public damping: number = 0;
    public bias: number = 0;
    public length: number = 0;
    public minLength: number = 0;
    public maxLength: number = 0;

    // Solver shared
    public readonly localAnchorA: Vec2 = new Vec2();
    public readonly localAnchorB: Vec2 = new Vec2();
    public gamma: number = 0;
    public impulse: number = 0;
    public lowerImpulse: number = 0;
    public upperImpulse: number = 0;

    // Solver temp
    public indexA: number = 0;
    public indexB: number = 0;
    public readonly u: Vec2 = new Vec2();
    public readonly rA: Vec2 = new Vec2();
    public readonly rB: Vec2 = new Vec2();
    public readonly localCenterA: Vec2 = new Vec2();
    public readonly localCenterB: Vec2 = new Vec2();
    public currentLength: number = 0;
    public invMassA: number = 0;
    public invMassB: number = 0;
    public invIA: number = 0;
    public invIB: number = 0;
    public softMass: number = 0;
    public mass: number = 0;

    public readonly qA: Rot = new Rot();
    public readonly qB: Rot = new Rot();
    public readonly lalcA: Vec2 = new Vec2();
    public readonly lalcB: Vec2 = new Vec2();

    constructor(def: IDistanceJointDef) {
      super(def);

      this.localAnchorA.Copy(Maybe(def.localAnchorA, Vec2.ZERO));
      this.localAnchorB.Copy(Maybe(def.localAnchorB, Vec2.ZERO));
      this.length = Max(Maybe(def.length, this.GetCurrentLength()), linearSlop);
      this.minLength = Max(Maybe(def.minLength, this.length), linearSlop);
      this.maxLength = Max(Maybe(def.maxLength, this.length), this.minLength);
      this.stiffness = Maybe(def.stiffness, 0);
      this.damping = Maybe(def.damping, 0);
    }

    public GetAnchorA<T extends XY>(out: T): T {
      return this.bodyA.GetWorldPoint(this.localAnchorA, out);
    }

    public GetAnchorB<T extends XY>(out: T): T {
      return this.bodyB.GetWorldPoint(this.localAnchorB, out);
    }

    public GetReactionForce<T extends XY>(inv_dt: number, out: T): T {
      // Vec2 F = inv_dt * (impulse + lowerImpulse - upperImpulse) * u;
      out.x = inv_dt * (this.impulse + this.lowerImpulse - this.upperImpulse) * this.u.x;
      out.y = inv_dt * (this.impulse + this.lowerImpulse - this.upperImpulse) * this.u.y;
      return out;
    }

    public GetReactionTorque(inv_dt: number): number {
      return 0;
    }

    public GetLocalAnchorA(): Vec2 { return this.localAnchorA; }

    public GetLocalAnchorB(): Vec2 { return this.localAnchorB; }

    public SetLength(length: number): number {
      this.impulse = 0;
      this.length = Max(linearSlop, length);
      return this.length;
    }

    public GetLength(): number {
      return this.length;
    }

    public SetMinLength(minLength: number): number {
      this.lowerImpulse = 0;
      this.minLength = Clamp(minLength, linearSlop, this.maxLength);
      return this.minLength;
    }

    public SetMaxLength(maxLength: number): number {
      this.upperImpulse = 0;
      this.maxLength = Max(maxLength, this.minLength);
      return this.maxLength;
    }

    public GetCurrentLength(): number {
      const pA: Vec2 = this.bodyA.GetWorldPoint(this.localAnchorA, new Vec2());
      const pB: Vec2 = this.bodyB.GetWorldPoint(this.localAnchorB, new Vec2());
      return Vec2.DistanceVV(pA, pB);
    }

    public SetStiffness(stiffness: number): void {
      this.stiffness = stiffness;
    }

    public GetStiffness() {
      return this.stiffness;
    }

    public SetDamping(damping: number): void {
      this.damping = damping;
    }

    public GetDamping() {
      return this.damping;
    }

    public Dump(log: (format: string, ...args: any[]) => void) {
      const indexA: number = this.bodyA.islandIndex;
      const indexB: number = this.bodyB.islandIndex;

      log("  const jd: DistanceJointDef = new DistanceJointDef();\n");
      log("  jd.bodyA = bodies[%d];\n", indexA);
      log("  jd.bodyB = bodies[%d];\n", indexB);
      log("  jd.collideConnected = %s;\n", (this.collideConnected) ? ("true") : ("false"));
      log("  jd.localAnchorA.Set(%.15f, %.15f);\n", this.localAnchorA.x, this.localAnchorA.y);
      log("  jd.localAnchorB.Set(%.15f, %.15f);\n", this.localAnchorB.x, this.localAnchorB.y);
      log("  jd.length = %.15f;\n", this.length);
      log("  jd.minLength = %.15f;\n", this.minLength);
      log("  jd.maxLength = %.15f;\n", this.maxLength);
      log("  jd.stiffness = %.15f;\n", this.stiffness);
      log("  jd.damping = %.15f;\n", this.damping);
      log("  joints[%d] = this.world.CreateJoint(jd);\n", this.index);
    }

    private static InitVelocityConstraints_s_P = new Vec2();
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

      // const qA: Rot = new Rot(aA), qB: Rot = new Rot(aB);
      const qA: Rot = this.qA.SetAngle(aA), qB: Rot = this.qB.SetAngle(aB);

      // rA = Mul(qA, localAnchorA - localCenterA);
      Vec2.SubVV(this.localAnchorA, this.localCenterA, this.lalcA);
      Rot.MulRV(qA, this.lalcA, this.rA);
      // rB = Mul(qB, localAnchorB - localCenterB);
      Vec2.SubVV(this.localAnchorB, this.localCenterB, this.lalcB);
      Rot.MulRV(qB, this.lalcB, this.rB);
      // u = cB + rB - cA - rA;
      this.u.x = cB.x + this.rB.x - cA.x - this.rA.x;
      this.u.y = cB.y + this.rB.y - cA.y - this.rA.y;

      // Handle singularity.
      this.currentLength = this.u.Length();
      if (this.currentLength > linearSlop) {
        this.u.SelfMul(1 / this.currentLength);
      } else {
        this.u.SetZero();
        this.mass = 0;
        this.impulse = 0;
        this.lowerImpulse = 0;
        this.upperImpulse = 0;
      }

      // float32 crAu = Cross(rA, u);
      const crAu: number = Vec2.CrossVV(this.rA, this.u);
      // float32 crBu = Cross(rB, u);
      const crBu: number = Vec2.CrossVV(this.rB, this.u);
      // float32 invMass = invMassA + invIA * crAu * crAu + invMassB + invIB * crBu * crBu;
      let invMass: number = this.invMassA + this.invIA * crAu * crAu + this.invMassB + this.invIB * crBu * crBu;
      this.mass = invMass !== 0 ? 1 / invMass : 0;

      if (this.stiffness > 0 && this.minLength < this.maxLength) {
        // soft
        const C: number = this.currentLength - this.length;

        const d: number = this.damping;
        const k: number = this.stiffness;

        // magic formulas
        const h: number = data.step.dt;

        // gamma = 1 / (h * (d + h * k))
        // the extra factor of h in the denominator is since the lambda is an impulse, not a force
        this.gamma = h * (d + h * k);
        this.gamma = this.gamma !== 0 ? 1 / this.gamma : 0;
        this.bias = C * h * k * this.gamma;

        invMass += this.gamma;
        this.softMass = invMass !== 0 ? 1 / invMass : 0;
      }
      else {
        // rigid
        this.gamma = 0;
        this.bias = 0;
        this.softMass = this.mass;
      }

      if (data.step.warmStarting) {
        // Scale the impulse to support a variable time step.
        this.impulse *= data.step.dtRatio;
        this.lowerImpulse *= data.step.dtRatio;
        this.upperImpulse *= data.step.dtRatio;

        const P: Vec2 = Vec2.MulSV(this.impulse + this.lowerImpulse - this.upperImpulse, this.u, DistanceJoint.InitVelocityConstraints_s_P);
        vA.SelfMulSub(this.invMassA, P);
        wA -= this.invIA * Vec2.CrossVV(this.rA, P);
        vB.SelfMulAdd(this.invMassB, P);
        wB += this.invIB * Vec2.CrossVV(this.rB, P);
      }
      else {
        this.impulse = 0;
      }

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    private static SolveVelocityConstraints_s_vpA = new Vec2();
    private static SolveVelocityConstraints_s_vpB = new Vec2();
    private static SolveVelocityConstraints_s_P = new Vec2();
    public SolveVelocityConstraints(data: SolverData): void {
      const vA: Vec2 = data.velocities[this.indexA].v;
      let wA: number = data.velocities[this.indexA].w;
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;

      if (this.minLength < this.maxLength) {
        if (this.stiffness > 0) {
          // Cdot = dot(u, v + cross(w, r))
          const vpA: Vec2 = Vec2.AddVCrossSV(vA, wA, this.rA, DistanceJoint.SolveVelocityConstraints_s_vpA);
          const vpB: Vec2 = Vec2.AddVCrossSV(vB, wB, this.rB, DistanceJoint.SolveVelocityConstraints_s_vpB);
          const Cdot: number = Vec2.DotVV(this.u, Vec2.SubVV(vpB, vpA, Vec2.s_t0));

          const impulse: number = -this.softMass * (Cdot + this.bias + this.gamma * this.impulse);
          this.impulse += impulse;

          const P: Vec2 = Vec2.MulSV(impulse, this.u, DistanceJoint.SolveVelocityConstraints_s_P);
          vA.SelfMulSub(this.invMassA, P);
          wA -= this.invIA * Vec2.CrossVV(this.rA, P);
          vB.SelfMulAdd(this.invMassB, P);
          wB += this.invIB * Vec2.CrossVV(this.rB, P);
        }

        // lower
        {
          const C: number = this.currentLength - this.minLength;
          const bias: number = Max(0, C) * data.step.inv_dt;

          const vpA: Vec2 = Vec2.AddVCrossSV(vA, wA, this.rA, DistanceJoint.SolveVelocityConstraints_s_vpA);
          const vpB: Vec2 = Vec2.AddVCrossSV(vB, wB, this.rB, DistanceJoint.SolveVelocityConstraints_s_vpB);
          const Cdot: number = Vec2.DotVV(this.u, Vec2.SubVV(vpB, vpA, Vec2.s_t0));

          let impulse: number = -this.mass * (Cdot + bias);
          const oldImpulse: number = this.lowerImpulse;
          this.lowerImpulse = Max(0, this.lowerImpulse + impulse);
          impulse = this.lowerImpulse - oldImpulse;
          const P: Vec2 = Vec2.MulSV(impulse, this.u, DistanceJoint.SolveVelocityConstraints_s_P);

          vA.SelfMulSub(this.invMassA, P);
          wA -= this.invIA * Vec2.CrossVV(this.rA, P);
          vB.SelfMulAdd(this.invMassB, P);
          wB += this.invIB * Vec2.CrossVV(this.rB, P);
        }

        // upper
        {
          const C: number = this.maxLength - this.currentLength;
          const bias: number = Max(0, C) * data.step.inv_dt;

          const vpA: Vec2 = Vec2.AddVCrossSV(vA, wA, this.rA, DistanceJoint.SolveVelocityConstraints_s_vpA);
          const vpB: Vec2 = Vec2.AddVCrossSV(vB, wB, this.rB, DistanceJoint.SolveVelocityConstraints_s_vpB);
          const Cdot: number = Vec2.DotVV(this.u, Vec2.SubVV(vpA, vpB, Vec2.s_t0));

          let impulse: number = -this.mass * (Cdot + bias);
          const oldImpulse: number = this.upperImpulse;
          this.upperImpulse = Max(0, this.upperImpulse + impulse);
          impulse = this.upperImpulse - oldImpulse;
          const P: Vec2 = Vec2.MulSV(-impulse, this.u, DistanceJoint.SolveVelocityConstraints_s_P);

          vA.SelfMulSub(this.invMassA, P);
          wA -= this.invIA * Vec2.CrossVV(this.rA, P);
          vB.SelfMulAdd(this.invMassB, P);
          wB += this.invIB * Vec2.CrossVV(this.rB, P);
        }
      }
      else {
        // Equal limits

        // Cdot = dot(u, v + cross(w, r))
        const vpA: Vec2 = Vec2.AddVCrossSV(vA, wA, this.rA, DistanceJoint.SolveVelocityConstraints_s_vpA);
        const vpB: Vec2 = Vec2.AddVCrossSV(vB, wB, this.rB, DistanceJoint.SolveVelocityConstraints_s_vpB);
        const Cdot: number = Vec2.DotVV(this.u, Vec2.SubVV(vpB, vpA, Vec2.s_t0));

        const impulse: number = -this.mass * Cdot;
        this.impulse += impulse;

        const P: Vec2 = Vec2.MulSV(impulse, this.u, DistanceJoint.SolveVelocityConstraints_s_P);
        vA.SelfMulSub(this.invMassA, P);
        wA -= this.invIA * Vec2.CrossVV(this.rA, P);
        vB.SelfMulAdd(this.invMassB, P);
        wB += this.invIB * Vec2.CrossVV(this.rB, P);
      }

      // data.velocities[this.indexA].v = vA;
      data.velocities[this.indexA].w = wA;
      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    private static SolvePositionConstraints_s_P = new Vec2();
    public SolvePositionConstraints(data: SolverData): boolean {
      const cA: Vec2 = data.positions[this.indexA].c;
      let aA: number = data.positions[this.indexA].a;
      const cB: Vec2 = data.positions[this.indexB].c;
      let aB: number = data.positions[this.indexB].a;

      // const qA: Rot = new Rot(aA), qB: Rot = new Rot(aB);
      const qA: Rot = this.qA.SetAngle(aA), qB: Rot = this.qB.SetAngle(aB);

      // Vec2 rA = Mul(qA, localAnchorA - localCenterA);
      const rA: Vec2 = Rot.MulRV(qA, this.lalcA, this.rA); // use rA
      // Vec2 rB = Mul(qB, localAnchorB - localCenterB);
      const rB: Vec2 = Rot.MulRV(qB, this.lalcB, this.rB); // use rB
      // Vec2 u = cB + rB - cA - rA;
      const u: Vec2 = this.u; // use u
      u.x = cB.x + rB.x - cA.x - rA.x;
      u.y = cB.y + rB.y - cA.y - rA.y;

      const length: number = this.u.Normalize();
      let C: number;
      if (this.minLength == this.maxLength)
      {
        C = length - this.minLength;
      }
      else if (length < this.minLength)
      {
        C = length - this.minLength;
      }
      else if (this.maxLength < length)
      {
        C = length - this.maxLength;
      }
      else
      {
        return true;
      }

      const impulse: number = -this.mass * C;
      const P: Vec2 = Vec2.MulSV(impulse, u, DistanceJoint.SolvePositionConstraints_s_P);

      cA.SelfMulSub(this.invMassA, P);
      aA -= this.invIA * Vec2.CrossVV(rA, P);
      cB.SelfMulAdd(this.invMassB, P);
      aB += this.invIB * Vec2.CrossVV(rB, P);

      // data.positions[this.indexA].c = cA;
      data.positions[this.indexA].a = aA;
      // data.positions[this.indexB].c = cB;
      data.positions[this.indexB].a = aB;

      return Abs(C) < linearSlop;
    }

    private static Draw_s_pA = new Vec2();
    private static Draw_s_pB = new Vec2();
    private static Draw_s_axis = new Vec2();
    private static Draw_s_c1 = new Color(0.7, 0.7, 0.7);
    private static Draw_s_c2 = new Color(0.3, 0.9, 0.3);
    private static Draw_s_c3 = new Color(0.9, 0.3, 0.3);
    private static Draw_s_c4 = new Color(0.4, 0.4, 0.4);
    private static Draw_s_pRest = new Vec2();
    private static Draw_s_pMin = new Vec2();
    private static Draw_s_pMax = new Vec2();
    public Draw(draw: Draw): void {
      const xfA: Transform = this.bodyA.GetTransform();
      const xfB: Transform = this.bodyB.GetTransform();
      const pA = Transform.MulXV(xfA, this.localAnchorA, DistanceJoint.Draw_s_pA);
      const pB = Transform.MulXV(xfB, this.localAnchorB, DistanceJoint.Draw_s_pB);

      const axis: Vec2 = Vec2.SubVV(pB, pA, DistanceJoint.Draw_s_axis);
      axis.Normalize();

      const c1 = DistanceJoint.Draw_s_c1; // Color c1(0.7f, 0.7f, 0.7f);
      const c2 = DistanceJoint.Draw_s_c2; // Color c2(0.3f, 0.9f, 0.3f);
      const c3 = DistanceJoint.Draw_s_c3; // Color c3(0.9f, 0.3f, 0.3f);
      const c4 = DistanceJoint.Draw_s_c4; // Color c4(0.4f, 0.4f, 0.4f);

      draw.DrawSegment(pA, pB, c4);

      // Vec2 pRest = pA + this.length * axis;
      const pRest: Vec2 = Vec2.AddVMulSV(pA, this.length, axis, DistanceJoint.Draw_s_pRest);
      draw.DrawPoint(pRest, 8.0, c1);

      if (this.minLength != this.maxLength) {
        if (this.minLength > linearSlop) {
          // Vec2 pMin = pA + this.minLength * axis;
          const pMin: Vec2 = Vec2.AddVMulSV(pA, this.minLength, axis, DistanceJoint.Draw_s_pMin);
          draw.DrawPoint(pMin, 4.0, c2);
        }

        if (this.maxLength < maxFloat) {
          // Vec2 pMax = pA + this.maxLength * axis;
          const pMax: Vec2 = Vec2.AddVMulSV(pA, this.maxLength, axis, DistanceJoint.Draw_s_pMax);
          draw.DrawPoint(pMax, 4.0, c3);
        }
      }
    }
  }

}
