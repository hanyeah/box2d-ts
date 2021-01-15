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
  export interface IMouseJointDef extends IJointDef {
    target?: XY;

    maxForce?: number;

    stiffness?: number;

    damping?: number;
  }

/// Mouse joint definition. This requires a world target point,
/// tuning parameters, and the time step.
  export class MouseJointDef extends JointDef implements IMouseJointDef {
    public readonly target: Vec2 = new Vec2();

    public maxForce: number = 0;

    public stiffness: number = 5;

    public damping: number = 0.7;

    constructor() {
      super(JointType.e_mouseJoint);
    }
  }

  export class MouseJoint extends Joint {
    public readonly localAnchorB: Vec2 = new Vec2();
    public readonly targetA: Vec2 = new Vec2();
    public stiffness: number = 0;
    public damping: number = 0;
    public beta: number = 0;

    // Solver shared
    public readonly impulse: Vec2 = new Vec2();
    public maxForce: number = 0;
    public gamma: number = 0;

    // Solver temp
    public indexA: number = 0;
    public indexB: number = 0;
    public readonly rB: Vec2 = new Vec2();
    public readonly localCenterB: Vec2 = new Vec2();
    public invMassB: number = 0;
    public invIB: number = 0;
    public readonly mass: Mat22 = new Mat22();
    public readonly C: Vec2 = new Vec2();
    public readonly qB: Rot = new Rot();
    public readonly lalcB: Vec2 = new Vec2();
    public readonly K: Mat22 = new Mat22();

    constructor(def: IMouseJointDef) {
      super(def);

      this.targetA.Copy(Maybe(def.target, Vec2.ZERO));
      // DEBUG: Assert(this.targetA.IsValid());
      Transform.MulTXV(this.bodyB.GetTransform(), this.targetA, this.localAnchorB);

      this.maxForce = Maybe(def.maxForce, 0);
      // DEBUG: Assert(IsValid(this.maxForce) && this.maxForce >= 0);
      this.impulse.SetZero();

      this.stiffness = Maybe(def.stiffness, 0);
      // DEBUG: Assert(IsValid(this.stiffness) && this.stiffness >= 0);
      this.damping = Maybe(def.damping, 0);
      // DEBUG: Assert(IsValid(this.damping) && this.damping >= 0);

      this.beta = 0;
      this.gamma = 0;
    }

    public SetTarget(target: Vec2): void {
      if (!this.bodyB.IsAwake()) {
        this.bodyB.SetAwake(true);
      }
      this.targetA.Copy(target);
    }

    public GetTarget() {
      return this.targetA;
    }

    public SetMaxForce(maxForce: number): void {
      this.maxForce = maxForce;
    }

    public GetMaxForce() {
      return this.maxForce;
    }

    public SetStiffness(stiffness: number): void {
      this.stiffness = stiffness;
    }

    public GetStiffness() {
      return this.stiffness;
    }

    public SetDamping(damping: number) {
      this.damping = damping;
    }

    public GetDamping() {
      return this.damping;
    }

    public InitVelocityConstraints(data: SolverData): void {
      this.indexB = this.bodyB.islandIndex;
      this.localCenterB.Copy(this.bodyB.sweep.localCenter);
      this.invMassB = this.bodyB.invMass;
      this.invIB = this.bodyB.invI;

      const cB: Vec2 = data.positions[this.indexB].c;
      const aB: number = data.positions[this.indexB].a;
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;

      const qB = this.qB.SetAngle(aB);

      const mass: number = this.bodyB.GetMass();

      // Frequency
      const omega: number = 2 * pi * this.stiffness;

      // Damping coefficient
      const d: number = 2 * mass * this.damping * omega;

      // Spring stiffness
      const k: number = mass * (omega * omega);

      // magic formulas
      // gamma has units of inverse mass.
      // beta has units of inverse time.
      const h: number = data.step.dt;
      this.gamma = h * (d + h * k);
      if (this.gamma !== 0) {
        this.gamma = 1 / this.gamma;
      }
      this.beta = h * k * this.gamma;

      // Compute the effective mass matrix.
      Vec2.SubVV(this.localAnchorB, this.localCenterB, this.lalcB);
      Rot.MulRV(qB, this.lalcB, this.rB);

      // K    = [(1/m1 + 1/m2) * eye(2) - skew(r1) * invI1 * skew(r1) - skew(r2) * invI2 * skew(r2)]
      //      = [1/m1+1/m2     0    ] + invI1 * [r1.y*r1.y -r1.x*r1.y] + invI2 * [r1.y*r1.y -r1.x*r1.y]
      //        [    0     1/m1+1/m2]           [-r1.x*r1.y r1.x*r1.x]           [-r1.x*r1.y r1.x*r1.x]
      const K = this.K;
      K.ex.x = this.invMassB + this.invIB * this.rB.y * this.rB.y + this.gamma;
      K.ex.y = -this.invIB * this.rB.x * this.rB.y;
      K.ey.x = K.ex.y;
      K.ey.y = this.invMassB + this.invIB * this.rB.x * this.rB.x + this.gamma;

      K.GetInverse(this.mass);

      // C = cB + rB - targetA;
      this.C.x = cB.x + this.rB.x - this.targetA.x;
      this.C.y = cB.y + this.rB.y - this.targetA.y;
      // C *= beta;
      this.C.SelfMul(this.beta);

      // Cheat with some damping
      wB *= 0.98;

      if (data.step.warmStarting) {
        this.impulse.SelfMul(data.step.dtRatio);
        // vB += invMassB * impulse;
        vB.x += this.invMassB * this.impulse.x;
        vB.y += this.invMassB * this.impulse.y;
        wB += this.invIB * Vec2.CrossVV(this.rB, this.impulse);
      } else {
        this.impulse.SetZero();
      }

      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    private static SolveVelocityConstraints_s_Cdot = new Vec2();
    private static SolveVelocityConstraints_s_impulse = new Vec2();
    private static SolveVelocityConstraints_s_oldImpulse = new Vec2();
    public SolveVelocityConstraints(data: SolverData): void {
      const vB: Vec2 = data.velocities[this.indexB].v;
      let wB: number = data.velocities[this.indexB].w;

      // Cdot = v + cross(w, r)
      // Vec2 Cdot = vB + Cross(wB, rB);
      const Cdot: Vec2 = Vec2.AddVCrossSV(vB, wB, this.rB, MouseJoint.SolveVelocityConstraints_s_Cdot);
      //  Vec2 impulse = Mul(mass, -(Cdot + C + gamma * impulse));
      const impulse: Vec2 = Mat22.MulMV(
        this.mass,
        Vec2.AddVV(
          Cdot,
          Vec2.AddVV(this.C,
            Vec2.MulSV(this.gamma, this.impulse, Vec2.s_t0),
            Vec2.s_t0),
          Vec2.s_t0).SelfNeg(),
        MouseJoint.SolveVelocityConstraints_s_impulse);

      // Vec2 oldImpulse = impulse;
      const oldImpulse = MouseJoint.SolveVelocityConstraints_s_oldImpulse.Copy(this.impulse);
      // impulse += impulse;
      this.impulse.SelfAdd(impulse);
      const maxImpulse: number = data.step.dt * this.maxForce;
      if (this.impulse.LengthSquared() > maxImpulse * maxImpulse) {
        this.impulse.SelfMul(maxImpulse / this.impulse.Length());
      }
      // impulse = impulse - oldImpulse;
      Vec2.SubVV(this.impulse, oldImpulse, impulse);

      // vB += invMassB * impulse;
      vB.SelfMulAdd(this.invMassB, impulse);
      wB += this.invIB * Vec2.CrossVV(this.rB, impulse);

      // data.velocities[this.indexB].v = vB;
      data.velocities[this.indexB].w = wB;
    }

    public SolvePositionConstraints(data: SolverData): boolean {
      return true;
    }

    public GetAnchorA<T extends XY>(out: T): T {
      out.x = this.targetA.x;
      out.y = this.targetA.y;
      return out;
    }

    public GetAnchorB<T extends XY>(out: T): T {
      return this.bodyB.GetWorldPoint(this.localAnchorB, out);
    }

    public GetReactionForce<T extends XY>(inv_dt: number, out: T): T {
      return Vec2.MulSV(inv_dt, this.impulse, out);
    }

    public GetReactionTorque(inv_dt: number): number {
      return 0;
    }

    public Dump(log: (format: string, ...args: any[]) => void) {
      log("Mouse joint dumping is not supported.\n");
    }

    public ShiftOrigin(newOrigin: Vec2) {
      this.targetA.SelfSub(newOrigin);
    }
  }

}
