// MIT License

// Copyright (c) 2019 Erin Catto

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
namespace b2 {
  export enum StretchingModel {
    pbdStretchingModel,
    xpbdStretchingModel,
  }

  export enum BendingModel {
    springAngleBendingModel = 0,
    pbdAngleBendingModel,
    xpbdAngleBendingModel,
    pbdDistanceBendingModel,
    pbdHeightBendingModel,
    pbdTriangleBendingModel,
  }

///
  export class RopeTuning {
    public stretchingModel: StretchingModel = StretchingModel.pbdStretchingModel;
    public bendingModel: BendingModel = BendingModel.pbdAngleBendingModel;
    public damping: number = 0.0;
    public stretchStiffness: number = 1.0;
    public stretchHertz: number = 0.0;
    public stretchDamping: number = 0.0;
    public bendStiffness: number = 0.5;
    public bendHertz: number = 1.0;
    public bendDamping: number = 0.0;
    public isometric: boolean = false;
    public fixedEffectiveMass: boolean = false;
    public warmStart: boolean = false;

    public Copy(other: RopeTuning): this {
      this.stretchingModel = other.stretchingModel;
      this.bendingModel = other.bendingModel;
      this.damping = other.damping;
      this.stretchStiffness = other.stretchStiffness;
      this.stretchHertz = other.stretchHertz;
      this.stretchDamping = other.stretchDamping;
      this.bendStiffness = other.bendStiffness;
      this.bendHertz = other.bendHertz;
      this.bendDamping = other.bendDamping;
      this.isometric = other.isometric;
      this.fixedEffectiveMass = other.fixedEffectiveMass;
      this.warmStart = other.warmStart;
      return this;
    }
  }

///
  export class RopeDef {
    public readonly position: Vec2 = new Vec2();
    // Vec2* vertices;
    public readonly vertices: Vec2[] = [];
    // int32 count;
    public count: number = 0;
    // float* masses;
    public readonly masses: number[] = [];
    // Vec2 gravity;
    public readonly gravity: Vec2 = new Vec2();
    // RopeTuning tuning;
    public readonly tuning: RopeTuning = new RopeTuning();
  }

  class RopeStretch {
    public i1: number = 0;
    public i2: number = 0;
    public invMass1: number = 0.0;
    public invMass2: number = 0.0;
    public L: number = 0.0;
    public lambda: number = 0.0;
    public spring: number = 0.0;
    public damper: number = 0.0;
  }

  class RopeBend {
    public i1: number = 0;
    public i2: number = 0;
    public i3: number = 0;
    public invMass1: number = 0.0;
    public invMass2: number = 0.0;
    public invMass3: number = 0.0;
    public invEffectiveMass: number = 0.0;
    public lambda: number = 0.0;
    public L1: number = 0.0;
    public L2: number = 0.0;
    public alpha1: number = 0.0;
    public alpha2: number = 0.0;
    public spring: number = 0.0;
    public damper: number = 0.0;
  }

///
  export class Rope {
    private readonly position: Vec2 = new Vec2();

    private count: number = 0;
    private stretchCount: number = 0;
    private bendCount: number = 0;

    // RopeStretch* stretchConstraints;
    private readonly stretchConstraints: RopeStretch[] = [];
    // RopeBend* bendConstraints;
    private readonly bendConstraints: RopeBend[] = [];

    // Vec2* bindPositions;
    private readonly bindPositions: Vec2[] = [];
    // Vec2* ps;
    private readonly ps: Vec2[] = [];
    // Vec2* p0s;
    private readonly p0s: Vec2[] = [];
    // Vec2* vs;
    private readonly vs: Vec2[] = [];

    // float* invMasses;
    private readonly invMasses: number[] = [];
    // Vec2 gravity;
    private readonly gravity: Vec2 = new Vec2();

    private readonly tuning: RopeTuning = new RopeTuning();

    public Create(def: RopeDef): void {
      // Assert(def.count >= 3);
      this.position.Copy(def.position);
      this.count = def.count;
      function make_array<T>(array: T[], count: number, make: (index: number) => T): void {
        for (let index = 0; index < count; ++index) {
          array[index] = make(index);
        }
      }
      // this.bindPositions = (Vec2*)Alloc(this.count * sizeof(Vec2));
      make_array(this.bindPositions, this.count, () => new Vec2());
      // this.ps = (Vec2*)Alloc(this.count * sizeof(Vec2));
      make_array(this.ps, this.count, () => new Vec2());
      // this.p0s = (Vec2*)Alloc(this.count * sizeof(Vec2));
      make_array(this.p0s, this.count, () => new Vec2());
      // this.vs = (Vec2*)Alloc(this.count * sizeof(Vec2));
      make_array(this.vs, this.count, () => new Vec2());
      // this.invMasses = (float*)Alloc(this.count * sizeof(float));
      make_array(this.invMasses, this.count, () => 0.0);

      for (let i = 0; i < this.count; ++i) {
        this.bindPositions[i].Copy(def.vertices[i]);
        // this.ps[i] = def.vertices[i] + this.position;
        this.ps[i].Copy(def.vertices[i]).SelfAdd(this.position);
        // this.p0s[i] = def.vertices[i] + this.position;
        this.p0s[i].Copy(def.vertices[i]).SelfAdd(this.position);
        this.vs[i].SetZero();

        const m: number = def.masses[i];
        if (m > 0.0) {
          this.invMasses[i] = 1.0 / m;
        } else {
          this.invMasses[i] = 0.0;
        }
      }

      this.stretchCount = this.count - 1;
      this.bendCount = this.count - 2;

      // this.stretchConstraints = (RopeStretch*)Alloc(this.stretchCount * sizeof(RopeStretch));
      make_array(this.stretchConstraints, this.stretchCount, () => new RopeStretch());
      // this.bendConstraints = (RopeBend*)Alloc(this.bendCount * sizeof(RopeBend));
      make_array(this.bendConstraints, this.bendCount, () => new RopeBend());

      for (let i = 0; i < this.stretchCount; ++i) {
        const c: RopeStretch = this.stretchConstraints[i];

        const p1: Vec2 = this.ps[i];
        const p2: Vec2 = this.ps[i + 1];

        c.i1 = i;
        c.i2 = i + 1;
        c.L = Vec2.DistanceVV(p1, p2);
        c.invMass1 = this.invMasses[i];
        c.invMass2 = this.invMasses[i + 1];
        c.lambda = 0.0;
        c.damper = 0.0;
        c.spring = 0.0;
      }

      for (let i = 0; i < this.bendCount; ++i) {
        const c: RopeBend = this.bendConstraints[i];

        const p1: Vec2 = this.ps[i];
        const p2: Vec2 = this.ps[i + 1];
        const p3: Vec2 = this.ps[i + 2];

        c.i1 = i;
        c.i2 = i + 1;
        c.i3 = i + 2;
        c.invMass1 = this.invMasses[i];
        c.invMass2 = this.invMasses[i + 1];
        c.invMass3 = this.invMasses[i + 2];
        c.invEffectiveMass = 0.0;
        c.L1 = Vec2.DistanceVV(p1, p2);
        c.L2 = Vec2.DistanceVV(p2, p3);
        c.lambda = 0.0;

        // Pre-compute effective mass (TODO use flattened config)
        const e1: Vec2 = Vec2.SubVV(p2, p1, new Vec2());
        const e2: Vec2 = Vec2.SubVV(p3, p2, new Vec2());
        const L1sqr: number = e1.LengthSquared();
        const L2sqr: number = e2.LengthSquared();

        if (L1sqr * L2sqr === 0.0) {
          continue;
        }

        // Vec2 Jd1 = (-1.0 / L1sqr) * e1.Skew();
        const Jd1: Vec2 = new Vec2().Copy(e1).SelfSkew().SelfMul(-1.0 / L1sqr);
        // Vec2 Jd2 = (1.0 / L2sqr) * e2.Skew();
        const Jd2: Vec2 = new Vec2().Copy(e2).SelfSkew().SelfMul(1.0 / L2sqr);

        // Vec2 J1 = -Jd1;
        const J1 = Jd1.Clone().SelfNeg();
        // Vec2 J2 = Jd1 - Jd2;
        const J2 = Jd1.Clone().SelfSub(Jd2);
        // Vec2 J3 = Jd2;
        const J3 = Jd2.Clone();

        c.invEffectiveMass = c.invMass1 * Vec2.DotVV(J1, J1) + c.invMass2 * Vec2.DotVV(J2, J2) + c.invMass3 * Vec2.DotVV(J3, J3);

        // Vec2 r = p3 - p1;
        const r: Vec2 = Vec2.SubVV(p3, p1, new Vec2());

        const rr: number = r.LengthSquared();
        if (rr === 0.0) {
          continue;
        }

        // a1 = h2 / (h1 + h2)
        // a2 = h1 / (h1 + h2)
        c.alpha1 = Vec2.DotVV(e2, r) / rr;
        c.alpha2 = Vec2.DotVV(e1, r) / rr;
      }

      this.gravity.Copy(def.gravity);

      this.SetTuning(def.tuning);
    }

    public SetTuning(tuning: RopeTuning): void {
      this.tuning.Copy(tuning);

      // Pre-compute spring and damper values based on tuning

      const bendOmega: number = 2.0 * pi * this.tuning.bendHertz;

      for (let i = 0; i < this.bendCount; ++i) {
        const c: RopeBend = this.bendConstraints[i];

        const L1sqr: number = c.L1 * c.L1;
        const L2sqr: number = c.L2 * c.L2;

        if (L1sqr * L2sqr === 0.0) {
          c.spring = 0.0;
          c.damper = 0.0;
          continue;
        }

        // Flatten the triangle formed by the two edges
        const J2: number = 1.0 / c.L1 + 1.0 / c.L2;
        const sum: number = c.invMass1 / L1sqr + c.invMass2 * J2 * J2 + c.invMass3 / L2sqr;
        if (sum === 0.0) {
          c.spring = 0.0;
          c.damper = 0.0;
          continue;
        }

        const mass: number = 1.0 / sum;

        c.spring = mass * bendOmega * bendOmega;
        c.damper = 2.0 * mass * this.tuning.bendDamping * bendOmega;
      }

      const stretchOmega: number = 2.0 * pi * this.tuning.stretchHertz;

      for (let i = 0; i < this.stretchCount; ++i) {
        const c: RopeStretch = this.stretchConstraints[i];

        const sum: number = c.invMass1 + c.invMass2;
        if (sum === 0.0) {
          continue;
        }

        const mass: number = 1.0 / sum;

        c.spring = mass * stretchOmega * stretchOmega;
        c.damper = 2.0 * mass * this.tuning.stretchDamping * stretchOmega;
      }
    }

    public Step(dt: number, iterations: number, position: Vec2): void {
      if (dt === 0.0) {
        return;
      }

      const inv_dt: number = 1.0 / dt;
      const d: number = Math.exp(- dt * this.tuning.damping);

      // Apply gravity and damping
      for (let i = 0; i < this.count; ++i) {
        if (this.invMasses[i] > 0.0) {
          // this.vs[i] *= d;
          this.vs[i].x *= d;
          this.vs[i].y *= d;
          // this.vs[i] += dt * this.gravity;
          this.vs[i].x += dt * this.gravity.x;
          this.vs[i].y += dt * this.gravity.y;
        }
        else {
          // this.vs[i] = inv_dt * (this.bindPositions[i] + position - this.p0s[i]);
          this.vs[i].x = inv_dt * (this.bindPositions[i].x + position.x - this.p0s[i].x);
          this.vs[i].y = inv_dt * (this.bindPositions[i].y + position.y - this.p0s[i].y);
        }
      }

      // Apply bending spring
      if (this.tuning.bendingModel === BendingModel.springAngleBendingModel) {
        this.ApplyBendForces(dt);
      }

      for (let i = 0; i < this.bendCount; ++i) {
        this.bendConstraints[i].lambda = 0.0;
      }

      for (let i = 0; i < this.stretchCount; ++i) {
        this.stretchConstraints[i].lambda = 0.0;
      }

      // Update position
      for (let i = 0; i < this.count; ++i) {
        // this.ps[i] += dt * this.vs[i];
        this.ps[i].x += dt * this.vs[i].x;
        this.ps[i].y += dt * this.vs[i].y;
      }

      // Solve constraints
      for (let i = 0; i < iterations; ++i) {
        if (this.tuning.bendingModel === BendingModel.pbdAngleBendingModel) {
          this.SolveBend_PBD_Angle();
        }
        else if (this.tuning.bendingModel === BendingModel.xpbdAngleBendingModel) {
          this.SolveBend_XPBD_Angle(dt);
        }
        else if (this.tuning.bendingModel === BendingModel.pbdDistanceBendingModel) {
          this.SolveBend_PBD_Distance();
        }
        else if (this.tuning.bendingModel === BendingModel.pbdHeightBendingModel) {
          this.SolveBend_PBD_Height();
        }
        else if (this.tuning.bendingModel === BendingModel.pbdTriangleBendingModel) {
          this.SolveBend_PBD_Triangle();
        }

        if (this.tuning.stretchingModel === StretchingModel.pbdStretchingModel) {
          this.SolveStretch_PBD();
        }
        else if (this.tuning.stretchingModel === StretchingModel.xpbdStretchingModel) {
          this.SolveStretch_XPBD(dt);
        }
      }

      // Constrain velocity
      for (let i = 0; i < this.count; ++i) {
        // this.vs[i] = inv_dt * (this.ps[i] - this.p0s[i]);
        this.vs[i].x = inv_dt * (this.ps[i].x - this.p0s[i].x);
        this.vs[i].y = inv_dt * (this.ps[i].y - this.p0s[i].y);
        this.p0s[i].Copy(this.ps[i]);
      }
    }

    public Reset(position: Vec2): void {
      this.position.Copy(position);

      for (let i = 0; i < this.count; ++i) {
        // this.ps[i] = this.bindPositions[i] + this.position;
        this.ps[i].x = this.bindPositions[i].x + this.position.x;
        this.ps[i].y = this.bindPositions[i].y + this.position.y;
        // this.p0s[i] = this.bindPositions[i] + this.position;
        this.p0s[i].x = this.bindPositions[i].x + this.position.x;
        this.p0s[i].y = this.bindPositions[i].y + this.position.y;
        this.vs[i].SetZero();
      }

      for (let i = 0; i < this.bendCount; ++i) {
        this.bendConstraints[i].lambda = 0.0;
      }

      for (let i = 0; i < this.stretchCount; ++i) {
        this.stretchConstraints[i].lambda = 0.0;
      }
    }

    public Draw(draw: Draw): void {
      const c: Color = new Color(0.4, 0.5, 0.7);
      const pg: Color = new Color(0.1, 0.8, 0.1);
      const pd: Color = new Color(0.7, 0.2, 0.4);

      for (let i = 0; i < this.count - 1; ++i) {
        draw.DrawSegment(this.ps[i], this.ps[i + 1], c);

        const pc: Color = this.invMasses[i] > 0.0 ? pd : pg;
        draw.DrawPoint(this.ps[i], 5.0, pc);
      }

      const pc: Color = this.invMasses[this.count - 1] > 0.0 ? pd : pg;
      draw.DrawPoint(this.ps[this.count - 1], 5.0, pc);
    }

    private SolveStretch_PBD(): void {
      const stiffness: number = this.tuning.stretchStiffness;

      for (let i = 0; i < this.stretchCount; ++i) {
        const c: RopeStretch = this.stretchConstraints[i];

        const p1: Vec2 = this.ps[c.i1].Clone();
        const p2: Vec2 = this.ps[c.i2].Clone();

        // Vec2 d = p2 - p1;
        const d: Vec2 = p2.Clone().SelfSub(p1);
        const L: number = d.Normalize();

        const sum: number = c.invMass1 + c.invMass2;
        if (sum === 0.0) {
          continue;
        }

        const s1: number = c.invMass1 / sum;
        const s2: number = c.invMass2 / sum;

        // p1 -= stiffness * s1 * (c.L - L) * d;
        p1.x -= stiffness * s1 * (c.L - L) * d.x;
        p1.y -= stiffness * s1 * (c.L - L) * d.y;
        // p2 += stiffness * s2 * (c.L - L) * d;
        p2.x += stiffness * s2 * (c.L - L) * d.x;
        p2.y += stiffness * s2 * (c.L - L) * d.y;

        this.ps[c.i1].Copy(p1);
        this.ps[c.i2].Copy(p2);
      }
    }

    private SolveStretch_XPBD(dt: number): void {
      // 	Assert(dt > 0.0);

      for (let i = 0; i < this.stretchCount; ++i) {
        const c: RopeStretch = this.stretchConstraints[i];

        const p1: Vec2 = this.ps[c.i1].Clone();
        const p2: Vec2 = this.ps[c.i2].Clone();

        const dp1: Vec2 = p1.Clone().SelfSub(this.p0s[c.i1]);
        const dp2: Vec2 = p2.Clone().SelfSub(this.p0s[c.i2]);

        // Vec2 u = p2 - p1;
        const u: Vec2 = p2.Clone().SelfSub(p1);
        const L: number = u.Normalize();

        // Vec2 J1 = -u;
        const J1: Vec2 = u.Clone().SelfNeg();
        // Vec2 J2 = u;
        const J2: Vec2 = u;

        const sum: number = c.invMass1 + c.invMass2;
        if (sum === 0.0) {
          continue;
        }

        const alpha: number = 1.0 / (c.spring * dt * dt);	// 1 / kg
        const beta: number = dt * dt * c.damper;				// kg * s
        const sigma: number = alpha * beta / dt;				// non-dimensional
        const C: number = L - c.L;

        // This is using the initial velocities
        const Cdot: number = Vec2.DotVV(J1, dp1) + Vec2.DotVV(J2, dp2);

        const B: number = C + alpha * c.lambda + sigma * Cdot;
        const sum2: number = (1.0 + sigma) * sum + alpha;

        const impulse: number = -B / sum2;

        // p1 += (c.invMass1 * impulse) * J1;
        p1.x += (c.invMass1 * impulse) * J1.x;
        p1.y += (c.invMass1 * impulse) * J1.y;
        // p2 += (c.invMass2 * impulse) * J2;
        p2.x += (c.invMass2 * impulse) * J2.x;
        p2.y += (c.invMass2 * impulse) * J2.y;

        this.ps[c.i1].Copy(p1);
        this.ps[c.i2].Copy(p2);
        c.lambda += impulse;
      }
    }

    private SolveBend_PBD_Angle(): void {
      const stiffness: number = this.tuning.bendStiffness;

      for (let i = 0; i < this.bendCount; ++i) {
        const c: RopeBend = this.bendConstraints[i];

        const p1: Vec2 = this.ps[c.i1];
        const p2: Vec2 = this.ps[c.i2];
        const p3: Vec2 = this.ps[c.i3];

        // Vec2 d1 = p2 - p1;
        const d1 = p2.Clone().SelfSub(p1);
        // Vec2 d2 = p3 - p2;
        const d2 = p3.Clone().SelfSub(p2);
        const a: number = Vec2.CrossVV(d1, d2);
        const b: number = Vec2.DotVV(d1, d2);

        const angle: number = Atan2(a, b);

        let L1sqr: number = 0.0, L2sqr: number = 0.0;

        if (this.tuning.isometric) {
          L1sqr = c.L1 * c.L1;
          L2sqr = c.L2 * c.L2;
        }
        else {
          L1sqr = d1.LengthSquared();
          L2sqr = d2.LengthSquared();
        }

        if (L1sqr * L2sqr === 0.0) {
          continue;
        }

        // Vec2 Jd1 = (-1.0 / L1sqr) * d1.Skew();
        const Jd1: Vec2 = new Vec2().Copy(d1).SelfSkew().SelfMul(-1.0 / L1sqr);
        // Vec2 Jd2 = (1.0 / L2sqr) * d2.Skew();
        const Jd2: Vec2 = new Vec2().Copy(d2).SelfSkew().SelfMul(1.0 / L2sqr);

        // Vec2 J1 = -Jd1;
        const J1 = Jd1.Clone().SelfNeg();
        // Vec2 J2 = Jd1 - Jd2;
        const J2 = Jd1.Clone().SelfSub(Jd2);
        // Vec2 J3 = Jd2;
        const J3 = Jd2;

        let sum: number = 0.0;
        if (this.tuning.fixedEffectiveMass) {
          sum = c.invEffectiveMass;
        }
        else {
          sum = c.invMass1 * Vec2.DotVV(J1, J1) + c.invMass2 * Vec2.DotVV(J2, J2) + c.invMass3 * Vec2.DotVV(J3, J3);
        }

        if (sum === 0.0) {
          sum = c.invEffectiveMass;
        }

        const impulse: number = -stiffness * angle / sum;

        // p1 += (c.invMass1 * impulse) * J1;
        p1.x += (c.invMass1 * impulse) * J1.x;
        p1.y += (c.invMass1 * impulse) * J1.y;
        // p2 += (c.invMass2 * impulse) * J2;
        p2.x += (c.invMass2 * impulse) * J2.x;
        p2.y += (c.invMass2 * impulse) * J2.y;
        // p3 += (c.invMass3 * impulse) * J3;
        p3.x += (c.invMass3 * impulse) * J3.x;
        p3.y += (c.invMass3 * impulse) * J3.y;

        this.ps[c.i1].Copy(p1);
        this.ps[c.i2].Copy(p2);
        this.ps[c.i3].Copy(p3);
      }
    }

    private SolveBend_XPBD_Angle(dt: number): void {
      // Assert(dt > 0.0);

      for (let i = 0; i < this.bendCount; ++i) {
        const c: RopeBend = this.bendConstraints[i];

        const p1: Vec2 = this.ps[c.i1];
        const p2: Vec2 = this.ps[c.i2];
        const p3: Vec2 = this.ps[c.i3];

        const dp1: Vec2 = p1.Clone().SelfSub(this.p0s[c.i1]);
        const dp2: Vec2 = p2.Clone().SelfSub(this.p0s[c.i2]);
        const dp3: Vec2 = p3.Clone().SelfSub(this.p0s[c.i3]);

        // Vec2 d1 = p2 - p1;
        const d1 = p2.Clone().SelfSub(p1);
        // Vec2 d2 = p3 - p2;
        const d2 = p3.Clone().SelfSub(p2);

        let L1sqr: number, L2sqr: number;

        if (this.tuning.isometric) {
          L1sqr = c.L1 * c.L1;
          L2sqr = c.L2 * c.L2;
        }
        else {
          L1sqr = d1.LengthSquared();
          L2sqr = d2.LengthSquared();
        }

        if (L1sqr * L2sqr === 0.0) {
          continue;
        }

        const a: number = Vec2.CrossVV(d1, d2);
        const b: number = Vec2.DotVV(d1, d2);

        const angle: number = Atan2(a, b);

        // Vec2 Jd1 = (-1.0 / L1sqr) * d1.Skew();
        // Vec2 Jd2 = (1.0 / L2sqr) * d2.Skew();

        // Vec2 J1 = -Jd1;
        // Vec2 J2 = Jd1 - Jd2;
        // Vec2 J3 = Jd2;

        // Vec2 Jd1 = (-1.0 / L1sqr) * d1.Skew();
        const Jd1: Vec2 = new Vec2().Copy(d1).SelfSkew().SelfMul(-1.0 / L1sqr);
        // Vec2 Jd2 = (1.0 / L2sqr) * d2.Skew();
        const Jd2: Vec2 = new Vec2().Copy(d2).SelfSkew().SelfMul(1.0 / L2sqr);

        // Vec2 J1 = -Jd1;
        const J1 = Jd1.Clone().SelfNeg();
        // Vec2 J2 = Jd1 - Jd2;
        const J2 = Jd1.Clone().SelfSub(Jd2);
        // Vec2 J3 = Jd2;
        const J3 = Jd2;

        let sum: number;
        if (this.tuning.fixedEffectiveMass) {
          sum = c.invEffectiveMass;
        }
        else {
          sum = c.invMass1 * Vec2.DotVV(J1, J1) + c.invMass2 * Vec2.DotVV(J2, J2) + c.invMass3 * Vec2.DotVV(J3, J3);
        }

        if (sum === 0.0) {
          continue;
        }

        const alpha: number = 1.0 / (c.spring * dt * dt);
        const beta: number = dt * dt * c.damper;
        const sigma: number = alpha * beta / dt;
        const C: number = angle;

        // This is using the initial velocities
        const Cdot: number = Vec2.DotVV(J1, dp1) + Vec2.DotVV(J2, dp2) + Vec2.DotVV(J3, dp3);

        const B: number = C + alpha * c.lambda + sigma * Cdot;
        const sum2: number = (1.0 + sigma) * sum + alpha;

        const impulse: number = -B / sum2;

        // p1 += (c.invMass1 * impulse) * J1;
        p1.x += (c.invMass1 * impulse) * J1.x;
        p1.y += (c.invMass1 * impulse) * J1.y;
        // p2 += (c.invMass2 * impulse) * J2;
        p2.x += (c.invMass2 * impulse) * J2.x;
        p2.y += (c.invMass2 * impulse) * J2.y;
        // p3 += (c.invMass3 * impulse) * J3;
        p3.x += (c.invMass3 * impulse) * J3.x;
        p3.y += (c.invMass3 * impulse) * J3.y;

        this.ps[c.i1].Copy(p1);
        this.ps[c.i2].Copy(p2);
        this.ps[c.i3].Copy(p3);
        c.lambda += impulse;
      }
    }

    private SolveBend_PBD_Distance(): void {
      const stiffness: number = this.tuning.bendStiffness;

      for (let i = 0; i < this.bendCount; ++i) {
        const c: RopeBend = this.bendConstraints[i];

        const i1: number = c.i1;
        const i2: number = c.i3;

        const p1: Vec2 = this.ps[i1].Clone();
        const p2: Vec2 = this.ps[i2].Clone();

        // Vec2 d = p2 - p1;
        const d = p2.Clone().SelfSub(p1);
        const L: number = d.Normalize();

        const sum: number = c.invMass1 + c.invMass3;
        if (sum === 0.0) {
          continue;
        }

        const s1: number = c.invMass1 / sum;
        const s2: number = c.invMass3 / sum;

        // p1 -= stiffness * s1 * (c.L1 + c.L2 - L) * d;
        p1.x -= stiffness * s1 * (c.L1 + c.L2 - L) * d.x;
        p1.y -= stiffness * s1 * (c.L1 + c.L2 - L) * d.y;
        // p2 += stiffness * s2 * (c.L1 + c.L2 - L) * d;
        p2.x += stiffness * s2 * (c.L1 + c.L2 - L) * d.x;
        p2.y += stiffness * s2 * (c.L1 + c.L2 - L) * d.y;

        this.ps[i1].Copy(p1);
        this.ps[i2].Copy(p2);
      }
    }

    private SolveBend_PBD_Height(): void {
      const stiffness: number = this.tuning.bendStiffness;

      for (let i = 0; i < this.bendCount; ++i) {
        const c: RopeBend = this.bendConstraints[i];

        const p1: Vec2 = this.ps[c.i1].Clone();
        const p2: Vec2 = this.ps[c.i2].Clone();
        const p3: Vec2 = this.ps[c.i3].Clone();

        // Barycentric coordinates are held constant
        const d = new Vec2();
        // Vec2 d = c.alpha1 * p1 + c.alpha2 * p3 - p2;
        d.x = c.alpha1 * p1.x + c.alpha2 * p3.x - p2.x;
        d.y = c.alpha1 * p1.y + c.alpha2 * p3.y - p2.y;
        const dLen: number = d.Length();

        if (dLen === 0.0) {
          continue;
        }

        // Vec2 dHat = (1.0 / dLen) * d;
        const dHat = d.Clone().SelfMul(1.0 / dLen);

        // Vec2 J1 = c.alpha1 * dHat;
        const J1 = dHat.Clone().SelfMul(c.alpha1);
        // Vec2 J2 = -dHat;
        const J2 = dHat.Clone().SelfNeg();
        // Vec2 J3 = c.alpha2 * dHat;
        const J3 = dHat.Clone().SelfMul(c.alpha2);

        const sum: number = c.invMass1 * c.alpha1 * c.alpha1 + c.invMass2 + c.invMass3 * c.alpha2 * c.alpha2;

        if (sum === 0.0) {
          continue;
        }

        const C: number = dLen;
        const mass: number = 1.0 / sum;
        const impulse: number = -stiffness * mass * C;

        // p1 += (c.invMass1 * impulse) * J1;
        p1.x += (c.invMass1 * impulse) * J1.x;
        p1.y += (c.invMass1 * impulse) * J1.y;
        // p2 += (c.invMass2 * impulse) * J2;
        p2.x += (c.invMass2 * impulse) * J2.x;
        p2.y += (c.invMass2 * impulse) * J2.y;
        // p3 += (c.invMass3 * impulse) * J3;
        p3.x += (c.invMass3 * impulse) * J3.x;
        p3.y += (c.invMass3 * impulse) * J3.y;

        this.ps[c.i1].Copy(p1);
        this.ps[c.i2].Copy(p2);
        this.ps[c.i3].Copy(p3);
      }
    }

    // M. Kelager: A Triangle Bending Constraint Model for PBD
    private SolveBend_PBD_Triangle(): void {
      const stiffness = this.tuning.bendStiffness;

      for (let i = 0; i < this.bendCount; ++i) {
        const c: RopeBend = this.bendConstraints[i];

        const b0 = this.ps[c.i1].Clone();
        const v = this.ps[c.i2].Clone();
        const b1 = this.ps[c.i3].Clone();

        const wb0 = c.invMass1;
        const wv = c.invMass2;
        const wb1 = c.invMass3;

        const W = wb0 + wb1 + 2.0 * wv;
        const invW = stiffness / W;

        const d = new Vec2();
        d.x = v.x - (1.0 / 3.0) * (b0.x + v.x + b1.x);
        d.y = v.y - (1.0 / 3.0) * (b0.y + v.y + b1.y);

        const db0 = new Vec2();
        db0.x = 2.0 * wb0 * invW * d.x;
        db0.y = 2.0 * wb0 * invW * d.y;
        const dv = new Vec2();
        dv.x = -4.0 * wv * invW * d.x;
        dv.y = -4.0 * wv * invW * d.y;
        const db1 = new Vec2();
        db1.x = 2.0 * wb1 * invW * d.x;
        db1.y = 2.0 * wb1 * invW * d.y;

        b0.SelfAdd(db0);
        v.SelfAdd(dv);
        b1.SelfAdd(db1);

        this.ps[c.i1].Copy(b0);
        this.ps[c.i2].Copy(v);
        this.ps[c.i3].Copy(b1);
      }
    }

    private ApplyBendForces(dt: number): void {
      // omega = 2 * pi * hz
      const omega: number = 2.0 * pi * this.tuning.bendHertz;

      for (let i = 0; i < this.bendCount; ++i) {
        const c: RopeBend = this.bendConstraints[i];

        const p1: Vec2 = this.ps[c.i1].Clone();
        const p2: Vec2 = this.ps[c.i2].Clone();
        const p3: Vec2 = this.ps[c.i3].Clone();

        const v1: Vec2 = this.vs[c.i1];
        const v2: Vec2 = this.vs[c.i2];
        const v3: Vec2 = this.vs[c.i3];

        // Vec2 d1 = p2 - p1;
        const d1 = p1.Clone().SelfSub(p1);
        // Vec2 d2 = p3 - p2;
        const d2 = p3.Clone().SelfSub(p2);

        let L1sqr: number, L2sqr: number;

        if (this.tuning.isometric) {
          L1sqr = c.L1 * c.L1;
          L2sqr = c.L2 * c.L2;
        }
        else {
          L1sqr = d1.LengthSquared();
          L2sqr = d2.LengthSquared();
        }

        if (L1sqr * L2sqr === 0.0) {
          continue;
        }

        const a: number = Vec2.CrossVV(d1, d2);
        const b: number = Vec2.DotVV(d1, d2);

        const angle: number = Atan2(a, b);

        // Vec2 Jd1 = (-1.0 / L1sqr) * d1.Skew();
        // Vec2 Jd2 = (1.0 / L2sqr) * d2.Skew();

        // Vec2 J1 = -Jd1;
        // Vec2 J2 = Jd1 - Jd2;
        // Vec2 J3 = Jd2;

        // Vec2 Jd1 = (-1.0 / L1sqr) * d1.Skew();
        const Jd1: Vec2 = new Vec2().Copy(d1).SelfSkew().SelfMul(-1.0 / L1sqr);
        // Vec2 Jd2 = (1.0 / L2sqr) * d2.Skew();
        const Jd2: Vec2 = new Vec2().Copy(d2).SelfSkew().SelfMul(1.0 / L2sqr);

        // Vec2 J1 = -Jd1;
        const J1 = Jd1.Clone().SelfNeg();
        // Vec2 J2 = Jd1 - Jd2;
        const J2 = Jd1.Clone().SelfSub(Jd2);
        // Vec2 J3 = Jd2;
        const J3 = Jd2;

        let sum: number = 0.0;
        if (this.tuning.fixedEffectiveMass) {
          sum = c.invEffectiveMass;
        }
        else {
          sum = c.invMass1 * Vec2.DotVV(J1, J1) + c.invMass2 * Vec2.DotVV(J2, J2) + c.invMass3 * Vec2.DotVV(J3, J3);
        }

        if (sum === 0.0) {
          continue;
        }

        const mass: number = 1.0 / sum;

        const spring: number = mass * omega * omega;
        const damper: number = 2.0 * mass * this.tuning.bendDamping * omega;

        const C: number = angle;
        const Cdot: number = Vec2.DotVV(J1, v1) + Vec2.DotVV(J2, v2) + Vec2.DotVV(J3, v3);

        const impulse: number = -dt * (spring * C + damper * Cdot);

        // this.vs[c.i1] += (c.invMass1 * impulse) * J1;
        this.vs[c.i1].x += (c.invMass1 * impulse) * J1.x;
        this.vs[c.i1].y += (c.invMass1 * impulse) * J1.y;
        // this.vs[c.i2] += (c.invMass2 * impulse) * J2;
        this.vs[c.i2].x += (c.invMass2 * impulse) * J2.x;
        this.vs[c.i2].y += (c.invMass2 * impulse) * J2.y;
        // this.vs[c.i3] += (c.invMass3 * impulse) * J3;
        this.vs[c.i3].x += (c.invMass3 * impulse) * J3.x;
        this.vs[c.i3].y += (c.invMass3 * impulse) * J3.y;
      }
    }
  }

}
