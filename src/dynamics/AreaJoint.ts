namespace b2 {
  export interface IAreaJointDef extends IJointDef {
    // world: World;

    bodies: Body[];

    stiffness?: number;

    damping?: number;
  }

  export class AreaJointDef extends JointDef implements IAreaJointDef {
    public bodies: Body[] = [];

    public stiffness: number = 0;

    public damping: number = 0;

    constructor() {
      super(JointType.AreaJoint);
    }

    public addBody(body: Body): void {
      this.bodies.push(body);

      if (this.bodies.length === 1) {
        this.bodyA = body;
      } else if (this.bodies.length === 2) {
        this.bodyB = body;
      }
    }
  }

  export class AreaJoint extends Joint {
    public bodies: Body[];
    public stiffness: number = 0;
    public damping: number = 0;

    // Solver shared
    public impulse: number = 0;

    // Solver temp
    public readonly targetLengths: number[];
    public targetArea: number = 0;
    public readonly normals: Vec2[];
    public readonly joints: DistanceJoint[];
    public readonly deltas: Vec2[];
    public readonly delta: Vec2 = new Vec2();

    constructor(def: IAreaJointDef) {
      super(def);

      // DEBUG: Assert(def.bodies.length >= 3, "You cannot create an area joint with less than three bodies.");

      this.bodies = def.bodies;
      this.stiffness = maybe(def.stiffness, 0);
      this.damping = maybe(def.damping, 0);

      this.targetLengths = MakeNumberArray(def.bodies.length);
      this.normals = Vec2.MakeArray(def.bodies.length);
      this.joints = []; // MakeNullArray(def.bodies.length);
      this.deltas = Vec2.MakeArray(def.bodies.length);

      const djd: DistanceJointDef = new DistanceJointDef();
      djd.stiffness = this.stiffness;
      djd.damping = this.damping;

      this.targetArea = 0;

      for (let i: number = 0; i < this.bodies.length; ++i) {
        const body: Body = this.bodies[i];
        const next: Body = this.bodies[(i + 1) % this.bodies.length];

        const body_c: Vec2 = body.getWorldCenter();
        const next_c: Vec2 = next.getWorldCenter();

        this.targetLengths[i] = Vec2.DistanceVV(body_c, next_c);

        this.targetArea += Vec2.CrossVV(body_c, next_c);

        djd.Initialize(body, next, body_c, next_c);
        this.joints[i] = body.getWorld().createJoint(djd);
      }

      this.targetArea *= 0.5;
    }

    public getAnchorA<T extends XY>(out: T): T {
      return out;
    }

    public getAnchorB<T extends XY>(out: T): T {
      return out;
    }

    public getReactionForce<T extends XY>(inv_dt: number, out: T): T {
      return out;
    }

    public getReactionTorque(inv_dt: number): number {
      return 0;
    }

    public setStiffness(stiffness: number): void {
      this.stiffness = stiffness;

      for (let i: number = 0; i < this.joints.length; ++i) {
        this.joints[i].SetStiffness(stiffness);
      }
    }

    public getStiffness() {
      return this.stiffness;
    }

    public setDamping(damping: number): void {
      this.damping = damping;

      for (let i: number = 0; i < this.joints.length; ++i) {
        this.joints[i].SetDamping(damping);
      }
    }

    public getDamping() {
      return this.damping;
    }

    public dump(log: (format: string, ...args: any[]) => void) {
      log("Area joint dumping is not supported.\n");
    }

    public initVelocityConstraints(data: SolverData): void {
      for (let i: number = 0; i < this.bodies.length; ++i) {
        const prev: Body = this.bodies[(i + this.bodies.length - 1) % this.bodies.length];
        const next: Body = this.bodies[(i + 1) % this.bodies.length];
        const prev_c: Vec2 = data.positions[prev.islandIndex].c;
        const next_c: Vec2 = data.positions[next.islandIndex].c;
        const delta: Vec2 = this.deltas[i];

        Vec2.SubVV(next_c, prev_c, delta);
      }

      if (data.step.warmStarting) {
        this.impulse *= data.step.dtRatio;

        for (let i: number = 0; i < this.bodies.length; ++i) {
          const body: Body = this.bodies[i];
          const body_v: Vec2 = data.velocities[body.islandIndex].v;
          const delta: Vec2 = this.deltas[i];

          body_v.x += body.invMass *  delta.y * 0.5 * this.impulse;
          body_v.y += body.invMass * -delta.x * 0.5 * this.impulse;
        }
      } else {
        this.impulse = 0;
      }
    }

    public solveVelocityConstraints(data: SolverData): void {
      let dotMassSum: number = 0;
      let crossMassSum: number = 0;

      for (let i: number = 0; i < this.bodies.length; ++i) {
        const body: Body = this.bodies[i];
        const body_v: Vec2 = data.velocities[body.islandIndex].v;
        const delta: Vec2 = this.deltas[i];

        dotMassSum += delta.lengthSquared() / body.getMass();
        crossMassSum += Vec2.CrossVV(body_v, delta);
      }

      const lambda: number = -2 * crossMassSum / dotMassSum;
      // lambda = Clamp(lambda, -maxLinearCorrection, maxLinearCorrection);

      this.impulse += lambda;

      for (let i: number = 0; i < this.bodies.length; ++i) {
        const body: Body = this.bodies[i];
        const body_v: Vec2 = data.velocities[body.islandIndex].v;
        const delta: Vec2 = this.deltas[i];

        body_v.x += body.invMass *  delta.y * 0.5 * lambda;
        body_v.y += body.invMass * -delta.x * 0.5 * lambda;
      }
    }

    public solvePositionConstraints(data: SolverData): boolean {
      let perimeter: number = 0;
      let area: number = 0;

      for (let i: number = 0; i < this.bodies.length; ++i) {
        const body: Body = this.bodies[i];
        const next: Body = this.bodies[(i + 1) % this.bodies.length];
        const body_c: Vec2 = data.positions[body.islandIndex].c;
        const next_c: Vec2 = data.positions[next.islandIndex].c;

        const delta: Vec2 = Vec2.SubVV(next_c, body_c, this.delta);

        let dist: number = delta.length();
        if (dist < epsilon) {
          dist = 1;
        }

        this.normals[i].x =  delta.y / dist;
        this.normals[i].y = -delta.x / dist;

        perimeter += dist;

        area += Vec2.CrossVV(body_c, next_c);
      }

      area *= 0.5;

      const deltaArea: number = this.targetArea - area;
      const toExtrude: number = 0.5 * deltaArea / perimeter;
      let done: boolean = true;

      for (let i: number = 0; i < this.bodies.length; ++i) {
        const body: Body = this.bodies[i];
        const body_c: Vec2 = data.positions[body.islandIndex].c;
        const next_i: number = (i + 1) % this.bodies.length;

        const delta: Vec2 = Vec2.AddVV(this.normals[i], this.normals[next_i], this.delta);
        delta.selfMul(toExtrude);

        const norsq: number = delta.lengthSquared();
        if (norsq > sqrt(maxLinearCorrection)) {
          delta.selfMul(maxLinearCorrection / Sqrt(norsq));
        }
        if (norsq > sqrt(linearSlop)) {
          done = false;
        }

        body_c.x += delta.x;
        body_c.y += delta.y;
      }

      return done;
    }
  }

}
