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
  export enum JointType {
    UnknownJoint = 0,
    RevoluteJoint = 1,
    PrismaticJoint = 2,
    DistanceJoint = 3,
    PulleyJoint = 4,
    MouseJoint = 5,
    GearJoint = 6,
    WheelJoint = 7,
    WeldJoint = 8,
    FrictionJoint = 9,
    RopeJoint = 10,
    MotorJoint = 11,
    AreaJoint = 12,
  }

  export class Jacobian {
    public readonly linear: Vec2 = new Vec2();
    public angularA: number = 0;
    public angularB: number = 0;

    public setZero(): Jacobian {
      this.linear.setZero();
      this.angularA = 0;
      this.angularB = 0;
      return this;
    }

    public set(x: XY, a1: number, a2: number): Jacobian {
      this.linear.copy(x);
      this.angularA = a1;
      this.angularB = a2;
      return this;
    }
  }

/// A joint edge is used to connect bodies and joints together
/// in a joint graph where each body is a node and each joint
/// is an edge. A joint edge belongs to a doubly linked list
/// maintained in each attached body. Each joint has two joint
/// nodes, one for each attached body.
  export class JointEdge {
    private _other: Body = null; ///< provides quick access to the other body attached.
    public get other(): Body {
      if (this._other === null) { throw new Error(); }
      return this._other;
    }
    public set other(value: Body) {
      if (this._other !== null) { throw new Error(); }
      this._other = value;
    }
    public readonly joint: Joint;    ///< the joint
    public prev: JointEdge = null;  ///< the previous joint edge in the body's joint list
    public next: JointEdge = null;  ///< the next joint edge in the body's joint list
    constructor(joint: Joint) {
      this.joint = joint;
    }
    public reset(): void {
      this._other = null;
      this.prev = null;
      this.next = null;
    }
  }

/// Joint definitions are used to construct joints.
  export interface IJointDef {
    /// The joint type is set automatically for concrete joint types.
    type: JointType;

    /// Use this to attach application specific data to your joints.
    userData?: any;

    /// The first attached body.
    bodyA: Body;

    /// The second attached body.
    bodyB: Body;

    /// Set this flag to true if the attached bodies should collide.
    collideConnected?: boolean;
  }

/// Joint definitions are used to construct joints.
  export abstract class JointDef implements IJointDef {
    /// The joint type is set automatically for concrete joint types.
    public readonly type: JointType = JointType.UnknownJoint;

    /// Use this to attach application specific data to your joints.
    public userData: any = null;

    /// The first attached body.
    public bodyA!: Body;

    /// The second attached body.
    public bodyB!: Body;

    /// Set this flag to true if the attached bodies should collide.
    public collideConnected: boolean = false;

    constructor(type: JointType) {
      this.type = type;
    }
  }

/// Utility to compute linear stiffness values from frequency and damping ratio
// void LinearStiffness(float& stiffness, float& damping,
// 	float frequencyHertz, float dampingRatio,
// 	const Body* bodyA, const Body* bodyB);
  export function linearStiffness(def: { stiffness: number, damping: number }, frequencyHertz: number, dampingRatio: number, bodyA: Body, bodyB: Body): void {
    const massA: number = bodyA.getMass();
    const massB: number = bodyB.getMass();
    let mass: number;
    if (massA > 0.0 && massB > 0.0) {
      mass = massA * massB / (massA + massB);
    } else if (massA > 0.0) {
      mass = massA;
    } else {
      mass = massB;
    }

    const omega: number = 2.0 * pi * frequencyHertz;
    def.stiffness = mass * omega * omega;
    def.damping = 2.0 * mass * dampingRatio * omega;
  }

/// Utility to compute rotational stiffness values frequency and damping ratio
// void AngularStiffness(float& stiffness, float& damping,
// 	float frequencyHertz, float dampingRatio,
// 	const Body* bodyA, const Body* bodyB);
  export function angularStiffness(def: { stiffness: number, damping: number }, frequencyHertz: number, dampingRatio: number, bodyA: Body, bodyB: Body): void {
    const IA: number = bodyA.getInertia();
    const IB: number = bodyB.getInertia();
    let I: number;
    if (IA > 0.0 && IB > 0.0) {
      I = IA * IB / (IA + IB);
    } else if (IA > 0.0) {
      I = IA;
    } else {
      I = IB;
    }

    const omega: number = 2.0 * pi * frequencyHertz;
    def.stiffness = I * omega * omega;
    def.damping = 2.0 * I * dampingRatio * omega;
  }

/// The base joint class. Joints are used to constraint two bodies together in
/// various fashions. Some joints also feature limits and motors.
  export abstract class Joint {
    public readonly type: JointType = JointType.UnknownJoint;
    public prev: Joint = null;
    public next: Joint = null;
    public readonly edgeA: JointEdge = new JointEdge(this);
    public readonly edgeB: JointEdge = new JointEdge(this);
    public bodyA: Body;
    public bodyB: Body;

    public index: number = 0;

    public islandFlag: boolean = false;
    public collideConnected: boolean = false;

    public userData: any = null;

    constructor(def: IJointDef) {
      // DEBUG: Assert(def.bodyA !== def.bodyB);

      this.type = def.type;
      this.edgeA.other = def.bodyB;
      this.edgeB.other = def.bodyA;
      this.bodyA = def.bodyA;
      this.bodyB = def.bodyB;

      this.collideConnected = maybe(def.collideConnected, false);

      this.userData = maybe(def.userData, null);
    }

    /// Get the type of the concrete joint.
    public getType(): JointType {
      return this.type;
    }

    /// Get the first body attached to this joint.
    public getBodyA(): Body {
      return this.bodyA;
    }

    /// Get the second body attached to this joint.
    public getBodyB(): Body {
      return this.bodyB;
    }

    /// Get the anchor point on bodyA in world coordinates.
    public abstract getAnchorA<T extends XY>(out: T): T;

    /// Get the anchor point on bodyB in world coordinates.
    public abstract getAnchorB<T extends XY>(out: T): T;

    /// Get the reaction force on bodyB at the joint anchor in Newtons.
    public abstract getReactionForce<T extends XY>(inv_dt: number, out: T): T;

    /// Get the reaction torque on bodyB in N*m.
    public abstract getReactionTorque(inv_dt: number): number;

    /// Short-cut function to determine if either body is inactive.
    public isEnabled(): boolean {
      return this.bodyA.isEnabled() && this.bodyB.isEnabled();
    }

    /// Get collide connected.
    /// Note: modifying the collide connect flag won't work correctly because
    /// the flag is only checked when fixture AABBs begin to overlap.
    public getCollideConnected(): boolean {
      return this.collideConnected;
    }

    /// Dump this joint to the log file.
    public dump(log: (format: string, ...args: any[]) => void): void {
      log("// Dump is not supported for this joint type.\n");
    }

    /// Shift the origin for any points stored in world coordinates.
    public shiftOrigin(newOrigin: XY): void { }

    /// Debug draw this joint
    private static draw_s_p1: Vec2 = new Vec2();
    private static draw_s_p2: Vec2 = new Vec2();
    private static draw_s_color: Color = new Color(0.5, 0.8, 0.8);
    private static draw_s_c: Color = new Color();
    public draw(draw: Draw): void {
      const xf1: Transform = this.bodyA.getTransform();
      const xf2: Transform = this.bodyB.getTransform();
      const x1: Vec2 = xf1.p;
      const x2: Vec2 = xf2.p;
      const p1: Vec2 = this.getAnchorA(Joint.draw_s_p1);
      const p2: Vec2 = this.getAnchorB(Joint.draw_s_p2);

      const color: Color = Joint.draw_s_color.setRGB(0.5, 0.8, 0.8);

      switch (this.type) {
        case JointType.DistanceJoint:
          draw.drawSegment(p1, p2, color);
          break;

        case JointType.PulleyJoint:
        {
          const pulley: PulleyJoint = this as unknown as PulleyJoint;
          const s1: Vec2 = pulley.getGroundAnchorA();
          const s2: Vec2 = pulley.getGroundAnchorB();
          draw.drawSegment(s1, p1, color);
          draw.drawSegment(s2, p2, color);
          draw.drawSegment(s1, s2, color);
        }
          break;

        case JointType.MouseJoint:
        {
          const c = Joint.draw_s_c;
          c.set(0.0, 1.0, 0.0);
          draw.drawPoint(p1, 4.0, c);
          draw.drawPoint(p2, 4.0, c);

          c.set(0.8, 0.8, 0.8);
          draw.drawSegment(p1, p2, c);
        }
          break;

        default:
          draw.drawSegment(x1, p1, color);
          draw.drawSegment(p1, p2, color);
          draw.drawSegment(x2, p2, color);
      }
    }

    public abstract initVelocityConstraints(data: SolverData): void;

    public abstract solveVelocityConstraints(data: SolverData): void;

    // This returns true if the position errors are within tolerance.
    public abstract solvePositionConstraints(data: SolverData): boolean;
  }

}
