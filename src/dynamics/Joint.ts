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
    e_unknownJoint = 0,
    e_revoluteJoint = 1,
    e_prismaticJoint = 2,
    e_distanceJoint = 3,
    e_pulleyJoint = 4,
    e_mouseJoint = 5,
    e_gearJoint = 6,
    e_wheelJoint = 7,
    e_weldJoint = 8,
    e_frictionJoint = 9,
    e_ropeJoint = 10,
    e_motorJoint = 11,
    e_areaJoint = 12,
  }

  export class Jacobian {
    public readonly linear: Vec2 = new Vec2();
    public angularA: number = 0;
    public angularB: number = 0;

    public SetZero(): Jacobian {
      this.linear.SetZero();
      this.angularA = 0;
      this.angularB = 0;
      return this;
    }

    public Set(x: XY, a1: number, a2: number): Jacobian {
      this.linear.Copy(x);
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
    private _other: Body | null = null; ///< provides quick access to the other body attached.
    public get other(): Body {
      if (this._other === null) { throw new Error(); }
      return this._other;
    }
    public set other(value: Body) {
      if (this._other !== null) { throw new Error(); }
      this._other = value;
    }
    public readonly joint: Joint;    ///< the joint
    public prev: JointEdge | null = null;  ///< the previous joint edge in the body's joint list
    public next: JointEdge | null = null;  ///< the next joint edge in the body's joint list
    constructor(joint: Joint) {
      this.joint = joint;
    }
    public Reset(): void {
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
    public readonly type: JointType = JointType.e_unknownJoint;

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
  export function LinearStiffness(def: { stiffness: number, damping: number }, frequencyHertz: number, dampingRatio: number, bodyA: Body, bodyB: Body): void {
    const massA: number = bodyA.GetMass();
    const massB: number = bodyB.GetMass();
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
  export function AngularStiffness(def: { stiffness: number, damping: number }, frequencyHertz: number, dampingRatio: number, bodyA: Body, bodyB: Body): void {
    const IA: number = bodyA.GetInertia();
    const IB: number = bodyB.GetInertia();
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
    public readonly type: JointType = JointType.e_unknownJoint;
    public prev: Joint | null = null;
    public next: Joint | null = null;
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

      this.collideConnected = Maybe(def.collideConnected, false);

      this.userData = Maybe(def.userData, null);
    }

    /// Get the type of the concrete joint.
    public GetType(): JointType {
      return this.type;
    }

    /// Get the first body attached to this joint.
    public GetBodyA(): Body {
      return this.bodyA;
    }

    /// Get the second body attached to this joint.
    public GetBodyB(): Body {
      return this.bodyB;
    }

    /// Get the anchor point on bodyA in world coordinates.
    public abstract GetAnchorA<T extends XY>(out: T): T;

    /// Get the anchor point on bodyB in world coordinates.
    public abstract GetAnchorB<T extends XY>(out: T): T;

    /// Get the reaction force on bodyB at the joint anchor in Newtons.
    public abstract GetReactionForce<T extends XY>(inv_dt: number, out: T): T;

    /// Get the reaction torque on bodyB in N*m.
    public abstract GetReactionTorque(inv_dt: number): number;

    /// Get the next joint the world joint list.
    public GetNext(): Joint | null {
      return this.next;
    }

    /// Get the user data pointer.
    public GetUserData(): any {
      return this.userData;
    }

    /// Set the user data pointer.
    public SetUserData(data: any): void {
      this.userData = data;
    }

    /// Short-cut function to determine if either body is inactive.
    public IsEnabled(): boolean {
      return this.bodyA.IsEnabled() && this.bodyB.IsEnabled();
    }

    /// Get collide connected.
    /// Note: modifying the collide connect flag won't work correctly because
    /// the flag is only checked when fixture AABBs begin to overlap.
    public GetCollideConnected(): boolean {
      return this.collideConnected;
    }

    /// Dump this joint to the log file.
    public Dump(log: (format: string, ...args: any[]) => void): void {
      log("// Dump is not supported for this joint type.\n");
    }

    /// Shift the origin for any points stored in world coordinates.
    public ShiftOrigin(newOrigin: XY): void { }

    /// Debug draw this joint
    private static Draw_s_p1: Vec2 = new Vec2();
    private static Draw_s_p2: Vec2 = new Vec2();
    private static Draw_s_color: Color = new Color(0.5, 0.8, 0.8);
    private static Draw_s_c: Color = new Color();
    public Draw(draw: Draw): void {
      const xf1: Transform = this.bodyA.GetTransform();
      const xf2: Transform = this.bodyB.GetTransform();
      const x1: Vec2 = xf1.p;
      const x2: Vec2 = xf2.p;
      const p1: Vec2 = this.GetAnchorA(Joint.Draw_s_p1);
      const p2: Vec2 = this.GetAnchorB(Joint.Draw_s_p2);

      const color: Color = Joint.Draw_s_color.SetRGB(0.5, 0.8, 0.8);

      switch (this.type) {
        case JointType.e_distanceJoint:
          draw.DrawSegment(p1, p2, color);
          break;

        case JointType.e_pulleyJoint:
        {
          const pulley: PulleyJoint = this as unknown as PulleyJoint;
          const s1: Vec2 = pulley.GetGroundAnchorA();
          const s2: Vec2 = pulley.GetGroundAnchorB();
          draw.DrawSegment(s1, p1, color);
          draw.DrawSegment(s2, p2, color);
          draw.DrawSegment(s1, s2, color);
        }
          break;

        case JointType.e_mouseJoint:
        {
          const c = Joint.Draw_s_c;
          c.Set(0.0, 1.0, 0.0);
          draw.DrawPoint(p1, 4.0, c);
          draw.DrawPoint(p2, 4.0, c);

          c.Set(0.8, 0.8, 0.8);
          draw.DrawSegment(p1, p2, c);
        }
          break;

        default:
          draw.DrawSegment(x1, p1, color);
          draw.DrawSegment(p1, p2, color);
          draw.DrawSegment(x2, p2, color);
      }
    }

    public abstract InitVelocityConstraints(data: SolverData): void;

    public abstract SolveVelocityConstraints(data: SolverData): void;

    // This returns true if the position errors are within tolerance.
    public abstract SolvePositionConstraints(data: SolverData): boolean;
  }

}
