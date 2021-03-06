/*
* Copyright (c) 2006-2010 Erin Catto http://www.box2d.org
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

/// A line segment (edge) shape. These can be connected in chains or loops
/// to other edge shapes. Edges created independently are two-sided and do
/// no provide smooth movement across junctions.
  export class EdgeShape extends Shape {
    public readonly vertex0: Vec2 = new Vec2();
    public readonly vertex1: Vec2 = new Vec2();
    public readonly vertex2: Vec2 = new Vec2();
    public readonly vertex3: Vec2 = new Vec2();

    /// Uses vertex0 and vertex3 to create smooth collision.
    public oneSided: boolean = false;

    constructor() {
      super(ShapeType.EdgeShape, polygonRadius);
    }

    /// Set this as a part of a sequence. Vertex v0 precedes the edge and vertex v3
    /// follows. These extra vertices are used to provide smooth movement
    /// across junctions. This also makes the collision one-sided. The edge
    /// normal points to the right looking from v1 to v2.
    // void SetOneSided(const Vec2& v0, const Vec2& v1,const Vec2& v2, const Vec2& v3);
    public setOneSided(v0: XY, v1: XY, v2: XY, v3: XY): EdgeShape {
      this.vertex0.copy(v0);
      this.vertex1.copy(v1);
      this.vertex2.copy(v2);
      this.vertex3.copy(v3);
      this.oneSided = true;
      return this;
    }

    /// Set this as an isolated edge. Collision is two-sided.
    public setTwoSided(v1: XY, v2: XY): EdgeShape {
      this.vertex1.copy(v1);
      this.vertex2.copy(v2);
      this.oneSided = false;
      return this;
    }

    /// Implement Shape.
    public clone(): EdgeShape {
      return new EdgeShape().copy(this);
    }

    public copy(other: EdgeShape): EdgeShape {
      super.copy(other);

      // DEBUG: Assert(other instanceof EdgeShape);

      this.vertex1.copy(other.vertex1);
      this.vertex2.copy(other.vertex2);
      this.vertex0.copy(other.vertex0);
      this.vertex3.copy(other.vertex3);
      this.oneSided = other.oneSided;

      return this;
    }

    /// @see Shape::GetChildCount
    public getChildCount(): number {
      return 1;
    }

    /// @see Shape::TestPoint
    public testPoint(xf: Transform, p: XY): boolean {
      return false;
    }

    // #if ENABLE_PARTICLE
    /// @see Shape::ComputeDistance
    private static ComputeDistance_s_v1 = new Vec2();
    private static ComputeDistance_s_v2 = new Vec2();
    private static ComputeDistance_s_d = new Vec2();
    private static ComputeDistance_s_s = new Vec2();
    public computeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number {
      const v1 = Transform.mulXV(xf, this.vertex1, EdgeShape.ComputeDistance_s_v1);
      const v2 = Transform.mulXV(xf, this.vertex2, EdgeShape.ComputeDistance_s_v2);

      const d = Vec2.SubVV(p, v1, EdgeShape.ComputeDistance_s_d);
      const s = Vec2.SubVV(v2, v1, EdgeShape.ComputeDistance_s_s);
      const ds = Vec2.DotVV(d, s);
      if (ds > 0) {
        const s2 = Vec2.DotVV(s, s);
        if (ds > s2) {
          Vec2.SubVV(p, v2, d);
        } else {
          d.selfMulSub(ds / s2, s);
        }
      }
      normal.copy(d);
      return normal.normalize();
    }
    // #endif

    /// Implement Shape.
    // p = p1 + t * d
    // v = v1 + s * e
    // p1 + t * d = v1 + s * e
    // s * e - t * d = p1 - v1
    private static rayCast_s_p1 = new Vec2();
    private static rayCast_s_p2 = new Vec2();
    private static rayCast_s_d = new Vec2();
    private static rayCast_s_e = new Vec2();
    private static rayCast_s_q = new Vec2();
    private static rayCast_s_r = new Vec2();
    public rayCast(output: RayCastOutput, input: RayCastInput, xf: Transform, childIndex: number): boolean {
      // Put the ray into the edge's frame of reference.
      const p1: Vec2 = Transform.mulTXV(xf, input.p1, EdgeShape.rayCast_s_p1);
      const p2: Vec2 = Transform.mulTXV(xf, input.p2, EdgeShape.rayCast_s_p2);
      const d: Vec2 = Vec2.SubVV(p2, p1, EdgeShape.rayCast_s_d);

      const v1: Vec2 = this.vertex1;
      const v2: Vec2 = this.vertex2;
      const e: Vec2 = Vec2.SubVV(v2, v1, EdgeShape.rayCast_s_e);

      // Normal points to the right, looking from v1 at v2
      const normal: Vec2 = output.normal.set(e.y, -e.x).selfNormalize();

      // q = p1 + t * d
      // dot(normal, q - v1) = 0
      // dot(normal, p1 - v1) + t * dot(normal, d) = 0
      const numerator: number = Vec2.DotVV(normal, Vec2.SubVV(v1, p1, Vec2.s_t0));
      if (this.oneSided && numerator > 0.0) {
        return false;
      }

      const denominator: number = Vec2.DotVV(normal, d);

      if (denominator === 0) {
        return false;
      }

      const t: number = numerator / denominator;
      if (t < 0 || input.maxFraction < t) {
        return false;
      }

      const q: Vec2 = Vec2.AddVMulSV(p1, t, d, EdgeShape.rayCast_s_q);

      // q = v1 + s * r
      // s = dot(q - v1, r) / dot(r, r)
      const r: Vec2 = Vec2.SubVV(v2, v1, EdgeShape.rayCast_s_r);
      const rr: number = Vec2.DotVV(r, r);
      if (rr === 0) {
        return false;
      }

      const s: number = Vec2.DotVV(Vec2.SubVV(q, v1, Vec2.s_t0), r) / rr;
      if (s < 0 || 1 < s) {
        return false;
      }

      output.fraction = t;
      Rot.mulRV(xf.q, output.normal, output.normal);
      if (numerator > 0) {
        output.normal.selfNeg();
      }
      return true;
    }

    /// @see Shape::ComputeAABB
    private static computeAABB_s_v1 = new Vec2();
    private static computeAABB_s_v2 = new Vec2();
    public computeAABB(aabb: AABB, xf: Transform, childIndex: number): void {
      const v1: Vec2 = Transform.mulXV(xf, this.vertex1, EdgeShape.computeAABB_s_v1);
      const v2: Vec2 = Transform.mulXV(xf, this.vertex2, EdgeShape.computeAABB_s_v2);

      Vec2.MinV(v1, v2, aabb.lowerBound);
      Vec2.MaxV(v1, v2, aabb.upperBound);

      const r: number = this.radius;
      aabb.lowerBound.selfSubXY(r, r);
      aabb.upperBound.selfAddXY(r, r);
    }

    /// @see Shape::ComputeMass
    public computeMass(massData: MassData, density: number): void {
      massData.mass = 0;
      Vec2.MidVV(this.vertex1, this.vertex2, massData.center);
      massData.I = 0;
    }

    public setupDistanceProxy(proxy: DistanceProxy, index: number): void {
      proxy.vertices = proxy.buffer;
      proxy.vertices[0].copy(this.vertex1);
      proxy.vertices[1].copy(this.vertex2);
      proxy.count = 2;
      proxy.radius = this.radius;
    }

    public computeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number {
      c.setZero();
      return 0;
    }

    public dump(log: (format: string, ...args: any[]) => void): void {
      log("    const shape: EdgeShape = new EdgeShape();\n");
      log("    shape.radius = %.15f;\n", this.radius);
      log("    shape.vertex0.Set(%.15f, %.15f);\n", this.vertex0.x, this.vertex0.y);
      log("    shape.vertex1.Set(%.15f, %.15f);\n", this.vertex1.x, this.vertex1.y);
      log("    shape.vertex2.Set(%.15f, %.15f);\n", this.vertex2.x, this.vertex2.y);
      log("    shape.vertex3.Set(%.15f, %.15f);\n", this.vertex3.x, this.vertex3.y);
      log("    shape.oneSided = %s;\n", this.oneSided);
    }
  }

}
