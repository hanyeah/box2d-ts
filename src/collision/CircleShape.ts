/*
* Copyright (c) 2006-2009 Erin Catto http://www.box2d.org
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

// DEBUG: import { Assert } from "../common/settings";

/// A solid circle shape
namespace b2 {
  export class CircleShape extends Shape {
    public readonly p: Vec2 = new Vec2();

    constructor(radius: number = 0) {
      super(ShapeType.e_circleShape, radius);
    }

    public Set(position: XY, radius: number = this.radius): this {
      this.p.Copy(position);
      this.radius = radius;
      return this;
    }

    /// Implement Shape.
    public Clone(): CircleShape {
      return new CircleShape().Copy(this);
    }

    public Copy(other: CircleShape): CircleShape {
      super.Copy(other);

      // DEBUG: Assert(other instanceof CircleShape);

      this.p.Copy(other.p);
      return this;
    }

    /// @see Shape::GetChildCount
    public GetChildCount(): number {
      return 1;
    }

    /// Implement Shape.
    private static TestPoint_s_center = new Vec2();
    private static TestPoint_s_d = new Vec2();
    public TestPoint(transform: Transform, p: XY): boolean {
      const center: Vec2 = Transform.MulXV(transform, this.p, CircleShape.TestPoint_s_center);
      const d: Vec2 = Vec2.SubVV(p, center, CircleShape.TestPoint_s_d);
      return Vec2.DotVV(d, d) <= Sq(this.radius);
    }

    // #if ENABLE_PARTICLE
    /// @see Shape::ComputeDistance
    private static ComputeDistance_s_center = new Vec2();
    public ComputeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number {
      const center = Transform.MulXV(xf, this.p, CircleShape.ComputeDistance_s_center);
      Vec2.SubVV(p, center, normal);
      return normal.Normalize() - this.radius;
    }
    // #endif

    /// Implement Shape.
    /// @note because the circle is solid, rays that start inside do not hit because the normal is
    /// not defined.
    // Collision Detection in Interactive 3D Environments by Gino van den Bergen
    // From Section 3.1.2
    // x = s + a * r
    // norm(x) = radius
    private static RayCast_s_position = new Vec2();
    private static RayCast_s_s = new Vec2();
    private static RayCast_s_r = new Vec2();
    public RayCast(output: RayCastOutput, input: RayCastInput, transform: Transform, childIndex: number): boolean {
      const position: Vec2 = Transform.MulXV(transform, this.p, CircleShape.RayCast_s_position);
      const s: Vec2 = Vec2.SubVV(input.p1, position, CircleShape.RayCast_s_s);
      const b: number = Vec2.DotVV(s, s) - Sq(this.radius);

      // Solve quadratic equation.
      const r: Vec2 = Vec2.SubVV(input.p2, input.p1, CircleShape.RayCast_s_r);
      const c: number = Vec2.DotVV(s, r);
      const rr: number = Vec2.DotVV(r, r);
      const sigma = c * c - rr * b;

      // Check for negative discriminant and short segment.
      if (sigma < 0 || rr < epsilon) {
        return false;
      }

      // Find the point of intersection of the line with the circle.
      let a: number = (-(c + Sqrt(sigma)));

      // Is the intersection point on the segment?
      if (0 <= a && a <= input.maxFraction * rr) {
        a /= rr;
        output.fraction = a;
        Vec2.AddVMulSV(s, a, r, output.normal).SelfNormalize();
        return true;
      }

      return false;
    }

    /// @see Shape::ComputeAABB
    private static ComputeAABB_s_p = new Vec2();
    public ComputeAABB(aabb: AABB, transform: Transform, childIndex: number): void {
      const p: Vec2 = Transform.MulXV(transform, this.p, CircleShape.ComputeAABB_s_p);
      aabb.lowerBound.Set(p.x - this.radius, p.y - this.radius);
      aabb.upperBound.Set(p.x + this.radius, p.y + this.radius);
    }

    /// @see Shape::ComputeMass
    public ComputeMass(massData: MassData, density: number): void {
      const radius_sq: number = Sq(this.radius);
      massData.mass = density * pi * radius_sq;
      massData.center.Copy(this.p);

      // inertia about the local origin
      massData.I = massData.mass * (0.5 * radius_sq + Vec2.DotVV(this.p, this.p));
    }

    public SetupDistanceProxy(proxy: DistanceProxy, index: number): void {
      proxy.vertices = proxy.buffer;
      proxy.vertices[0].Copy(this.p);
      proxy.count = 1;
      proxy.radius = this.radius;
    }

    public ComputeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number {
      const p: Vec2 = Transform.MulXV(xf, this.p, new Vec2());
      const l: number = (-(Vec2.DotVV(normal, p) - offset));

      if (l < (-this.radius) + epsilon) {
        // Completely dry
        return 0;
      }
      if (l > this.radius) {
        // Completely wet
        c.Copy(p);
        return pi * this.radius * this.radius;
      }

      // Magic
      const r2: number = this.radius * this.radius;
      const l2: number = l * l;
      const area: number = r2 * (Asin(l / this.radius) + pi / 2) + l * Sqrt(r2 - l2);
      const com: number = (-2 / 3 * Pow(r2 - l2, 1.5) / area);

      c.x = p.x + normal.x * com;
      c.y = p.y + normal.y * com;

      return area;
    }

    public Dump(log: (format: string, ...args: any[]) => void): void {
      log("    const shape: CircleShape = new CircleShape();\n");
      log("    shape.radius = %.15f;\n", this.radius);
      log("    shape.p.Set(%.15f, %.15f);\n", this.p.x, this.p.y);
    }
  }
}

