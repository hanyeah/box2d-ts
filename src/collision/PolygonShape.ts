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

namespace b2 {

/// A solid convex polygon. It is assumed that the interior of the polygon is to
/// the left of each edge.
/// In most cases you should not need many vertices for a convex polygon.
  export class PolygonShape extends Shape {
    public readonly centroid: Vec2 = new Vec2(0, 0);
    public vertices: Vec2[] = [];
    public normals: Vec2[] = [];
    public count: number = 0;

    constructor() {
      super(ShapeType.PolygonShape, polygonRadius);
    }

    /// Implement Shape.
    public clone(): PolygonShape {
      return new PolygonShape().copy(this);
    }

    public copy(other: PolygonShape): PolygonShape {
      super.copy(other);

      // DEBUG: Assert(other instanceof PolygonShape);

      this.centroid.copy(other.centroid);
      this.count = other.count;
      this.vertices = Vec2.MakeArray(this.count);
      this.normals = Vec2.MakeArray(this.count);
      for (let i: number = 0; i < this.count; ++i) {
        this.vertices[i].copy(other.vertices[i]);
        this.normals[i].copy(other.normals[i]);
      }
      return this;
    }

    /// @see Shape::GetChildCount
    public getChildCount(): number {
      return 1;
    }

    /// Create a convex hull from the given array of points.
    /// @warning the points may be re-ordered, even if they form a convex polygon
    /// @warning collinear points are handled but not removed. Collinear points
    /// may lead to poor stacking behavior.
    private static Set_s_r = new Vec2();
    private static Set_s_v = new Vec2();
    public set(vertices: XY[]): PolygonShape;
    public set(vertices: XY[], count: number): PolygonShape;
    public set(vertices: number[]): PolygonShape;
    public set(...args: any[]): PolygonShape {
      if (typeof args[0][0] === "number") {
        const vertices: number[] = args[0];
        if (vertices.length % 2 !== 0) { throw new Error(); }
        return this._set((index: number): XY => ({ x: vertices[index * 2], y: vertices[index * 2 + 1] }), vertices.length / 2);
      } else {
        const vertices: XY[] = args[0];
        const count: number = args[1] || vertices.length;
        return this._set((index: number): XY => vertices[index], count);
      }
    }
    public _set(vertices: (index: number) => XY, count: number): PolygonShape {

      // DEBUG: Assert(3 <= count);
      if (count < 3) {
        return this.setAsBox(1, 1);
      }

      let n: number = count;

      // Perform welding and copy vertices into local buffer.
      const ps: XY[] = [];
      for (let i = 0; i < n; ++i) {
        const /*Vec2*/ v = vertices(i);

        let /*bool*/ unique = true;
        for (let /*int32*/ j = 0; j < ps.length; ++j) {
          if (Vec2.DistanceSquaredVV(v, ps[j]) < ((0.5 * linearSlop) * (0.5 * linearSlop))) {
            unique = false;
            break;
          }
        }

        if (unique) {
          ps.push(v);
        }
      }

      n = ps.length;
      if (n < 3) {
        // Polygon is degenerate.
        // DEBUG: Assert(false);
        return this.setAsBox(1.0, 1.0);
      }

      // Create the convex hull using the Gift wrapping algorithm
      // http://en.wikipedia.org/wiki/Gift_wrapping_algorithm

      // Find the right most point on the hull
      let i0: number = 0;
      let x0: number = ps[0].x;
      for (let i: number = 1; i < n; ++i) {
        const x: number = ps[i].x;
        if (x > x0 || (x === x0 && ps[i].y < ps[i0].y)) {
          i0 = i;
          x0 = x;
        }
      }

      const hull: number[] = [];
      let m: number = 0;
      let ih: number = i0;

      for (; ;) {
        hull[m] = ih;

        let ie: number = 0;
        for (let j: number = 1; j < n; ++j) {
          if (ie === ih) {
            ie = j;
            continue;
          }

          const r: Vec2 = Vec2.SubVV(ps[ie], ps[hull[m]], PolygonShape.Set_s_r);
          const v: Vec2 = Vec2.SubVV(ps[j], ps[hull[m]], PolygonShape.Set_s_v);
          const c: number = Vec2.CrossVV(r, v);
          if (c < 0) {
            ie = j;
          }

          // Collinearity check
          if (c === 0 && v.lengthSquared() > r.lengthSquared()) {
            ie = j;
          }
        }

        ++m;
        ih = ie;

        if (ie === i0) {
          break;
        }
      }

      this.count = m;
      this.vertices = Vec2.MakeArray(this.count);
      this.normals = Vec2.MakeArray(this.count);

      // Copy vertices.
      for (let i: number = 0; i < m; ++i) {
        this.vertices[i].copy(ps[hull[i]]);
      }

      // Compute normals. Ensure the edges have non-zero length.
      for (let i: number = 0; i < m; ++i) {
        const vertexi1: Vec2 = this.vertices[i];
        const vertexi2: Vec2 = this.vertices[(i + 1) % m];
        const edge: Vec2 = Vec2.SubVV(vertexi2, vertexi1, Vec2.s_t0); // edge uses s_t0
        // DEBUG: Assert(edge.LengthSquared() > epsilon_sq);
        Vec2.CrossVOne(edge, this.normals[i]).selfNormalize();
      }

      // Compute the polygon centroid.
      PolygonShape.computeCentroid(this.vertices, m, this.centroid);

      return this;
    }

    /// Build vertices to represent an axis-aligned box or an oriented box.
    /// @param hx the half-width.
    /// @param hy the half-height.
    /// @param center the center of the box in local coordinates.
    /// @param angle the rotation of the box in local coordinates.
    public setAsBox(hx: number, hy: number, center?: XY, angle: number = 0): PolygonShape {
      this.count = 4;
      this.vertices = Vec2.MakeArray(this.count);
      this.normals = Vec2.MakeArray(this.count);
      this.vertices[0].set((-hx), (-hy));
      this.vertices[1].set(hx, (-hy));
      this.vertices[2].set(hx, hy);
      this.vertices[3].set((-hx), hy);
      this.normals[0].set(0, (-1));
      this.normals[1].set(1, 0);
      this.normals[2].set(0, 1);
      this.normals[3].set((-1), 0);
      this.centroid.setZero();

      if (center) {
        this.centroid.copy(center);

        const xf: Transform = new Transform();
        xf.setPosition(center);
        xf.setRotationAngle(angle);

        // Transform vertices and normals.
        for (let i: number = 0; i < this.count; ++i) {
          Transform.mulXV(xf, this.vertices[i], this.vertices[i]);
          Rot.mulRV(xf.q, this.normals[i], this.normals[i]);
        }
      }

      return this;
    }

    /// @see Shape::TestPoint
    private static testPointSPLocal = new Vec2();
    public testPoint(xf: Transform, p: XY): boolean {
      const pLocal: Vec2 = Transform.mulTXV(xf, p, PolygonShape.testPointSPLocal);

      for (let i: number = 0; i < this.count; ++i) {
        const dot: number = Vec2.DotVV(this.normals[i], Vec2.SubVV(pLocal, this.vertices[i], Vec2.s_t0));
        if (dot > 0) {
          return false;
        }
      }

      return true;
    }

    // #if ENABLE_PARTICLE
    /// @see Shape::ComputeDistance
    private static computeDistance_s_pLocal = new Vec2();
    private static computeDistance_s_normalForMaxDistance = new Vec2();
    private static computeDistance_s_minDistance = new Vec2();
    private static computeDistance_s_distance = new Vec2();
    public computeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number {
      const pLocal = Transform.mulTXV(xf, p, PolygonShape.computeDistance_s_pLocal);
      let maxDistance = -maxFloat;
      const normalForMaxDistance = PolygonShape.computeDistance_s_normalForMaxDistance.copy(pLocal);

      for (let i = 0; i < this.count; ++i) {
        const dot = Vec2.DotVV(this.normals[i], Vec2.SubVV(pLocal, this.vertices[i], Vec2.s_t0));
        if (dot > maxDistance) {
          maxDistance = dot;
          normalForMaxDistance.copy(this.normals[i]);
        }
      }

      if (maxDistance > 0) {
        const minDistance = PolygonShape.computeDistance_s_minDistance.copy(normalForMaxDistance);
        let minDistance2 = maxDistance * maxDistance;
        for (let i = 0; i < this.count; ++i) {
          const distance = Vec2.SubVV(pLocal, this.vertices[i], PolygonShape.computeDistance_s_distance);
          const distance2 = distance.lengthSquared();
          if (minDistance2 > distance2) {
            minDistance.copy(distance);
            minDistance2 = distance2;
          }
        }

        Rot.mulRV(xf.q, minDistance, normal);
        normal.normalize();
        return Math.sqrt(minDistance2);
      } else {
        Rot.mulRV(xf.q, normalForMaxDistance, normal);
        return maxDistance;
      }
    }
    // #endif

    /// Implement Shape.
    /// @note because the polygon is solid, rays that start inside do not hit because the normal is
    /// not defined.
    private static rayCast_s_p1 = new Vec2();
    private static rayCast_s_p2 = new Vec2();
    private static rayCast_s_d = new Vec2();
    public rayCast(output: RayCastOutput, input: RayCastInput, xf: Transform, childIndex: number): boolean {
      // Put the ray into the polygon's frame of reference.
      const p1: Vec2 = Transform.mulTXV(xf, input.p1, PolygonShape.rayCast_s_p1);
      const p2: Vec2 = Transform.mulTXV(xf, input.p2, PolygonShape.rayCast_s_p2);
      const d: Vec2 = Vec2.SubVV(p2, p1, PolygonShape.rayCast_s_d);

      let lower: number = 0, upper = input.maxFraction;

      let index: number = -1;

      for (let i: number = 0; i < this.count; ++i) {
        // p = p1 + a * d
        // dot(normal, p - v) = 0
        // dot(normal, p1 - v) + a * dot(normal, d) = 0
        const numerator: number = Vec2.DotVV(this.normals[i], Vec2.SubVV(this.vertices[i], p1, Vec2.s_t0));
        const denominator: number = Vec2.DotVV(this.normals[i], d);

        if (denominator === 0) {
          if (numerator < 0) {
            return false;
          }
        } else {
          // Note: we want this predicate without division:
          // lower < numerator / denominator, where denominator < 0
          // Since denominator < 0, we have to flip the inequality:
          // lower < numerator / denominator <==> denominator * lower > numerator.
          if (denominator < 0 && numerator < lower * denominator) {
            // Increase lower.
            // The segment enters this half-space.
            lower = numerator / denominator;
            index = i;
          } else if (denominator > 0 && numerator < upper * denominator) {
            // Decrease upper.
            // The segment exits this half-space.
            upper = numerator / denominator;
          }
        }

        // The use of epsilon here causes the assert on lower to trip
        // in some cases. Apparently the use of epsilon was to make edge
        // shapes work, but now those are handled separately.
        // if (upper < lower - epsilon)
        if (upper < lower) {
          return false;
        }
      }

      // DEBUG: Assert(0 <= lower && lower <= input.maxFraction);

      if (index >= 0) {
        output.fraction = lower;
        Rot.mulRV(xf.q, this.normals[index], output.normal);
        return true;
      }

      return false;
    }

    /// @see Shape::ComputeAABB
    private static computeAABB_s_v = new Vec2();
    public computeAABB(aabb: AABB, xf: Transform, childIndex: number): void {
      const lower: Vec2 = Transform.mulXV(xf, this.vertices[0], aabb.lowerBound);
      const upper: Vec2 = aabb.upperBound.copy(lower);

      for (let i: number = 0; i < this.count; ++i) {
        const v: Vec2 = Transform.mulXV(xf, this.vertices[i], PolygonShape.computeAABB_s_v);
        Vec2.MinV(v, lower, lower);
        Vec2.MaxV(v, upper, upper);
      }

      const r: number = this.radius;
      lower.selfSubXY(r, r);
      upper.selfAddXY(r, r);
    }

    /// @see Shape::ComputeMass
    private static computeMass_s_center = new Vec2();
    private static computeMass_s_s = new Vec2();
    private static computeMass_s_e1 = new Vec2();
    private static computeMass_s_e2 = new Vec2();
    public computeMass(massData: MassData, density: number): void {
      // Polygon mass, centroid, and inertia.
      // Let rho be the polygon density in mass per unit area.
      // Then:
      // mass = rho * int(dA)
      // centroid.x = (1/mass) * rho * int(x * dA)
      // centroid.y = (1/mass) * rho * int(y * dA)
      // I = rho * int((x*x + y*y) * dA)
      //
      // We can compute these integrals by summing all the integrals
      // for each triangle of the polygon. To evaluate the integral
      // for a single triangle, we make a change of variables to
      // the (u,v) coordinates of the triangle:
      // x = x0 + e1x * u + e2x * v
      // y = y0 + e1y * u + e2y * v
      // where 0 <= u && 0 <= v && u + v <= 1.
      //
      // We integrate u from [0,1-v] and then v from [0,1].
      // We also need to use the Jacobian of the transformation:
      // D = cross(e1, e2)
      //
      // Simplification: triangle centroid = (1/3) * (p1 + p2 + p3)
      //
      // The rest of the derivation is handled by computer algebra.

      // DEBUG: Assert(this.count >= 3);

      const center: Vec2 = PolygonShape.computeMass_s_center.setZero();
      let area: number = 0;
      let I: number = 0;

      // Get a reference point for forming triangles.
      // Use the first vertex to reduce round-off errors.
      const s: Vec2 = PolygonShape.computeMass_s_s.copy(this.vertices[0]);

      const k_inv3: number = 1 / 3;

      for (let i: number = 0; i < this.count; ++i) {
        // Triangle vertices.
        const e1: Vec2 = Vec2.SubVV(this.vertices[i], s, PolygonShape.computeMass_s_e1);
        const e2: Vec2 = Vec2.SubVV(this.vertices[(i + 1) % this.count], s, PolygonShape.computeMass_s_e2);

        const D: number = Vec2.CrossVV(e1, e2);

        const triangleArea: number = 0.5 * D;
        area += triangleArea;

        // Area weighted centroid
        center.selfAdd(Vec2.MulSV(triangleArea * k_inv3, Vec2.AddVV(e1, e2, Vec2.s_t0), Vec2.s_t1));

        const ex1: number = e1.x;
        const ey1: number = e1.y;
        const ex2: number = e2.x;
        const ey2: number = e2.y;

        const intx2: number = ex1 * ex1 + ex2 * ex1 + ex2 * ex2;
        const inty2: number = ey1 * ey1 + ey2 * ey1 + ey2 * ey2;

        I += (0.25 * k_inv3 * D) * (intx2 + inty2);
      }

      // Total mass
      massData.mass = density * area;

      // Center of mass
      // DEBUG: Assert(area > epsilon);
      center.selfMul(1 / area);
      Vec2.AddVV(center, s, massData.center);

      // Inertia tensor relative to the local origin (point s).
      massData.I = density * I;

      // Shift to center of mass then to original body origin.
      massData.I += massData.mass * (Vec2.DotVV(massData.center, massData.center) - Vec2.DotVV(center, center));
    }

    private static validate_s_e = new Vec2();
    private static validate_s_v = new Vec2();
    public validate(): boolean {
      for (let i: number = 0; i < this.count; ++i) {
        const i1 = i;
        const i2 = (i + 1) % this.count;
        const p: Vec2 = this.vertices[i1];
        const e: Vec2 = Vec2.SubVV(this.vertices[i2], p, PolygonShape.validate_s_e);

        for (let j: number = 0; j < this.count; ++j) {
          if (j === i1 || j === i2) {
            continue;
          }

          const v: Vec2 = Vec2.SubVV(this.vertices[j], p, PolygonShape.validate_s_v);
          const c: number = Vec2.CrossVV(e, v);
          if (c < 0) {
            return false;
          }
        }
      }

      return true;
    }

    public setupDistanceProxy(proxy: DistanceProxy, index: number): void {
      proxy.vertices = this.vertices;
      proxy.count = this.count;
      proxy.radius = this.radius;
    }

    private static ComputeSubmergedArea_s_normalL = new Vec2();
    private static ComputeSubmergedArea_s_md = new MassData();
    private static ComputeSubmergedArea_s_intoVec = new Vec2();
    private static ComputeSubmergedArea_s_outoVec = new Vec2();
    private static ComputeSubmergedArea_s_center = new Vec2();
    public computeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number {
      // Transform plane into shape co-ordinates
      const normalL: Vec2 = Rot.mulTRV(xf.q, normal, PolygonShape.ComputeSubmergedArea_s_normalL);
      const offsetL: number = offset - Vec2.DotVV(normal, xf.p);

      const depths: number[] = [];
      let diveCount: number = 0;
      let intoIndex: number = -1;
      let outoIndex: number = -1;

      let lastSubmerged: boolean = false;
      for (let i: number = 0; i < this.count; ++i) {
        depths[i] = Vec2.DotVV(normalL, this.vertices[i]) - offsetL;
        const isSubmerged: boolean = depths[i] < (-epsilon);
        if (i > 0) {
          if (isSubmerged) {
            if (!lastSubmerged) {
              intoIndex = i - 1;
              diveCount++;
            }
          } else {
            if (lastSubmerged) {
              outoIndex = i - 1;
              diveCount++;
            }
          }
        }
        lastSubmerged = isSubmerged;
      }
      switch (diveCount) {
        case 0:
          if (lastSubmerged) {
            // Completely submerged
            const md: MassData = PolygonShape.ComputeSubmergedArea_s_md;
            this.computeMass(md, 1);
            Transform.mulXV(xf, md.center, c);
            return md.mass;
          } else {
            // Completely dry
            return 0;
          }
        case 1:
          if (intoIndex === (-1)) {
            intoIndex = this.count - 1;
          } else {
            outoIndex = this.count - 1;
          }
          break;
      }
      const intoIndex2: number = ((intoIndex + 1) % this.count);
      const outoIndex2: number = ((outoIndex + 1) % this.count);
      const intoLamdda: number = (0 - depths[intoIndex]) / (depths[intoIndex2] - depths[intoIndex]);
      const outoLamdda: number = (0 - depths[outoIndex]) / (depths[outoIndex2] - depths[outoIndex]);

      const intoVec: Vec2 = PolygonShape.ComputeSubmergedArea_s_intoVec.set(
        this.vertices[intoIndex].x * (1 - intoLamdda) + this.vertices[intoIndex2].x * intoLamdda,
        this.vertices[intoIndex].y * (1 - intoLamdda) + this.vertices[intoIndex2].y * intoLamdda);
      const outoVec: Vec2 = PolygonShape.ComputeSubmergedArea_s_outoVec.set(
        this.vertices[outoIndex].x * (1 - outoLamdda) + this.vertices[outoIndex2].x * outoLamdda,
        this.vertices[outoIndex].y * (1 - outoLamdda) + this.vertices[outoIndex2].y * outoLamdda);

      // Initialize accumulator
      let area: number = 0;
      const center: Vec2 = PolygonShape.ComputeSubmergedArea_s_center.setZero();
      let p2: Vec2 = this.vertices[intoIndex2];
      let p3: Vec2;

      // An awkward loop from intoIndex2+1 to outIndex2
      let i: number = intoIndex2;
      while (i !== outoIndex2) {
        i = (i + 1) % this.count;
        if (i === outoIndex2) {
          p3 = outoVec;
        } else {
          p3 = this.vertices[i];
        }

        const triangleArea: number = 0.5 * ((p2.x - intoVec.x) * (p3.y - intoVec.y) - (p2.y - intoVec.y) * (p3.x - intoVec.x));
        area += triangleArea;
        // Area weighted centroid
        center.x += triangleArea * (intoVec.x + p2.x + p3.x) / 3;
        center.y += triangleArea * (intoVec.y + p2.y + p3.y) / 3;

        p2 = p3;
      }

      // Normalize and transform centroid
      center.selfMul(1 / area);
      Transform.mulXV(xf, center, c);

      return area;
    }

    public dump(log: (format: string, ...args: any[]) => void): void {
      log("    const shape: PolygonShape = new PolygonShape();\n");
      log("    const vs: Vec2[] = [];\n");
      for (let i: number = 0; i < this.count; ++i) {
        log("    vs[%d] = new Vec2(%.15f, %.15f);\n", i, this.vertices[i].x, this.vertices[i].y);
      }
      log("    shape.Set(vs, %d);\n", this.count);
    }

    private static computeCentroid_s_s = new Vec2();
    private static computeCentroid_s_p1 = new Vec2();
    private static computeCentroid_s_p2 = new Vec2();
    private static computeCentroid_s_p3 = new Vec2();
    private static computeCentroid_s_e1 = new Vec2();
    private static computeCentroid_s_e2 = new Vec2();
    public static computeCentroid(vs: Vec2[], count: number, out: Vec2): Vec2 {
      // DEBUG: Assert(count >= 3);

      const c: Vec2 = out; c.setZero();
      let area: number = 0;

      // Get a reference point for forming triangles.
      // Use the first vertex to reduce round-off errors.
      const s: Vec2 = PolygonShape.computeCentroid_s_s.copy(vs[0]);

      const inv3: number = 1 / 3;

      for (let i: number = 0; i < count; ++i) {
        // Triangle vertices.
        const p1: Vec2 = Vec2.SubVV(vs[0], s, PolygonShape.computeCentroid_s_p1);
        const p2: Vec2 = Vec2.SubVV(vs[i], s, PolygonShape.computeCentroid_s_p2);
        const p3: Vec2 = Vec2.SubVV(vs[(i + 1) % count], s, PolygonShape.computeCentroid_s_p3);

        const e1: Vec2 = Vec2.SubVV(p2, p1, PolygonShape.computeCentroid_s_e1);
        const e2: Vec2 = Vec2.SubVV(p3, p1, PolygonShape.computeCentroid_s_e2);

        const D: number = Vec2.CrossVV(e1, e2);

        const triangleArea: number = 0.5 * D;
        area += triangleArea;

        // Area weighted centroid
        c.x += triangleArea * inv3 * (p1.x + p2.x + p3.x);
        c.y += triangleArea * inv3 * (p1.y + p2.y + p3.y);
      }

      // Centroid
      // DEBUG: Assert(area > epsilon);
      // c = (1.0f / area) * c + s;
      c.x = (1 / area) * c.x + s.x;
      c.y = (1 / area) * c.y + s.y;
      return c;
    }

    /*
    public static ComputeOBB(obb, vs, count) {
      const i: number = 0;
      const p: Array = [count + 1];
      for (i = 0; i < count; ++i) {
        p[i] = vs[i];
      }
      p[count] = p[0];
      const minArea = maxFloat;
      for (i = 1; i <= count; ++i) {
        const root = p[i - 1];
        const uxX = p[i].x - root.x;
        const uxY = p[i].y - root.y;
        const length = Sqrt(uxX * uxX + uxY * uxY);
        uxX /= length;
        uxY /= length;
        const uyX = (-uxY);
        const uyY = uxX;
        const lowerX = maxFloat;
        const lowerY = maxFloat;
        const upperX = (-maxFloat);
        const upperY = (-maxFloat);
        for (let j: number = 0; j < count; ++j) {
          const dX = p[j].x - root.x;
          const dY = p[j].y - root.y;
          const rX = (uxX * dX + uxY * dY);
          const rY = (uyX * dX + uyY * dY);
          if (rX < lowerX) lowerX = rX;
          if (rY < lowerY) lowerY = rY;
          if (rX > upperX) upperX = rX;
          if (rY > upperY) upperY = rY;
        }
        const area = (upperX - lowerX) * (upperY - lowerY);
        if (area < 0.95 * minArea) {
          minArea = area;
          obb.R.ex.x = uxX;
          obb.R.ex.y = uxY;
          obb.R.ey.x = uyX;
          obb.R.ey.y = uyY;
          const center_x: number = 0.5 * (lowerX + upperX);
          const center_y: number = 0.5 * (lowerY + upperY);
          const tMat = obb.R;
          obb.center.x = root.x + (tMat.ex.x * center_x + tMat.ey.x * center_y);
          obb.center.y = root.y + (tMat.ex.y * center_x + tMat.ey.y * center_y);
          obb.extents.x = 0.5 * (upperX - lowerX);
          obb.extents.y = 0.5 * (upperY - lowerY);
        }
      }
    }
    */
  }
}

