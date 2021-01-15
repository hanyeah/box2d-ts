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

/// A distance proxy is used by the GJK algorithm.
/// It encapsulates any shape.
  export class DistanceProxy {
    public readonly buffer: Vec2[] = Vec2.MakeArray(2);
    public vertices: Vec2[] = this.buffer;
    public count: number = 0;
    public radius: number = 0;

    public copy(other: DistanceProxy): this {
      if (other.vertices === other.buffer) {
        this.vertices = this.buffer;
        this.buffer[0].copy(other.buffer[0]);
        this.buffer[1].copy(other.buffer[1]);
      } else {
        this.vertices = other.vertices;
      }
      this.count = other.count;
      this.radius = other.radius;
      return this;
    }

    public reset(): DistanceProxy {
      this.vertices = this.buffer;
      this.count = 0;
      this.radius = 0;
      return this;
    }

    public setShape(shape: Shape, index: number): void {
      shape.setupDistanceProxy(this, index);
    }

    public setVerticesRadius(vertices: Vec2[], count: number, radius: number): void {
      this.vertices = vertices;
      this.count = count;
      this.radius = radius;
    }

    public getSupport(d: Vec2): number {
      let bestIndex: number = 0;
      let bestValue: number = Vec2.DotVV(this.vertices[0], d);
      for (let i: number = 1; i < this.count; ++i) {
        const value: number = Vec2.DotVV(this.vertices[i], d);
        if (value > bestValue) {
          bestIndex = i;
          bestValue = value;
        }
      }

      return bestIndex;
    }

    public getSupportVertex(d: Vec2): Vec2 {
      let bestIndex: number = 0;
      let bestValue: number = Vec2.DotVV(this.vertices[0], d);
      for (let i: number = 1; i < this.count; ++i) {
        const value: number = Vec2.DotVV(this.vertices[i], d);
        if (value > bestValue) {
          bestIndex = i;
          bestValue = value;
        }
      }

      return this.vertices[bestIndex];
    }

    public getVertexCount(): number {
      return this.count;
    }

    public getVertex(index: number): Vec2 {
      // DEBUG: Assert(0 <= index && index < this.count);
      return this.vertices[index];
    }
  }

  export class SimplexCache {
    public metric: number = 0;
    public count: number = 0;
    public readonly indexA: [number, number, number] = [0, 0, 0];
    public readonly indexB: [number, number, number] = [0, 0, 0];

    public reset(): SimplexCache {
      this.metric = 0;
      this.count = 0;
      return this;
    }
  }

  export class DistanceInput {
    public readonly proxyA: DistanceProxy = new DistanceProxy();
    public readonly proxyB: DistanceProxy = new DistanceProxy();
    public readonly transformA: Transform = new Transform();
    public readonly transformB: Transform = new Transform();
    public useRadii: boolean = false;

    public reset(): DistanceInput {
      this.proxyA.reset();
      this.proxyB.reset();
      this.transformA.setIdentity();
      this.transformB.setIdentity();
      this.useRadii = false;
      return this;
    }
  }

  export class DistanceOutput {
    public readonly pointA: Vec2 = new Vec2();
    public readonly pointB: Vec2 = new Vec2();
    public distance: number = 0;
    public iterations: number = 0; ///< number of GJK iterations used

    public reset(): DistanceOutput {
      this.pointA.setZero();
      this.pointB.setZero();
      this.distance = 0;
      this.iterations = 0;
      return this;
    }
  }

/// Input parameters for ShapeCast
  export class ShapeCastInput {
    public readonly proxyA: DistanceProxy = new DistanceProxy();
    public readonly proxyB: DistanceProxy = new DistanceProxy();
    public readonly transformA: Transform = new Transform();
    public readonly transformB: Transform = new Transform();
    public readonly translationB: Vec2 = new Vec2();
  }

/// Output results for ShapeCast
  export class ShapeCastOutput {
    public readonly point: Vec2 = new Vec2();
    public readonly normal: Vec2 = new Vec2();
    public lambda: number = 0.0;
    public iterations: number = 0;
  }

  export let gjkCalls: number = 0;
  export let gjkIters: number = 0;
  export let gjkMaxIters: number = 0;
  export function gjkReset(): void {
    gjkCalls = 0;
    gjkIters = 0;
    gjkMaxIters = 0;
  }

  export class SimplexVertex {
    public readonly wA: Vec2 = new Vec2(); // support point in proxyA
    public readonly wB: Vec2 = new Vec2(); // support point in proxyB
    public readonly w: Vec2 = new Vec2(); // wB - wA
    public a: number = 0; // barycentric coordinate for closest point
    public indexA: number = 0; // wA index
    public indexB: number = 0; // wB index

    public copy(other: SimplexVertex): SimplexVertex {
      this.wA.copy(other.wA);     // support point in proxyA
      this.wB.copy(other.wB);     // support point in proxyB
      this.w.copy(other.w);       // wB - wA
      this.a = other.a;           // barycentric coordinate for closest point
      this.indexA = other.indexA; // wA index
      this.indexB = other.indexB; // wB index
      return this;
    }
  }

  export class Simplex {
    public readonly v1: SimplexVertex = new SimplexVertex();
    public readonly v2: SimplexVertex = new SimplexVertex();
    public readonly v3: SimplexVertex = new SimplexVertex();
    public readonly vertices: SimplexVertex[] = [/*3*/];
    public count: number = 0;

    constructor() {
      this.vertices[0] = this.v1;
      this.vertices[1] = this.v2;
      this.vertices[2] = this.v3;
    }

    public readCache(cache: SimplexCache, proxyA: DistanceProxy, transformA: Transform, proxyB: DistanceProxy, transformB: Transform): void {
      // DEBUG: Assert(0 <= cache.count && cache.count <= 3);

      // Copy data from cache.
      this.count = cache.count;
      const vertices: SimplexVertex[] = this.vertices;
      for (let i: number = 0; i < this.count; ++i) {
        const v: SimplexVertex = vertices[i];
        v.indexA = cache.indexA[i];
        v.indexB = cache.indexB[i];
        const wALocal: Vec2 = proxyA.getVertex(v.indexA);
        const wBLocal: Vec2 = proxyB.getVertex(v.indexB);
        Transform.mulXV(transformA, wALocal, v.wA);
        Transform.mulXV(transformB, wBLocal, v.wB);
        Vec2.SubVV(v.wB, v.wA, v.w);
        v.a = 0;
      }

      // Compute the new simplex metric, if it is substantially different than
      // old metric then flush the simplex.
      if (this.count > 1) {
        const metric1: number = cache.metric;
        const metric2: number = this.getMetric();
        if (metric2 < 0.5 * metric1 || 2 * metric1 < metric2 || metric2 < epsilon) {
          // Reset the simplex.
          this.count = 0;
        }
      }

      // If the cache is empty or invalid ...
      if (this.count === 0) {
        const v: SimplexVertex = vertices[0];
        v.indexA = 0;
        v.indexB = 0;
        const wALocal: Vec2 = proxyA.getVertex(0);
        const wBLocal: Vec2 = proxyB.getVertex(0);
        Transform.mulXV(transformA, wALocal, v.wA);
        Transform.mulXV(transformB, wBLocal, v.wB);
        Vec2.SubVV(v.wB, v.wA, v.w);
        v.a = 1;
        this.count = 1;
      }
    }

    public writeCache(cache: SimplexCache): void {
      cache.metric = this.getMetric();
      cache.count = this.count;
      const vertices: SimplexVertex[] = this.vertices;
      for (let i: number = 0; i < this.count; ++i) {
        cache.indexA[i] = vertices[i].indexA;
        cache.indexB[i] = vertices[i].indexB;
      }
    }

    public getSearchDirection(out: Vec2): Vec2 {
      switch (this.count) {
        case 1:
          return Vec2.NegV(this.v1.w, out);

        case 2: {
          const e12: Vec2 = Vec2.SubVV(this.v2.w, this.v1.w, out);
          const sgn: number = Vec2.CrossVV(e12, Vec2.NegV(this.v1.w, Vec2.s_t0));
          if (sgn > 0) {
            // Origin is left of e12.
            return Vec2.CrossOneV(e12, out);
          } else {
            // Origin is right of e12.
            return Vec2.CrossVOne(e12, out);
          }
        }

        default:
          // DEBUG: Assert(false);
          return out.setZero();
      }
    }

    public getClosestPoint(out: Vec2): Vec2 {
      switch (this.count) {
        case 0:
          // DEBUG: Assert(false);
          return out.setZero();

        case 1:
          return out.copy(this.v1.w);

        case 2:
          return out.set(
            this.v1.a * this.v1.w.x + this.v2.a * this.v2.w.x,
            this.v1.a * this.v1.w.y + this.v2.a * this.v2.w.y);

        case 3:
          return out.setZero();

        default:
          // DEBUG: Assert(false);
          return out.setZero();
      }
    }

    public getWitnessPoints(pA: Vec2, pB: Vec2): void {
      switch (this.count) {
        case 0:
          // DEBUG: Assert(false);
          break;

        case 1:
          pA.copy(this.v1.wA);
          pB.copy(this.v1.wB);
          break;

        case 2:
          pA.x = this.v1.a * this.v1.wA.x + this.v2.a * this.v2.wA.x;
          pA.y = this.v1.a * this.v1.wA.y + this.v2.a * this.v2.wA.y;
          pB.x = this.v1.a * this.v1.wB.x + this.v2.a * this.v2.wB.x;
          pB.y = this.v1.a * this.v1.wB.y + this.v2.a * this.v2.wB.y;
          break;

        case 3:
          pB.x = pA.x = this.v1.a * this.v1.wA.x + this.v2.a * this.v2.wA.x + this.v3.a * this.v3.wA.x;
          pB.y = pA.y = this.v1.a * this.v1.wA.y + this.v2.a * this.v2.wA.y + this.v3.a * this.v3.wA.y;
          break;

        default:
          // DEBUG: Assert(false);
          break;
      }
    }

    public getMetric(): number {
      switch (this.count) {
        case 0:
          // DEBUG: Assert(false);
          return 0;

        case 1:
          return 0;

        case 2:
          return Vec2.DistanceVV(this.v1.w, this.v2.w);

        case 3:
          return Vec2.CrossVV(Vec2.SubVV(this.v2.w, this.v1.w, Vec2.s_t0), Vec2.SubVV(this.v3.w, this.v1.w, Vec2.s_t1));

        default:
          // DEBUG: Assert(false);
          return 0;
      }
    }

    public solve2(): void {
      const w1: Vec2 = this.v1.w;
      const w2: Vec2 = this.v2.w;
      const e12: Vec2 = Vec2.SubVV(w2, w1, Simplex.s_e12);

      // w1 region
      const d12_2: number = (-Vec2.DotVV(w1, e12));
      if (d12_2 <= 0) {
        // a2 <= 0, so we clamp it to 0
        this.v1.a = 1;
        this.count = 1;
        return;
      }

      // w2 region
      const d12_1: number = Vec2.DotVV(w2, e12);
      if (d12_1 <= 0) {
        // a1 <= 0, so we clamp it to 0
        this.v2.a = 1;
        this.count = 1;
        this.v1.copy(this.v2);
        return;
      }

      // Must be in e12 region.
      const inv_d12: number = 1 / (d12_1 + d12_2);
      this.v1.a = d12_1 * inv_d12;
      this.v2.a = d12_2 * inv_d12;
      this.count = 2;
    }

    public solve3(): void {
      const w1: Vec2 = this.v1.w;
      const w2: Vec2 = this.v2.w;
      const w3: Vec2 = this.v3.w;

      // Edge12
      // [1      1     ][a1] = [1]
      // [w1.e12 w2.e12][a2] = [0]
      // a3 = 0
      const e12: Vec2 = Vec2.SubVV(w2, w1, Simplex.s_e12);
      const w1e12: number = Vec2.DotVV(w1, e12);
      const w2e12: number = Vec2.DotVV(w2, e12);
      const d12_1: number = w2e12;
      const d12_2: number = (-w1e12);

      // Edge13
      // [1      1     ][a1] = [1]
      // [w1.e13 w3.e13][a3] = [0]
      // a2 = 0
      const e13: Vec2 = Vec2.SubVV(w3, w1, Simplex.s_e13);
      const w1e13: number = Vec2.DotVV(w1, e13);
      const w3e13: number = Vec2.DotVV(w3, e13);
      const d13_1: number = w3e13;
      const d13_2: number = (-w1e13);

      // Edge23
      // [1      1     ][a2] = [1]
      // [w2.e23 w3.e23][a3] = [0]
      // a1 = 0
      const e23: Vec2 = Vec2.SubVV(w3, w2, Simplex.s_e23);
      const w2e23: number = Vec2.DotVV(w2, e23);
      const w3e23: number = Vec2.DotVV(w3, e23);
      const d23_1: number = w3e23;
      const d23_2: number = (-w2e23);

      // Triangle123
      const n123: number = Vec2.CrossVV(e12, e13);

      const d123_1: number = n123 * Vec2.CrossVV(w2, w3);
      const d123_2: number = n123 * Vec2.CrossVV(w3, w1);
      const d123_3: number = n123 * Vec2.CrossVV(w1, w2);

      // w1 region
      if (d12_2 <= 0 && d13_2 <= 0) {
        this.v1.a = 1;
        this.count = 1;
        return;
      }

      // e12
      if (d12_1 > 0 && d12_2 > 0 && d123_3 <= 0) {
        const inv_d12: number = 1 / (d12_1 + d12_2);
        this.v1.a = d12_1 * inv_d12;
        this.v2.a = d12_2 * inv_d12;
        this.count = 2;
        return;
      }

      // e13
      if (d13_1 > 0 && d13_2 > 0 && d123_2 <= 0) {
        const inv_d13: number = 1 / (d13_1 + d13_2);
        this.v1.a = d13_1 * inv_d13;
        this.v3.a = d13_2 * inv_d13;
        this.count = 2;
        this.v2.copy(this.v3);
        return;
      }

      // w2 region
      if (d12_1 <= 0 && d23_2 <= 0) {
        this.v2.a = 1;
        this.count = 1;
        this.v1.copy(this.v2);
        return;
      }

      // w3 region
      if (d13_1 <= 0 && d23_1 <= 0) {
        this.v3.a = 1;
        this.count = 1;
        this.v1.copy(this.v3);
        return;
      }

      // e23
      if (d23_1 > 0 && d23_2 > 0 && d123_1 <= 0) {
        const inv_d23: number = 1 / (d23_1 + d23_2);
        this.v2.a = d23_1 * inv_d23;
        this.v3.a = d23_2 * inv_d23;
        this.count = 2;
        this.v1.copy(this.v3);
        return;
      }

      // Must be in triangle123
      const inv_d123: number = 1 / (d123_1 + d123_2 + d123_3);
      this.v1.a = d123_1 * inv_d123;
      this.v2.a = d123_2 * inv_d123;
      this.v3.a = d123_3 * inv_d123;
      this.count = 3;
    }
    private static s_e12: Vec2 = new Vec2();
    private static s_e13: Vec2 = new Vec2();
    private static s_e23: Vec2 = new Vec2();
  }

  const distance_s_simplex: Simplex = new Simplex();
  const distance_s_saveA: [number, number, number] = [0, 0, 0];
  const distance_s_saveB: [number, number, number] = [0, 0, 0];
  const distance_s_p: Vec2 = new Vec2();
  const distance_s_d: Vec2 = new Vec2();
  const distance_s_normal: Vec2 = new Vec2();
  const distance_s_supportA: Vec2 = new Vec2();
  const distance_s_supportB: Vec2 = new Vec2();
  export function distance(output: DistanceOutput, cache: SimplexCache, input: DistanceInput): void {
    ++gjkCalls;

    const proxyA: DistanceProxy = input.proxyA;
    const proxyB: DistanceProxy = input.proxyB;

    const transformA: Transform = input.transformA;
    const transformB: Transform = input.transformB;

    // Initialize the simplex.
    const simplex: Simplex = distance_s_simplex;
    simplex.readCache(cache, proxyA, transformA, proxyB, transformB);

    // Get simplex vertices as an array.
    const vertices: SimplexVertex[] = simplex.vertices;
    const k_maxIters: number = 20;

    // These store the vertices of the last simplex so that we
    // can check for duplicates and prevent cycling.
    const saveA: [number, number, number] = distance_s_saveA;
    const saveB: [number, number, number] = distance_s_saveB;
    let saveCount: number = 0;

    // Main iteration loop.
    let iter: number = 0;
    while (iter < k_maxIters) {
      // Copy simplex so we can identify duplicates.
      saveCount = simplex.count;
      for (let i: number = 0; i < saveCount; ++i) {
        saveA[i] = vertices[i].indexA;
        saveB[i] = vertices[i].indexB;
      }

      switch (simplex.count) {
        case 1:
          break;

        case 2:
          simplex.solve2();
          break;

        case 3:
          simplex.solve3();
          break;

        default:
          // DEBUG: Assert(false);
          break;
      }

      // If we have 3 points, then the origin is in the corresponding triangle.
      if (simplex.count === 3) {
        break;
      }

      // Get search direction.
      const d: Vec2 = simplex.getSearchDirection(distance_s_d);

      // Ensure the search direction is numerically fit.
      if (d.lengthSquared() < epsilonSq) {
        // The origin is probably contained by a line segment
        // or triangle. Thus the shapes are overlapped.

        // We can't return zero here even though there may be overlap.
        // In case the simplex is a point, segment, or triangle it is difficult
        // to determine if the origin is contained in the CSO or very close to it.
        break;
      }

      // Compute a tentative new simplex vertex using support points.
      const vertex: SimplexVertex = vertices[simplex.count];
      vertex.indexA = proxyA.getSupport(Rot.mulTRV(transformA.q, Vec2.NegV(d, Vec2.s_t0), distance_s_supportA));
      Transform.mulXV(transformA, proxyA.getVertex(vertex.indexA), vertex.wA);
      vertex.indexB = proxyB.getSupport(Rot.mulTRV(transformB.q, d, distance_s_supportB));
      Transform.mulXV(transformB, proxyB.getVertex(vertex.indexB), vertex.wB);
      Vec2.SubVV(vertex.wB, vertex.wA, vertex.w);

      // Iteration count is equated to the number of support point calls.
      ++iter;
      ++gjkIters;

      // Check for duplicate support points. This is the main termination criteria.
      let duplicate: boolean = false;
      for (let i: number = 0; i < saveCount; ++i) {
        if (vertex.indexA === saveA[i] && vertex.indexB === saveB[i]) {
          duplicate = true;
          break;
        }
      }

      // If we found a duplicate support point we must exit to avoid cycling.
      if (duplicate) {
        break;
      }

      // New vertex is ok and needed.
      ++simplex.count;
    }

    gjkMaxIters = Max(gjkMaxIters, iter);

    // Prepare output.
    simplex.getWitnessPoints(output.pointA, output.pointB);
    output.distance = Vec2.DistanceVV(output.pointA, output.pointB);
    output.iterations = iter;

    // Cache the simplex.
    simplex.writeCache(cache);

    // Apply radii if requested.
    if (input.useRadii) {
      const rA: number = proxyA.radius;
      const rB: number = proxyB.radius;

      if (output.distance > (rA + rB) && output.distance > epsilon) {
        // Shapes are still no overlapped.
        // Move the witness points to the outer surface.
        output.distance -= rA + rB;
        const normal: Vec2 = Vec2.SubVV(output.pointB, output.pointA, distance_s_normal);
        normal.normalize();
        output.pointA.selfMulAdd(rA, normal);
        output.pointB.selfMulSub(rB, normal);
      } else {
        // Shapes are overlapped when radii are considered.
        // Move the witness points to the middle.
        const p: Vec2 = Vec2.MidVV(output.pointA, output.pointB, distance_s_p);
        output.pointA.copy(p);
        output.pointB.copy(p);
        output.distance = 0;
      }
    }
  }

/// Perform a linear shape cast of shape B moving and shape A fixed. Determines the hit point, normal, and translation fraction.

// GJK-raycast
// Algorithm by Gino van den Bergen.
// "Smooth Mesh Contacts with GJK" in Game Physics Pearls. 2010
// bool ShapeCast(ShapeCastOutput* output, const ShapeCastInput* input);
  const shapeCast_s_n = new Vec2();
  const shapeCast_s_simplex = new Simplex();
  const shapeCast_s_wA = new Vec2();
  const shapeCast_s_wB = new Vec2();
  const shapeCast_s_v = new Vec2();
  const shapeCast_s_p = new Vec2();
  const shapeCast_s_pointA = new Vec2();
  const shapeCast_s_pointB = new Vec2();
  export function shapeCast(output: ShapeCastOutput, input: ShapeCastInput): boolean {
    output.iterations = 0;
    output.lambda = 1.0;
    output.normal.setZero();
    output.point.setZero();

    // const DistanceProxy* proxyA = &input.proxyA;
    const proxyA = input.proxyA;
    // const DistanceProxy* proxyB = &input.proxyB;
    const proxyB = input.proxyB;

    // float32 radiusA = Max(proxyA.radius, polygonRadius);
    const radiusA = Max(proxyA.radius, polygonRadius);
    // float32 radiusB = Max(proxyB.radius, polygonRadius);
    const radiusB = Max(proxyB.radius, polygonRadius);
    // float32 radius = radiusA + radiusB;
    const radius = radiusA + radiusB;

    // Transform xfA = input.transformA;
    const xfA = input.transformA;
    // Transform xfB = input.transformB;
    const xfB = input.transformB;

    // Vec2 r = input.translationB;
    const r = input.translationB;
    // Vec2 n(0.0f, 0.0f);
    const n = shapeCast_s_n.set(0.0, 0.0);
    // float32 lambda = 0.0f;
    let lambda = 0.0;

    // Initial simplex
    const simplex = shapeCast_s_simplex;
    simplex.count = 0;

    // Get simplex vertices as an array.
    // SimplexVertex* vertices = &simplex.v1;
    const vertices = simplex.vertices;

    // Get support point in -r direction
    // int32 indexA = proxyA.GetSupport(MulT(xfA.q, -r));
    let indexA = proxyA.getSupport(Rot.mulTRV(xfA.q, Vec2.NegV(r, Vec2.s_t1), Vec2.s_t0));
    // Vec2 wA = Mul(xfA, proxyA.GetVertex(indexA));
    let wA = Transform.mulXV(xfA, proxyA.getVertex(indexA), shapeCast_s_wA);
    // int32 indexB = proxyB.GetSupport(MulT(xfB.q, r));
    let indexB = proxyB.getSupport(Rot.mulTRV(xfB.q, r, Vec2.s_t0));
    // Vec2 wB = Mul(xfB, proxyB.GetVertex(indexB));
    let wB = Transform.mulXV(xfB, proxyB.getVertex(indexB), shapeCast_s_wB);
    // Vec2 v = wA - wB;
    const v = Vec2.SubVV(wA, wB, shapeCast_s_v);

    // Sigma is the target distance between polygons
    // float32 sigma = Max(polygonRadius, radius - polygonRadius);
    const sigma = Max(polygonRadius, radius - polygonRadius);
    // const float32 tolerance = 0.5f * linearSlop;
    const tolerance = 0.5 * linearSlop;

    // Main iteration loop.
    // const int32 k_maxIters = 20;
    const k_maxIters = 20;
    // int32 iter = 0;
    let iter = 0;
    // while (iter < k_maxIters && v.Length() - sigma > tolerance)
    while (iter < k_maxIters && v.length() - sigma > tolerance) {
      // DEBUG: Assert(simplex.count < 3);

      output.iterations += 1;

      // Support in direction -v (A - B)
      // indexA = proxyA.GetSupport(MulT(xfA.q, -v));
      indexA = proxyA.getSupport(Rot.mulTRV(xfA.q, Vec2.NegV(v, Vec2.s_t1), Vec2.s_t0));
      // wA = Mul(xfA, proxyA.GetVertex(indexA));
      wA = Transform.mulXV(xfA, proxyA.getVertex(indexA), shapeCast_s_wA);
      // indexB = proxyB.GetSupport(MulT(xfB.q, v));
      indexB = proxyB.getSupport(Rot.mulTRV(xfB.q, v, Vec2.s_t0));
      // wB = Mul(xfB, proxyB.GetVertex(indexB));
      wB = Transform.mulXV(xfB, proxyB.getVertex(indexB), shapeCast_s_wB);
      // Vec2 p = wA - wB;
      const p = Vec2.SubVV(wA, wB, shapeCast_s_p);

      // -v is a normal at p
      v.normalize();

      // Intersect ray with plane
      const vp = Vec2.DotVV(v, p);
      const vr = Vec2.DotVV(v, r);
      if (vp - sigma > lambda * vr) {
        if (vr <= 0.0) {
          return false;
        }

        lambda = (vp - sigma) / vr;
        if (lambda > 1.0) {
          return false;
        }

        // n = -v;
        n.copy(v).selfNeg();
        simplex.count = 0;
      }

      // Reverse simplex since it works with B - A.
      // Shift by lambda * r because we want the closest point to the current clip point.
      // Note that the support point p is not shifted because we want the plane equation
      // to be formed in unshifted space.
      // SimplexVertex* vertex = vertices + simplex.count;
      const vertex: SimplexVertex = vertices[simplex.count];
      vertex.indexA = indexB;
      // vertex.wA = wB + lambda * r;
      vertex.wA.copy(wB).selfMulAdd(lambda, r);
      vertex.indexB = indexA;
      // vertex.wB = wA;
      vertex.wB.copy(wA);
      // vertex.w = vertex.wB - vertex.wA;
      vertex.w.copy(vertex.wB).selfSub(vertex.wA);
      vertex.a = 1.0;
      simplex.count += 1;

      switch (simplex.count) {
        case 1:
          break;

        case 2:
          simplex.solve2();
          break;

        case 3:
          simplex.solve3();
          break;

        default:
        // DEBUG: Assert(false);
      }

      // If we have 3 points, then the origin is in the corresponding triangle.
      if (simplex.count === 3) {
        // Overlap
        return false;
      }

      // Get search direction.
      // v = simplex.GetClosestPoint();
      simplex.getClosestPoint(v);

      // Iteration count is equated to the number of support point calls.
      ++iter;
    }

    if (iter === 0) {
      // Initial overlap
      return false;
    }

    // Prepare output.
    const pointA = shapeCast_s_pointA;
    const pointB = shapeCast_s_pointB;
    simplex.getWitnessPoints(pointA, pointB);

    if (v.lengthSquared() > 0.0) {
      // n = -v;
      n.copy(v).selfNeg();
      n.normalize();
    }

    // output.point = pointA + radiusA * n;
    output.normal.copy(n);
    output.lambda = lambda;
    output.iterations = iter;
    return true;
  }

}
