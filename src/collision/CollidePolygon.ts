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

// Find the max separation between poly1 and poly2 using edge normals from poly1.

namespace b2 {
  const FindMaxSeparation_s_xf: Transform = new Transform();
  const FindMaxSeparation_s_n: Vec2 = new Vec2();
  const FindMaxSeparation_s_v1: Vec2 = new Vec2();
  function FindMaxSeparation(edgeIndex: [number], poly1: PolygonShape, xf1: Transform, poly2: PolygonShape, xf2: Transform): number {
    const count1: number = poly1.count;
    const count2: number = poly2.count;
    const n1s: Vec2[] = poly1.normals;
    const v1s: Vec2[] = poly1.vertices;
    const v2s: Vec2[] = poly2.vertices;
    const xf: Transform = Transform.MulTXX(xf2, xf1, FindMaxSeparation_s_xf);

    let bestIndex: number = 0;
    let maxSeparation: number = -maxFloat;

    for (let i: number = 0; i < count1; ++i) {
      // Get poly1 normal in frame2.
      const n: Vec2 = Rot.MulRV(xf.q, n1s[i], FindMaxSeparation_s_n);
      const v1: Vec2 = Transform.MulXV(xf, v1s[i], FindMaxSeparation_s_v1);

      // Find deepest point for normal i.
      let si: number = maxFloat;
      for (let j: number = 0; j < count2; ++j) {
        const sij = Vec2.DotVV(n, Vec2.SubVV(v2s[j], v1, Vec2.s_t0));
        if (sij < si) {
          si = sij;
        }
      }

      if (si > maxSeparation) {
        maxSeparation = si;
        bestIndex = i;
      }
    }

    edgeIndex[0] = bestIndex;
    return maxSeparation;
  }

  const FindIncidentEdge_s_normal1: Vec2 = new Vec2();
  function FindIncidentEdge(c: [ClipVertex, ClipVertex], poly1: PolygonShape, xf1: Transform, edge1: number, poly2: PolygonShape, xf2: Transform): void {
    const normals1: Vec2[] = poly1.normals;

    const count2: number = poly2.count;
    const vertices2: Vec2[] = poly2.vertices;
    const normals2: Vec2[] = poly2.normals;

    // DEBUG: Assert(0 <= edge1 && edge1 < poly1.count);

    // Get the normal of the reference edge in poly2's frame.
    const normal1: Vec2 = Rot.MulTRV(xf2.q, Rot.MulRV(xf1.q, normals1[edge1], Vec2.s_t0), FindIncidentEdge_s_normal1);

    // Find the incident edge on poly2.
    let index: number = 0;
    let minDot: number = maxFloat;
    for (let i: number = 0; i < count2; ++i) {
      const dot: number = Vec2.DotVV(normal1, normals2[i]);
      if (dot < minDot) {
        minDot = dot;
        index = i;
      }
    }

    // Build the clip vertices for the incident edge.
    const i1: number = index;
    const i2: number = i1 + 1 < count2 ? i1 + 1 : 0;

    const c0: ClipVertex = c[0];
    Transform.MulXV(xf2, vertices2[i1], c0.v);
    const cf0: ContactFeature = c0.id.cf;
    cf0.indexA = edge1;
    cf0.indexB = i1;
    cf0.typeA = ContactFeatureType.e_face;
    cf0.typeB = ContactFeatureType.e_vertex;

    const c1: ClipVertex = c[1];
    Transform.MulXV(xf2, vertices2[i2], c1.v);
    const cf1: ContactFeature = c1.id.cf;
    cf1.indexA = edge1;
    cf1.indexB = i2;
    cf1.typeA = ContactFeatureType.e_face;
    cf1.typeB = ContactFeatureType.e_vertex;
  }

// Find edge normal of max separation on A - return if separating axis is found
// Find edge normal of max separation on B - return if separation axis is found
// Choose reference edge as min(minA, minB)
// Find incident edge
// Clip

// The normal points from 1 to 2
  const CollidePolygons_s_incidentEdge: [ClipVertex, ClipVertex] = [ new ClipVertex(), new ClipVertex() ];
  const CollidePolygons_s_clipPoints1: [ClipVertex, ClipVertex] = [ new ClipVertex(), new ClipVertex() ];
  const CollidePolygons_s_clipPoints2: [ClipVertex, ClipVertex] = [ new ClipVertex(), new ClipVertex() ];
  const CollidePolygons_s_edgeA: [number] = [ 0 ];
  const CollidePolygons_s_edgeB: [number] = [ 0 ];
  const CollidePolygons_s_localTangent: Vec2 = new Vec2();
  const CollidePolygons_s_localNormal: Vec2 = new Vec2();
  const CollidePolygons_s_planePoint: Vec2 = new Vec2();
  const CollidePolygons_s_normal: Vec2 = new Vec2();
  const CollidePolygons_s_tangent: Vec2 = new Vec2();
  const CollidePolygons_s_ntangent: Vec2 = new Vec2();
  const CollidePolygons_s_v11: Vec2 = new Vec2();
  const CollidePolygons_s_v12: Vec2 = new Vec2();
  export function CollidePolygons(manifold: Manifold, polyA: PolygonShape, xfA: Transform, polyB: PolygonShape, xfB: Transform): void {
    manifold.pointCount = 0;
    const totalRadius: number = polyA.radius + polyB.radius;

    const edgeA: [number] = CollidePolygons_s_edgeA; edgeA[0] = 0;
    const separationA: number = FindMaxSeparation(edgeA, polyA, xfA, polyB, xfB);
    if (separationA > totalRadius) {
      return;
    }

    const edgeB: [number] = CollidePolygons_s_edgeB; edgeB[0] = 0;
    const separationB: number = FindMaxSeparation(edgeB, polyB, xfB, polyA, xfA);
    if (separationB > totalRadius) {
      return;
    }

    let poly1: PolygonShape; // reference polygon
    let poly2: PolygonShape; // incident polygon
    let xf1: Transform, xf2: Transform;
    let edge1: number = 0; // reference edge
    let flip: number = 0;
    const k_tol: number = 0.1 * linearSlop;

    if (separationB > separationA + k_tol) {
      poly1 = polyB;
      poly2 = polyA;
      xf1 = xfB;
      xf2 = xfA;
      edge1 = edgeB[0];
      manifold.type = ManifoldType.e_faceB;
      flip = 1;
    } else {
      poly1 = polyA;
      poly2 = polyB;
      xf1 = xfA;
      xf2 = xfB;
      edge1 = edgeA[0];
      manifold.type = ManifoldType.e_faceA;
      flip = 0;
    }

    const incidentEdge: [ClipVertex, ClipVertex] = CollidePolygons_s_incidentEdge;
    FindIncidentEdge(incidentEdge, poly1, xf1, edge1, poly2, xf2);

    const count1: number = poly1.count;
    const vertices1: Vec2[] = poly1.vertices;

    const iv1: number = edge1;
    const iv2: number = edge1 + 1 < count1 ? edge1 + 1 : 0;

    const local_v11: Vec2 = vertices1[iv1];
    const local_v12: Vec2 = vertices1[iv2];

    const localTangent: Vec2 = Vec2.SubVV(local_v12, local_v11, CollidePolygons_s_localTangent);
    localTangent.Normalize();

    const localNormal: Vec2 = Vec2.CrossVOne(localTangent, CollidePolygons_s_localNormal);
    const planePoint: Vec2 = Vec2.MidVV(local_v11, local_v12, CollidePolygons_s_planePoint);

    const tangent: Vec2 = Rot.MulRV(xf1.q, localTangent, CollidePolygons_s_tangent);
    const normal: Vec2 = Vec2.CrossVOne(tangent, CollidePolygons_s_normal);

    const v11: Vec2 = Transform.MulXV(xf1, local_v11, CollidePolygons_s_v11);
    const v12: Vec2 = Transform.MulXV(xf1, local_v12, CollidePolygons_s_v12);

    // Face offset.
    const frontOffset: number = Vec2.DotVV(normal, v11);

    // Side offsets, extended by polytope skin thickness.
    const sideOffset1: number = -Vec2.DotVV(tangent, v11) + totalRadius;
    const sideOffset2: number = Vec2.DotVV(tangent, v12) + totalRadius;

    // Clip incident edge against extruded edge1 side edges.
    const clipPoints1: [ClipVertex, ClipVertex] = CollidePolygons_s_clipPoints1;
    const clipPoints2: [ClipVertex, ClipVertex] = CollidePolygons_s_clipPoints2;
    let np: number;

    // Clip to box side 1
    const ntangent: Vec2 = Vec2.NegV(tangent, CollidePolygons_s_ntangent);
    np = ClipSegmentToLine(clipPoints1, incidentEdge, ntangent, sideOffset1, iv1);

    if (np < 2) {
      return;
    }

    // Clip to negative box side 1
    np = ClipSegmentToLine(clipPoints2, clipPoints1, tangent, sideOffset2, iv2);

    if (np < 2) {
      return;
    }

    // Now clipPoints2 contains the clipped points.
    manifold.localNormal.Copy(localNormal);
    manifold.localPoint.Copy(planePoint);

    let pointCount: number = 0;
    for (let i: number = 0; i < maxManifoldPoints; ++i) {
      const cv: ClipVertex = clipPoints2[i];
      const separation: number = Vec2.DotVV(normal, cv.v) - frontOffset;

      if (separation <= totalRadius) {
        const cp: ManifoldPoint = manifold.points[pointCount];
        Transform.MulTXV(xf2, cv.v, cp.localPoint);
        cp.id.Copy(cv.id);
        if (flip) {
          // Swap features
          const cf: ContactFeature = cp.id.cf;
          cp.id.cf.indexA = cf.indexB;
          cp.id.cf.indexB = cf.indexA;
          cp.id.cf.typeA = cf.typeB;
          cp.id.cf.typeB = cf.typeA;
        }
        ++pointCount;
      }
    }

    manifold.pointCount = pointCount;
  }
}

