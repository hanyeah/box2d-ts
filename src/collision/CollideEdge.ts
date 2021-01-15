// DEBUG: import { Assert } from "../common/settings";

namespace b2 {
  const CollideEdgeAndCircle_s_Q: Vec2 = new Vec2();
  const CollideEdgeAndCircle_s_e: Vec2 = new Vec2();
  const CollideEdgeAndCircle_s_d: Vec2 = new Vec2();
  const CollideEdgeAndCircle_s_e1: Vec2 = new Vec2();
  const CollideEdgeAndCircle_s_e2: Vec2 = new Vec2();
  const CollideEdgeAndCircle_s_P: Vec2 = new Vec2();
  const CollideEdgeAndCircle_s_n: Vec2 = new Vec2();
  const CollideEdgeAndCircle_s_id: ContactID = new ContactID();
  export function CollideEdgeAndCircle(manifold: Manifold, edgeA: EdgeShape, xfA: Transform, circleB: CircleShape, xfB: Transform): void {
    manifold.pointCount = 0;

    // Compute circle in frame of edge
    const Q: Vec2 = Transform.MulTXV(xfA, Transform.MulXV(xfB, circleB.p, Vec2.s_t0), CollideEdgeAndCircle_s_Q);

    const A: Vec2 = edgeA.vertex1;
    const B: Vec2 = edgeA.vertex2;
    const e: Vec2 = Vec2.SubVV(B, A, CollideEdgeAndCircle_s_e);

    // Normal points to the right for a CCW winding
    // Vec2 n(e.y, -e.x);
    // const n: Vec2 = CollideEdgeAndCircle_s_n.Set(-e.y, e.x);
    const n: Vec2 = CollideEdgeAndCircle_s_n.Set(e.y, -e.x);
    // float offset = Dot(n, Q - A);
    const offset: number = Vec2.DotVV(n, Vec2.SubVV(Q, A, Vec2.s_t0));

    const oneSided: boolean = edgeA.oneSided;
    if (oneSided && offset < 0.0) {
      return;
    }

    // Barycentric coordinates
    const u: number = Vec2.DotVV(e, Vec2.SubVV(B, Q, Vec2.s_t0));
    const v: number = Vec2.DotVV(e, Vec2.SubVV(Q, A, Vec2.s_t0));

    const radius: number = edgeA.radius + circleB.radius;

    // const cf: ContactFeature = new ContactFeature();
    const id: ContactID = CollideEdgeAndCircle_s_id;
    id.cf.indexB = 0;
    id.cf.typeB = ContactFeatureType.e_vertex;

    // Region A
    if (v <= 0) {
      const P: Vec2 = A;
      const d: Vec2 = Vec2.SubVV(Q, P, CollideEdgeAndCircle_s_d);
      const dd: number = Vec2.DotVV(d, d);
      if (dd > radius * radius) {
        return;
      }

      // Is there an edge connected to A?
      if (edgeA.oneSided) {
        const A1: Vec2 = edgeA.vertex0;
        const B1: Vec2 = A;
        const e1: Vec2 = Vec2.SubVV(B1, A1, CollideEdgeAndCircle_s_e1);
        const u1: number = Vec2.DotVV(e1, Vec2.SubVV(B1, Q, Vec2.s_t0));

        // Is the circle in Region AB of the previous edge?
        if (u1 > 0) {
          return;
        }
      }

      id.cf.indexA = 0;
      id.cf.typeA = ContactFeatureType.e_vertex;
      manifold.pointCount = 1;
      manifold.type = ManifoldType.e_circles;
      manifold.localNormal.SetZero();
      manifold.localPoint.Copy(P);
      manifold.points[0].id.Copy(id);
      // manifold.points[0].id.key = 0;
      // manifold.points[0].id.cf = cf;
      manifold.points[0].localPoint.Copy(circleB.p);
      return;
    }

    // Region B
    if (u <= 0) {
      const P: Vec2 = B;
      const d: Vec2 = Vec2.SubVV(Q, P, CollideEdgeAndCircle_s_d);
      const dd: number = Vec2.DotVV(d, d);
      if (dd > radius * radius) {
        return;
      }

      // Is there an edge connected to B?
      if (edgeA.oneSided) {
        const v: Vec2 = edgeA.vertex3;
        const A2: Vec2 = B;
        const e2: Vec2 = Vec2.SubVV(v, A2, CollideEdgeAndCircle_s_e2);
        const v2: number = Vec2.DotVV(e2, Vec2.SubVV(Q, A2, Vec2.s_t0));

        // Is the circle in Region AB of the next edge?
        if (v2 > 0) {
          return;
        }
      }

      id.cf.indexA = 1;
      id.cf.typeA = ContactFeatureType.e_vertex;
      manifold.pointCount = 1;
      manifold.type = ManifoldType.e_circles;
      manifold.localNormal.SetZero();
      manifold.localPoint.Copy(P);
      manifold.points[0].id.Copy(id);
      // manifold.points[0].id.key = 0;
      // manifold.points[0].id.cf = cf;
      manifold.points[0].localPoint.Copy(circleB.p);
      return;
    }

    // Region AB
    const den: number = Vec2.DotVV(e, e);
    // DEBUG: Assert(den > 0);
    const P: Vec2 = CollideEdgeAndCircle_s_P;
    P.x = (1 / den) * (u * A.x + v * B.x);
    P.y = (1 / den) * (u * A.y + v * B.y);
    const d: Vec2 = Vec2.SubVV(Q, P, CollideEdgeAndCircle_s_d);
    const dd: number = Vec2.DotVV(d, d);
    if (dd > radius * radius) {
      return;
    }

    if (offset < 0) {
      n.Set(-n.x, -n.y);
    }
    n.Normalize();

    id.cf.indexA = 0;
    id.cf.typeA = ContactFeatureType.e_face;
    manifold.pointCount = 1;
    manifold.type = ManifoldType.e_faceA;
    manifold.localNormal.Copy(n);
    manifold.localPoint.Copy(A);
    manifold.points[0].id.Copy(id);
    // manifold.points[0].id.key = 0;
    // manifold.points[0].id.cf = cf;
    manifold.points[0].localPoint.Copy(circleB.p);
  }

  enum EPAxisType {
    e_unknown = 0,
    e_edgeA = 1,
    e_edgeB = 2,
  }

  class EPAxis {
    public normal: Vec2 = new Vec2();
    public type: EPAxisType = EPAxisType.e_unknown;
    public index: number = 0;
    public separation: number = 0;
  }

  class TempPolygon {
    public vertices: Vec2[] = [];
    public normals: Vec2[] = [];
    public count: number = 0;
  }

  class ReferenceFace {
    public i1: number = 0;
    public i2: number = 0;
    public readonly v1: Vec2 = new Vec2();
    public readonly v2: Vec2 = new Vec2();
    public readonly normal: Vec2 = new Vec2();
    public readonly sideNormal1: Vec2 = new Vec2();
    public sideOffset1: number = 0;
    public readonly sideNormal2: Vec2 = new Vec2();
    public sideOffset2: number = 0;
  }

// static EPAxis ComputeEdgeSeparation(const TempPolygon& polygonB, const Vec2& v1, const Vec2& normal1)
  const ComputeEdgeSeparation_s_axis = new EPAxis();
  const ComputeEdgeSeparation_s_axes: [ Vec2, Vec2 ] = [ new Vec2(), new Vec2() ];
  function ComputeEdgeSeparation(polygonB: TempPolygon, v1: Vec2, normal1: Vec2): EPAxis {
    // EPAxis axis;
    const axis: EPAxis = ComputeEdgeSeparation_s_axis;
    axis.type = EPAxisType.e_edgeA;
    axis.index = -1;
    axis.separation = -Number.MAX_VALUE; // -FLT_MAX;
    axis.normal.SetZero();

    // Vec2 axes[2] = { normal1, -normal1 };
    const axes: [Vec2, Vec2] = ComputeEdgeSeparation_s_axes;
    axes[0].Copy(normal1);
    axes[1].Copy(normal1).SelfNeg();

    // Find axis with least overlap (min-max problem)
    for (let j = 0; j < 2; ++j) {
      let sj: number = Number.MAX_VALUE; // FLT_MAX;

      // Find deepest polygon vertex along axis j
      for (let i = 0; i < polygonB.count; ++i) {
        // float si = Dot(axes[j], polygonB.vertices[i] - v1);
        const si: number = Vec2.DotVV(axes[j], Vec2.SubVV(polygonB.vertices[i], v1, Vec2.s_t0));
        if (si < sj) {
          sj = si;
        }
      }

      if (sj > axis.separation) {
        axis.index = j;
        axis.separation = sj;
        axis.normal.Copy(axes[j]);
      }
    }

    return axis;
  }

// static EPAxis ComputePolygonSeparation(const TempPolygon& polygonB, const Vec2& v1, const Vec2& v2)
  const ComputePolygonSeparation_s_axis = new EPAxis();
  const ComputePolygonSeparation_s_n = new Vec2();
  function ComputePolygonSeparation(polygonB: TempPolygon, v1: Vec2, v2: Vec2): EPAxis {
    const axis: EPAxis = ComputePolygonSeparation_s_axis;
    axis.type = EPAxisType.e_unknown;
    axis.index = -1;
    axis.separation = -Number.MAX_VALUE; // -FLT_MAX;
    axis.normal.SetZero();

    for (let i = 0; i < polygonB.count; ++i) {
      // Vec2 n = -polygonB.normals[i];
      const n: Vec2 = Vec2.NegV(polygonB.normals[i], ComputePolygonSeparation_s_n);

      // float s1 = Dot(n, polygonB.vertices[i] - v1);
      const s1: number = Vec2.DotVV(n, Vec2.SubVV(polygonB.vertices[i], v1, Vec2.s_t0));
      // float s2 = Dot(n, polygonB.vertices[i] - v2);
      const s2: number = Vec2.DotVV(n, Vec2.SubVV(polygonB.vertices[i], v2, Vec2.s_t0));
      // float s = Min(s1, s2);
      const s: number = Min(s1, s2);

      if (s > axis.separation) {
        axis.type = EPAxisType.e_edgeB;
        axis.index = i;
        axis.separation = s;
        axis.normal.Copy(n);
      }
    }

    return axis;
  }

  const CollideEdgeAndPolygon_s_xf = new Transform();
  const CollideEdgeAndPolygon_s_centroidB = new Vec2();
  const CollideEdgeAndPolygon_s_edge1 = new Vec2();
  const CollideEdgeAndPolygon_s_normal1 = new Vec2();
  const CollideEdgeAndPolygon_s_edge0 = new Vec2();
  const CollideEdgeAndPolygon_s_normal0 = new Vec2();
  const CollideEdgeAndPolygon_s_edge2 = new Vec2();
  const CollideEdgeAndPolygon_s_normal2 = new Vec2();
  const CollideEdgeAndPolygon_s_tempPolygonB = new TempPolygon();
  const CollideEdgeAndPolygon_s_ref = new ReferenceFace();
  const CollideEdgeAndPolygon_s_clipPoints: [ClipVertex, ClipVertex] = [ new ClipVertex(), new ClipVertex() ];
  const CollideEdgeAndPolygon_s_clipPoints1: [ClipVertex, ClipVertex] = [ new ClipVertex(), new ClipVertex() ];
  const CollideEdgeAndPolygon_s_clipPoints2: [ClipVertex, ClipVertex] = [ new ClipVertex(), new ClipVertex() ];
  export function CollideEdgeAndPolygon(manifold: Manifold, edgeA: EdgeShape, xfA: Transform, polygonB: PolygonShape, xfB: Transform): void {
    manifold.pointCount = 0;

    // Transform xf = MulT(xfA, xfB);
    const xf = Transform.MulTXX(xfA, xfB, CollideEdgeAndPolygon_s_xf);

    // Vec2 centroidB = Mul(xf, polygonB.centroid);
    const centroidB: Vec2 = Transform.MulXV(xf, polygonB.centroid, CollideEdgeAndPolygon_s_centroidB);

    // Vec2 v1 = edgeA.vertex1;
    const v1: Vec2 = edgeA.vertex1;
    // Vec2 v2 = edgeA.vertex2;
    const v2: Vec2 = edgeA.vertex2;

    // Vec2 edge1 = v2 - v1;
    const edge1: Vec2 = Vec2.SubVV(v2, v1, CollideEdgeAndPolygon_s_edge1);
    edge1.Normalize();

    // Normal points to the right for a CCW winding
    // Vec2 normal1(edge1.y, -edge1.x);
    const normal1 = CollideEdgeAndPolygon_s_normal1.Set(edge1.y, -edge1.x);
    // float offset1 = Dot(normal1, centroidB - v1);
    const offset1: number = Vec2.DotVV(normal1, Vec2.SubVV(centroidB, v1, Vec2.s_t0));

    const oneSided: boolean = edgeA.oneSided;
    if (oneSided && offset1 < 0.0) {
      return;
    }

    // Get polygonB in frameA
    // TempPolygon tempPolygonB;
    const tempPolygonB: TempPolygon = CollideEdgeAndPolygon_s_tempPolygonB;
    tempPolygonB.count = polygonB.count;
    for (let i = 0; i < polygonB.count; ++i) {
      if (tempPolygonB.vertices.length <= i) { tempPolygonB.vertices.push(new Vec2()); }
      if (tempPolygonB.normals.length <= i) { tempPolygonB.normals.push(new Vec2()); }
      // tempPolygonB.vertices[i] = Mul(xf, polygonB.vertices[i]);
      Transform.MulXV(xf, polygonB.vertices[i], tempPolygonB.vertices[i]);
      // tempPolygonB.normals[i] = Mul(xf.q, polygonB.normals[i]);
      Rot.MulRV(xf.q, polygonB.normals[i], tempPolygonB.normals[i]);
    }

    const radius: number = polygonB.radius + edgeA.radius;

    // EPAxis edgeAxis = ComputeEdgeSeparation(tempPolygonB, v1, normal1);
    const edgeAxis: EPAxis = ComputeEdgeSeparation(tempPolygonB, v1, normal1);
    if (edgeAxis.separation > radius) {
      return;
    }

    // EPAxis polygonAxis = ComputePolygonSeparation(tedge0.y, -edge0.xempPolygonB, v1, v2);
    const polygonAxis: EPAxis = ComputePolygonSeparation(tempPolygonB, v1, v2);
    if (polygonAxis.separation > radius) {
      return;
    }

    // Use hysteresis for jitter reduction.
    const k_relativeTol: number = 0.98;
    const k_absoluteTol: number = 0.001;

    // EPAxis primaryAxis;
    let primaryAxis: EPAxis;
    if (polygonAxis.separation - radius > k_relativeTol * (edgeAxis.separation - radius) + k_absoluteTol) {
      primaryAxis = polygonAxis;
    } else {
      primaryAxis = edgeAxis;
    }

    if (oneSided) {
      // Smooth collision
      // See https://box2d.org/posts/2020/06/ghost-collisions/

      // Vec2 edge0 = v1 - edgeA.vertex0;
      const edge0: Vec2 = Vec2.SubVV(v1, edgeA.vertex0, CollideEdgeAndPolygon_s_edge0);
      edge0.Normalize();
      // Vec2 normal0(edge0.y, -edge0.x);
      const normal0: Vec2 = CollideEdgeAndPolygon_s_normal0.Set(edge0.y, -edge0.x);
      const convex1: boolean = Vec2.CrossVV(edge0, edge1) >= 0.0;

      // Vec2 edge2 = edgeA.vertex3 - v2;
      const edge2: Vec2 = Vec2.SubVV(edgeA.vertex3, v2, CollideEdgeAndPolygon_s_edge2);
      edge2.Normalize();
      // Vec2 normal2(edge2.y, -edge2.x);
      const normal2: Vec2 = CollideEdgeAndPolygon_s_normal2.Set(edge2.y, -edge2.x);
      const convex2: boolean = Vec2.CrossVV(edge1, edge2) >= 0.0;

      const sinTol: number = 0.1;
      const side1: boolean = Vec2.DotVV(primaryAxis.normal, edge1) <= 0.0;

      // Check Gauss Map
      if (side1) {
        if (convex1) {
          if (Vec2.CrossVV(primaryAxis.normal, normal0) > sinTol) {
            // Skip region
            return;
          }

          // Admit region
        } else {
          // Snap region
          primaryAxis = edgeAxis;
        }
      } else {
        if (convex2) {
          if (Vec2.CrossVV(normal2, primaryAxis.normal) > sinTol) {
            // Skip region
            return;
          }

          // Admit region
        } else {
          // Snap region
          primaryAxis = edgeAxis;
        }
      }
    }

    // ClipVertex clipPoints[2];
    const clipPoints: [ClipVertex, ClipVertex] = CollideEdgeAndPolygon_s_clipPoints;
    // ReferenceFace ref;
    const ref: ReferenceFace = CollideEdgeAndPolygon_s_ref;
    if (primaryAxis.type === EPAxisType.e_edgeA) {
      manifold.type = ManifoldType.e_faceA;

      // Search for the polygon normal that is most anti-parallel to the edge normal.
      let bestIndex: number = 0;
      let bestValue: number = Vec2.DotVV(primaryAxis.normal, tempPolygonB.normals[0]);
      for (let i = 1; i < tempPolygonB.count; ++i) {
        const value: number = Vec2.DotVV(primaryAxis.normal, tempPolygonB.normals[i]);
        if (value < bestValue) {
          bestValue = value;
          bestIndex = i;
        }
      }

      const i1: number = bestIndex;
      const i2: number = i1 + 1 < tempPolygonB.count ? i1 + 1 : 0;

      clipPoints[0].v.Copy(tempPolygonB.vertices[i1]);
      clipPoints[0].id.cf.indexA = 0;
      clipPoints[0].id.cf.indexB = i1;
      clipPoints[0].id.cf.typeA = ContactFeatureType.e_face;
      clipPoints[0].id.cf.typeB = ContactFeatureType.e_vertex;

      clipPoints[1].v.Copy(tempPolygonB.vertices[i2]);
      clipPoints[1].id.cf.indexA = 0;
      clipPoints[1].id.cf.indexB = i2;
      clipPoints[1].id.cf.typeA = ContactFeatureType.e_face;
      clipPoints[1].id.cf.typeB = ContactFeatureType.e_vertex;

      ref.i1 = 0;
      ref.i2 = 1;
      ref.v1.Copy(v1);
      ref.v2.Copy(v2);
      ref.normal.Copy(primaryAxis.normal);
      ref.sideNormal1.Copy(edge1).SelfNeg(); // ref.sideNormal1 = -edge1;
      ref.sideNormal2.Copy(edge1);
    } else {
      manifold.type = ManifoldType.e_faceB;

      clipPoints[0].v.Copy(v2);
      clipPoints[0].id.cf.indexA = 1;
      clipPoints[0].id.cf.indexB = primaryAxis.index;
      clipPoints[0].id.cf.typeA = ContactFeatureType.e_vertex;
      clipPoints[0].id.cf.typeB = ContactFeatureType.e_face;

      clipPoints[1].v.Copy(v1);
      clipPoints[1].id.cf.indexA = 0;
      clipPoints[1].id.cf.indexB = primaryAxis.index;
      clipPoints[1].id.cf.typeA = ContactFeatureType.e_vertex;
      clipPoints[1].id.cf.typeB = ContactFeatureType.e_face;

      ref.i1 = primaryAxis.index;
      ref.i2 = ref.i1 + 1 < tempPolygonB.count ? ref.i1 + 1 : 0;
      ref.v1.Copy(tempPolygonB.vertices[ref.i1]);
      ref.v2.Copy(tempPolygonB.vertices[ref.i2]);
      ref.normal.Copy(tempPolygonB.normals[ref.i1]);

      // CCW winding
      ref.sideNormal1.Set(ref.normal.y, -ref.normal.x);
      ref.sideNormal2.Copy(ref.sideNormal1).SelfNeg(); // ref.sideNormal2 = -ref.sideNormal1;
    }

    ref.sideOffset1 = Vec2.DotVV(ref.sideNormal1, ref.v1);
    ref.sideOffset2 = Vec2.DotVV(ref.sideNormal2, ref.v2);

    // Clip incident edge against reference face side planes
    // ClipVertex clipPoints1[2];
    const clipPoints1: [ClipVertex, ClipVertex] = CollideEdgeAndPolygon_s_clipPoints1; // [new ClipVertex(), new ClipVertex()];
    // ClipVertex clipPoints2[2];
    const clipPoints2: [ClipVertex, ClipVertex] = CollideEdgeAndPolygon_s_clipPoints2; // [new ClipVertex(), new ClipVertex()];
    // int32 np;
    let np: number;

    // Clip to side 1
    np = ClipSegmentToLine(clipPoints1, clipPoints, ref.sideNormal1, ref.sideOffset1, ref.i1);

    if (np < maxManifoldPoints) {
      return;
    }

    // Clip to side 2
    np = ClipSegmentToLine(clipPoints2, clipPoints1, ref.sideNormal2, ref.sideOffset2, ref.i2);

    if (np < maxManifoldPoints) {
      return;
    }

    // Now clipPoints2 contains the clipped points.
    if (primaryAxis.type === EPAxisType.e_edgeA) {
      manifold.localNormal.Copy(ref.normal);
      manifold.localPoint.Copy(ref.v1);
    } else {
      manifold.localNormal.Copy(polygonB.normals[ref.i1]);
      manifold.localPoint.Copy(polygonB.vertices[ref.i1]);
    }

    let pointCount = 0;
    for (let i = 0; i < maxManifoldPoints; ++i) {
      const separation: number = Vec2.DotVV(ref.normal, Vec2.SubVV(clipPoints2[i].v, ref.v1, Vec2.s_t0));

      if (separation <= radius) {
        const cp: ManifoldPoint = manifold.points[pointCount];

        if (primaryAxis.type === EPAxisType.e_edgeA) {
          Transform.MulTXV(xf, clipPoints2[i].v, cp.localPoint); // cp.localPoint = MulT(xf, clipPoints2[i].v);
          cp.id.Copy(clipPoints2[i].id);
        } else {
          cp.localPoint.Copy(clipPoints2[i].v);
          cp.id.cf.typeA = clipPoints2[i].id.cf.typeB;
          cp.id.cf.typeB = clipPoints2[i].id.cf.typeA;
          cp.id.cf.indexA = clipPoints2[i].id.cf.indexB;
          cp.id.cf.indexB = clipPoints2[i].id.cf.indexA;
        }

        ++pointCount;
      }
    }

    manifold.pointCount = pointCount;
  }

}
