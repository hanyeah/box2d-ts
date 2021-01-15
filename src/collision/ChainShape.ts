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

// DEBUG: import { Assert, linearSlop } from "../common/settings";

/// A chain shape is a free form sequence of line segments.
/// The chain has one-sided collision, with the surface normal pointing to the right of the edge.
/// This provides a counter-clockwise winding like the polygon shape.
/// Connectivity information is used to create smooth collisions.
/// @warning the chain will not collide properly if there are self-intersections.

namespace b2 {
  export class ChainShape extends Shape {
    public vertices: Vec2[] = [];
    public count: number = 0;
    public readonly prevVertex: Vec2 = new Vec2();
    public readonly nextVertex: Vec2 = new Vec2();

    constructor() {
      super(ShapeType.e_chainShape, polygonRadius);
    }

    /// Create a loop. This automatically adjusts connectivity.
    /// @param vertices an array of vertices, these are copied
    /// @param count the vertex count
    public CreateLoop(vertices: XY[]): ChainShape;
    public CreateLoop(vertices: XY[], count: number): ChainShape;
    public CreateLoop(vertices: number[]): ChainShape;
    public CreateLoop(...args: any[]): ChainShape {
      if (typeof args[0][0] === "number") {
        const vertices: number[] = args[0];
        if (vertices.length % 2 !== 0) { throw new Error(); }
        return this._CreateLoop((index: number): XY => ({ x: vertices[index * 2], y: vertices[index * 2 + 1] }), vertices.length / 2);
      } else {
        const vertices: XY[] = args[0];
        const count: number = args[1] || vertices.length;
        return this._CreateLoop((index: number): XY => vertices[index], count);
      }
    }
    private _CreateLoop(vertices: (index: number) => XY, count: number): ChainShape {
      // DEBUG: Assert(count >= 3);
      if (count < 3) {
        return this;
      }
      // DEBUG: for (let i: number = 1; i < count; ++i) {
      // DEBUG:   const v1 = vertices[start + i - 1];
      // DEBUG:   const v2 = vertices[start + i];
      // DEBUG:   // If the code crashes here, it means your vertices are too close together.
      // DEBUG:   Assert(Vec2.DistanceSquaredVV(v1, v2) > linearSlop * linearSlop);
      // DEBUG: }

      this.count = count + 1;
      this.vertices = Vec2.MakeArray(this.count);
      for (let i: number = 0; i < count; ++i) {
        this.vertices[i].Copy(vertices(i));
      }
      this.vertices[count].Copy(this.vertices[0]);
      this.prevVertex.Copy(this.vertices[this.count - 2]);
      this.nextVertex.Copy(this.vertices[1]);
      return this;
    }

    /// Create a chain with ghost vertices to connect multiple chains together.
    /// @param vertices an array of vertices, these are copied
    /// @param count the vertex count
    /// @param prevVertex previous vertex from chain that connects to the start
    /// @param nextVertex next vertex from chain that connects to the end
    public CreateChain(vertices: XY[], prevVertex: XY, nextVertex: XY): ChainShape;
    public CreateChain(vertices: XY[], count: number, prevVertex: XY, nextVertex: XY): ChainShape;
    public CreateChain(vertices: number[], prevVertex: XY, nextVertex: XY): ChainShape;
    public CreateChain(...args: any[]): ChainShape {
      if (typeof args[0][0] === "number") {
        const vertices: number[] = args[0];
        const prevVertex: XY = args[1];
        const nextVertex: XY = args[2];
        if (vertices.length % 2 !== 0) { throw new Error(); }
        return this._CreateChain((index: number): XY => ({ x: vertices[index * 2], y: vertices[index * 2 + 1] }), vertices.length / 2, prevVertex, nextVertex);
      } else {
        const vertices: XY[] = args[0];
        const count: number = args[1] || vertices.length;
        const prevVertex: XY = args[2];
        const nextVertex: XY = args[3];
        return this._CreateChain((index: number): XY => vertices[index], count, prevVertex, nextVertex);
      }
    }
    private _CreateChain(vertices: (index: number) => XY, count: number, prevVertex: XY, nextVertex: XY): ChainShape {
      // DEBUG: Assert(count >= 2);
      // DEBUG: for (let i: number = 1; i < count; ++i) {
      // DEBUG:   const v1 = vertices[start + i - 1];
      // DEBUG:   const v2 = vertices[start + i];
      // DEBUG:   // If the code crashes here, it means your vertices are too close together.
      // DEBUG:   Assert(Vec2.DistanceSquaredVV(v1, v2) > linearSlop * linearSlop);
      // DEBUG: }

      this.count = count;
      this.vertices = Vec2.MakeArray(count);
      for (let i: number = 0; i < count; ++i) {
        this.vertices[i].Copy(vertices(i));
      }

      this.prevVertex.Copy(prevVertex);
      this.nextVertex.Copy(nextVertex);

      return this;
    }

    /// Implement Shape. Vertices are cloned using Alloc.
    public Clone(): ChainShape {
      return new ChainShape().Copy(this);
    }

    public Copy(other: ChainShape): ChainShape {
      super.Copy(other);

      // DEBUG: Assert(other instanceof ChainShape);

      this._CreateChain((index: number): XY => other.vertices[index], other.count, other.prevVertex, other.nextVertex);
      this.prevVertex.Copy(other.prevVertex);
      this.nextVertex.Copy(other.nextVertex);

      return this;
    }

    /// @see Shape::GetChildCount
    public GetChildCount(): number {
      // edge count = vertex count - 1
      return this.count - 1;
    }

    /// Get a child edge.
    public GetChildEdge(edge: EdgeShape, index: number): void {
      // DEBUG: Assert(0 <= index && index < this.count - 1);
      edge.radius = this.radius;

      edge.vertex1.Copy(this.vertices[index]);
      edge.vertex2.Copy(this.vertices[index + 1]);
      edge.oneSided = true;

      if (index > 0) {
        edge.vertex0.Copy(this.vertices[index - 1]);
      } else {
        edge.vertex0.Copy(this.prevVertex);
      }

      if (index < this.count - 2) {
        edge.vertex3.Copy(this.vertices[index + 2]);
      } else {
        edge.vertex3.Copy(this.nextVertex);
      }
    }

    /// This always return false.
    /// @see Shape::TestPoint
    public TestPoint(xf: Transform, p: XY): boolean {
      return false;
    }

    // #if ENABLE_PARTICLE
    /// @see Shape::ComputeDistance
    private static ComputeDistance_s_edgeShape = new EdgeShape();
    public ComputeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number {
      const edge = ChainShape.ComputeDistance_s_edgeShape;
      this.GetChildEdge(edge, childIndex);
      return edge.ComputeDistance(xf, p, normal, 0);
    }
    // #endif

    /// Implement Shape.
    private static RayCast_s_edgeShape = new EdgeShape();
    public RayCast(output: RayCastOutput, input: RayCastInput, xf: Transform, childIndex: number): boolean {
      // DEBUG: Assert(childIndex < this.count);

      const edgeShape: EdgeShape = ChainShape.RayCast_s_edgeShape;

      edgeShape.vertex1.Copy(this.vertices[childIndex]);
      edgeShape.vertex2.Copy(this.vertices[(childIndex + 1) % this.count]);

      return edgeShape.RayCast(output, input, xf, 0);
    }

    /// @see Shape::ComputeAABB
    private static ComputeAABB_s_v1 = new Vec2();
    private static ComputeAABB_s_v2 = new Vec2();
    private static ComputeAABB_s_lower = new Vec2();
    private static ComputeAABB_s_upper = new Vec2();
    public ComputeAABB(aabb: AABB, xf: Transform, childIndex: number): void {
      // DEBUG: Assert(childIndex < this.count);

      const vertexi1: Vec2 = this.vertices[childIndex];
      const vertexi2: Vec2 = this.vertices[(childIndex + 1) % this.count];

      const v1: Vec2 = Transform.MulXV(xf, vertexi1, ChainShape.ComputeAABB_s_v1);
      const v2: Vec2 = Transform.MulXV(xf, vertexi2, ChainShape.ComputeAABB_s_v2);

      const lower: Vec2 = Vec2.MinV(v1, v2, ChainShape.ComputeAABB_s_lower);
      const upper: Vec2 = Vec2.MaxV(v1, v2, ChainShape.ComputeAABB_s_upper);

      aabb.lowerBound.x = lower.x - this.radius;
      aabb.lowerBound.y = lower.y - this.radius;
      aabb.upperBound.x = upper.x + this.radius;
      aabb.upperBound.y = upper.y + this.radius;
    }

    /// Chains have zero mass.
    /// @see Shape::ComputeMass
    public ComputeMass(massData: MassData, density: number): void {
      massData.mass = 0;
      massData.center.SetZero();
      massData.I = 0;
    }

    public SetupDistanceProxy(proxy: DistanceProxy, index: number): void {
      // DEBUG: Assert(0 <= index && index < this.count);

      proxy.vertices = proxy.buffer;
      proxy.vertices[0].Copy(this.vertices[index]);
      if (index + 1 < this.count) {
        proxy.vertices[1].Copy(this.vertices[index + 1]);
      } else {
        proxy.vertices[1].Copy(this.vertices[0]);
      }
      proxy.count = 2;
      proxy.radius = this.radius;
    }

    public ComputeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number {
      c.SetZero();
      return 0;
    }

    public Dump(log: (format: string, ...args: any[]) => void): void {
      log("    const shape: ChainShape = new ChainShape();\n");
      log("    const vs: Vec2[] = [];\n");
      for (let i: number = 0; i < this.count; ++i) {
        log("    vs[%d] = new bVec2(%.15f, %.15f);\n", i, this.vertices[i].x, this.vertices[i].y);
      }
      log("    shape.CreateChain(vs, %d);\n", this.count);
      log("    shape.prevVertex.Set(%.15f, %.15f);\n", this.prevVertex.x, this.prevVertex.y);
      log("    shape.nextVertex.Set(%.15f, %.15f);\n", this.nextVertex.x, this.nextVertex.y);
    }
  }
}

