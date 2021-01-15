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
      super(ShapeType.ChainShape, polygonRadius);
    }

    /// Create a loop. This automatically adjusts connectivity.
    /// @param vertices an array of vertices, these are copied
    /// @param count the vertex count
    public createLoop(vertices: XY[]): ChainShape;
    public createLoop(vertices: XY[], count: number): ChainShape;
    public createLoop(vertices: number[]): ChainShape;
    public createLoop(...args: any[]): ChainShape {
      if (typeof args[0][0] === "number") {
        const vertices: number[] = args[0];
        if (vertices.length % 2 !== 0) { throw new Error(); }
        return this._createLoop((index: number): XY => ({ x: vertices[index * 2], y: vertices[index * 2 + 1] }), vertices.length / 2);
      } else {
        const vertices: XY[] = args[0];
        const count: number = args[1] || vertices.length;
        return this._createLoop((index: number): XY => vertices[index], count);
      }
    }
    private _createLoop(vertices: (index: number) => XY, count: number): ChainShape {
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
        this.vertices[i].copy(vertices(i));
      }
      this.vertices[count].copy(this.vertices[0]);
      this.prevVertex.copy(this.vertices[this.count - 2]);
      this.nextVertex.copy(this.vertices[1]);
      return this;
    }

    /// Create a chain with ghost vertices to connect multiple chains together.
    /// @param vertices an array of vertices, these are copied
    /// @param count the vertex count
    /// @param prevVertex previous vertex from chain that connects to the start
    /// @param nextVertex next vertex from chain that connects to the end
    public createChain(vertices: XY[], prevVertex: XY, nextVertex: XY): ChainShape;
    public createChain(vertices: XY[], count: number, prevVertex: XY, nextVertex: XY): ChainShape;
    public createChain(vertices: number[], prevVertex: XY, nextVertex: XY): ChainShape;
    public createChain(...args: any[]): ChainShape {
      if (typeof args[0][0] === "number") {
        const vertices: number[] = args[0];
        const prevVertex: XY = args[1];
        const nextVertex: XY = args[2];
        if (vertices.length % 2 !== 0) { throw new Error(); }
        return this._createChain((index: number): XY => ({ x: vertices[index * 2], y: vertices[index * 2 + 1] }), vertices.length / 2, prevVertex, nextVertex);
      } else {
        const vertices: XY[] = args[0];
        const count: number = args[1] || vertices.length;
        const prevVertex: XY = args[2];
        const nextVertex: XY = args[3];
        return this._createChain((index: number): XY => vertices[index], count, prevVertex, nextVertex);
      }
    }
    private _createChain(vertices: (index: number) => XY, count: number, prevVertex: XY, nextVertex: XY): ChainShape {
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
        this.vertices[i].copy(vertices(i));
      }

      this.prevVertex.copy(prevVertex);
      this.nextVertex.copy(nextVertex);

      return this;
    }

    /// Implement Shape. Vertices are cloned using Alloc.
    public clone(): ChainShape {
      return new ChainShape().copy(this);
    }

    public copy(other: ChainShape): ChainShape {
      super.copy(other);

      // DEBUG: Assert(other instanceof ChainShape);

      this._createChain((index: number): XY => other.vertices[index], other.count, other.prevVertex, other.nextVertex);
      this.prevVertex.copy(other.prevVertex);
      this.nextVertex.copy(other.nextVertex);

      return this;
    }

    /// @see Shape::GetChildCount
    public getChildCount(): number {
      // edge count = vertex count - 1
      return this.count - 1;
    }

    /// Get a child edge.
    public getChildEdge(edge: EdgeShape, index: number): void {
      // DEBUG: Assert(0 <= index && index < this.count - 1);
      edge.radius = this.radius;

      edge.vertex1.copy(this.vertices[index]);
      edge.vertex2.copy(this.vertices[index + 1]);
      edge.oneSided = true;

      if (index > 0) {
        edge.vertex0.copy(this.vertices[index - 1]);
      } else {
        edge.vertex0.copy(this.prevVertex);
      }

      if (index < this.count - 2) {
        edge.vertex3.copy(this.vertices[index + 2]);
      } else {
        edge.vertex3.copy(this.nextVertex);
      }
    }

    /// This always return false.
    /// @see Shape::TestPoint
    public testPoint(xf: Transform, p: XY): boolean {
      return false;
    }

    // #if ENABLE_PARTICLE
    /// @see Shape::ComputeDistance
    private static ComputeDistance_s_edgeShape = new EdgeShape();
    public computeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number {
      const edge = ChainShape.ComputeDistance_s_edgeShape;
      this.getChildEdge(edge, childIndex);
      return edge.computeDistance(xf, p, normal, 0);
    }
    // #endif

    /// Implement Shape.
    private static rayCast_s_edgeShape = new EdgeShape();
    public rayCast(output: RayCastOutput, input: RayCastInput, xf: Transform, childIndex: number): boolean {
      // DEBUG: Assert(childIndex < this.count);

      const edgeShape: EdgeShape = ChainShape.rayCast_s_edgeShape;

      edgeShape.vertex1.copy(this.vertices[childIndex]);
      edgeShape.vertex2.copy(this.vertices[(childIndex + 1) % this.count]);

      return edgeShape.rayCast(output, input, xf, 0);
    }

    /// @see Shape::ComputeAABB
    private static computeAABB_s_v1 = new Vec2();
    private static computeAABB_s_v2 = new Vec2();
    private static computeAABB_s_lower = new Vec2();
    private static computeAABB_s_upper = new Vec2();
    public computeAABB(aabb: AABB, xf: Transform, childIndex: number): void {
      // DEBUG: Assert(childIndex < this.count);

      const vertexi1: Vec2 = this.vertices[childIndex];
      const vertexi2: Vec2 = this.vertices[(childIndex + 1) % this.count];

      const v1: Vec2 = Transform.mulXV(xf, vertexi1, ChainShape.computeAABB_s_v1);
      const v2: Vec2 = Transform.mulXV(xf, vertexi2, ChainShape.computeAABB_s_v2);

      const lower: Vec2 = Vec2.MinV(v1, v2, ChainShape.computeAABB_s_lower);
      const upper: Vec2 = Vec2.MaxV(v1, v2, ChainShape.computeAABB_s_upper);

      aabb.lowerBound.x = lower.x - this.radius;
      aabb.lowerBound.y = lower.y - this.radius;
      aabb.upperBound.x = upper.x + this.radius;
      aabb.upperBound.y = upper.y + this.radius;
    }

    /// Chains have zero mass.
    /// @see Shape::ComputeMass
    public computeMass(massData: MassData, density: number): void {
      massData.mass = 0;
      massData.center.setZero();
      massData.I = 0;
    }

    public setupDistanceProxy(proxy: DistanceProxy, index: number): void {
      // DEBUG: Assert(0 <= index && index < this.count);

      proxy.vertices = proxy.buffer;
      proxy.vertices[0].copy(this.vertices[index]);
      if (index + 1 < this.count) {
        proxy.vertices[1].copy(this.vertices[index + 1]);
      } else {
        proxy.vertices[1].copy(this.vertices[0]);
      }
      proxy.count = 2;
      proxy.radius = this.radius;
    }

    public computeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number {
      c.setZero();
      return 0;
    }

    public dump(log: (format: string, ...args: any[]) => void): void {
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

