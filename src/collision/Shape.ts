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

/// This holds the mass data computed for a shape.
  export class MassData {
    /// The mass of the shape, usually in kilograms.
    public mass: number = 0;

    /// The position of the shape's centroid relative to the shape's origin.
    public readonly center: Vec2 = new Vec2(0, 0);

    /// The rotational inertia of the shape about the local origin.
    public I: number = 0;
  }

  export enum ShapeType {
    Unknown = -1,
    CircleShape = 0,
    EdgeShape = 1,
    PolygonShape = 2,
    ChainShape = 3,
    ShapeTypeCount = 4,
  }

/// A shape is used for collision detection. You can create a shape however you like.
/// Shapes used for simulation in World are created automatically when a Fixture
/// is created. Shapes may encapsulate a one or more child shapes.
  export abstract class Shape {
    public readonly type: ShapeType = ShapeType.Unknown;

    /// Radius of a shape. For polygonal shapes this must be polygonRadius. There is no support for
    /// making rounded polygons.
    public radius: number = 0;

    constructor(type: ShapeType, radius: number) {
      this.type = type;
      this.radius = radius;
    }

    /// Clone the concrete shape.
    public abstract clone(): Shape;

    public copy(other: Shape): Shape {
      // DEBUG: Assert(this.type === other.type);
      this.radius = other.radius;
      return this;
    }

    /// Get the type of this shape. You can use this to down cast to the concrete shape.
    /// @return the shape type.
    public getType(): ShapeType {
      return this.type;
    }

    /// Get the number of child primitives.
    public abstract getChildCount(): number;

    /// Test a point for containment in this shape. This only works for convex shapes.
    /// @param xf the shape world transform.
    /// @param p a point in world coordinates.
    public abstract testPoint(xf: Transform, p: XY): boolean;

    // #if ENABLE_PARTICLE
    /// Compute the distance from the current shape to the specified point. This only works for convex shapes.
    /// @param xf the shape world transform.
    /// @param p a point in world coordinates.
    /// @param distance returns the distance from the current shape.
    /// @param normal returns the direction in which the distance increases.
    public abstract computeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number;
    // #endif

    /// Cast a ray against a child shape.
    /// @param output the ray-cast results.
    /// @param input the ray-cast input parameters.
    /// @param transform the transform to be applied to the shape.
    /// @param childIndex the child shape index
    public abstract rayCast(output: RayCastOutput, input: RayCastInput, transform: Transform, childIndex: number): boolean;

    /// Given a transform, compute the associated axis aligned bounding box for a child shape.
    /// @param aabb returns the axis aligned box.
    /// @param xf the world transform of the shape.
    /// @param childIndex the child shape
    public abstract computeAABB(aabb: AABB, xf: Transform, childIndex: number): void;

    /// Compute the mass properties of this shape using its dimensions and density.
    /// The inertia tensor is computed about the local origin.
    /// @param massData returns the mass data for this shape.
    /// @param density the density in kilograms per meter squared.
    public abstract computeMass(massData: MassData, density: number): void;

    public abstract setupDistanceProxy(proxy: DistanceProxy, index: number): void;

    public abstract computeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number;

    public abstract dump(log: (format: string, ...args: any[]) => void): void;
  }
}

