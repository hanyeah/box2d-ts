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

  export let toiTime: number = 0;
  export let toiMaxTime: number = 0;
  export let toiCalls: number = 0;
  export let toiIters: number = 0;
  export let toiMaxIters: number = 0;
  export let toiRootIters: number = 0;
  export let toiMaxRootIters: number = 0;
  export function toiReset(): void {
    toiTime = 0;
    toiMaxTime = 0;
    toiCalls = 0;
    toiIters = 0;
    toiMaxIters = 0;
    toiRootIters = 0;
    toiMaxRootIters = 0;
  }

  const timeOfImpact_s_xfA: Transform = new Transform();
  const timeOfImpact_s_xfB: Transform = new Transform();
  const timeOfImpact_s_pointA: Vec2 = new Vec2();
  const timeOfImpact_s_pointB: Vec2 = new Vec2();
  const timeOfImpact_s_normal: Vec2 = new Vec2();
  const timeOfImpact_s_axisA: Vec2 = new Vec2();
  const timeOfImpact_s_axisB: Vec2 = new Vec2();

/// Input parameters for TimeOfImpact
  export class TOIInput {
    public readonly proxyA: DistanceProxy = new DistanceProxy();
    public readonly proxyB: DistanceProxy = new DistanceProxy();
    public readonly sweepA: Sweep = new Sweep();
    public readonly sweepB: Sweep = new Sweep();
    public tMax: number = 0; // defines sweep interval [0, tMax]
  }

/// Output parameters for TimeOfImpact.
  export enum TOIOutputState {
    Unknown = 0,
    Failed = 1,
    Overlapped = 2,
    Touching = 3,
    Separated = 4,
  }

  export class TOIOutput {
    public state = TOIOutputState.Unknown;
    public t: number = 0;
  }

  export enum SeparationFunctionType {
    Unknown = -1,
    Points = 0,
    FaceA = 1,
    FaceB = 2,
  }

  export class SeparationFunction {
    public proxyA!: DistanceProxy;
    public proxyB!: DistanceProxy;
    public readonly sweepA: Sweep = new Sweep();
    public readonly sweepB: Sweep = new Sweep();
    public type: SeparationFunctionType = SeparationFunctionType.Unknown;
    public readonly localPoint: Vec2 = new Vec2();
    public readonly axis: Vec2 = new Vec2();

    public initialize(cache: SimplexCache, proxyA: DistanceProxy, sweepA: Sweep, proxyB: DistanceProxy, sweepB: Sweep, t1: number): number {
      this.proxyA = proxyA;
      this.proxyB = proxyB;
      const count: number = cache.count;
      // DEBUG: Assert(0 < count && count < 3);

      this.sweepA.copy(sweepA);
      this.sweepB.copy(sweepB);

      const xfA: Transform = timeOfImpact_s_xfA;
      const xfB: Transform = timeOfImpact_s_xfB;
      this.sweepA.getTransform(xfA, t1);
      this.sweepB.getTransform(xfB, t1);

      if (count === 1) {
        this.type = SeparationFunctionType.Points;
        const localPointA: Vec2 = this.proxyA.getVertex(cache.indexA[0]);
        const localPointB: Vec2 = this.proxyB.getVertex(cache.indexB[0]);
        const pointA: Vec2 = Transform.mulXV(xfA, localPointA, timeOfImpact_s_pointA);
        const pointB: Vec2 = Transform.mulXV(xfB, localPointB, timeOfImpact_s_pointB);
        Vec2.SubVV(pointB, pointA, this.axis);
        const s: number = this.axis.normalize();
        // #if ENABLE_PARTICLE
        this.localPoint.setZero();
        // #endif
        return s;
      } else if (cache.indexA[0] === cache.indexA[1]) {
        // Two points on B and one on A.
        this.type = SeparationFunctionType.FaceB;
        const localPointB1: Vec2 = this.proxyB.getVertex(cache.indexB[0]);
        const localPoint: Vec2 = this.proxyB.getVertex(cache.indexB[1]);

        Vec2.CrossVOne(Vec2.SubVV(localPoint, localPointB1, Vec2.s_t0), this.axis).selfNormalize();
        const normal: Vec2 = Rot.mulRV(xfB.q, this.axis, timeOfImpact_s_normal);

        Vec2.MidVV(localPointB1, localPoint, this.localPoint);
        const pointB: Vec2 = Transform.mulXV(xfB, this.localPoint, timeOfImpact_s_pointB);

        const localPointA: Vec2 = this.proxyA.getVertex(cache.indexA[0]);
        const pointA: Vec2 = Transform.mulXV(xfA, localPointA, timeOfImpact_s_pointA);

        let s: number = Vec2.DotVV(Vec2.SubVV(pointA, pointB, Vec2.s_t0), normal);
        if (s < 0) {
          this.axis.selfNeg();
          s = -s;
        }
        return s;
      } else {
        // Two points on A and one or two points on B.
        this.type = SeparationFunctionType.FaceA;
        const localPointA1: Vec2 = this.proxyA.getVertex(cache.indexA[0]);
        const localPointA2: Vec2 = this.proxyA.getVertex(cache.indexA[1]);

        Vec2.CrossVOne(Vec2.SubVV(localPointA2, localPointA1, Vec2.s_t0), this.axis).selfNormalize();
        const normal: Vec2 = Rot.mulRV(xfA.q, this.axis, timeOfImpact_s_normal);

        Vec2.MidVV(localPointA1, localPointA2, this.localPoint);
        const pointA: Vec2 = Transform.mulXV(xfA, this.localPoint, timeOfImpact_s_pointA);

        const localPointB: Vec2 = this.proxyB.getVertex(cache.indexB[0]);
        const pointB: Vec2 = Transform.mulXV(xfB, localPointB, timeOfImpact_s_pointB);

        let s: number = Vec2.DotVV(Vec2.SubVV(pointB, pointA, Vec2.s_t0), normal);
        if (s < 0) {
          this.axis.selfNeg();
          s = -s;
        }
        return s;
      }
    }

    public findMinSeparation(indexA: [number], indexB: [number], t: number): number {
      const xfA: Transform = timeOfImpact_s_xfA;
      const xfB: Transform = timeOfImpact_s_xfB;
      this.sweepA.getTransform(xfA, t);
      this.sweepB.getTransform(xfB, t);

      switch (this.type) {
        case SeparationFunctionType.Points: {
          const axisA: Vec2 = Rot.mulTRV(xfA.q, this.axis, timeOfImpact_s_axisA);
          const axisB: Vec2 = Rot.mulTRV(xfB.q, Vec2.NegV(this.axis, Vec2.s_t0), timeOfImpact_s_axisB);

          indexA[0] = this.proxyA.getSupport(axisA);
          indexB[0] = this.proxyB.getSupport(axisB);

          const localPointA: Vec2 = this.proxyA.getVertex(indexA[0]);
          const localPointB: Vec2 = this.proxyB.getVertex(indexB[0]);

          const pointA: Vec2 = Transform.mulXV(xfA, localPointA, timeOfImpact_s_pointA);
          const pointB: Vec2 = Transform.mulXV(xfB, localPointB, timeOfImpact_s_pointB);

          const separation: number = Vec2.DotVV(Vec2.SubVV(pointB, pointA, Vec2.s_t0), this.axis);
          return separation;
        }

        case SeparationFunctionType.FaceA: {
          const normal: Vec2 = Rot.mulRV(xfA.q, this.axis, timeOfImpact_s_normal);
          const pointA: Vec2 = Transform.mulXV(xfA, this.localPoint, timeOfImpact_s_pointA);

          const axisB: Vec2 = Rot.mulTRV(xfB.q, Vec2.NegV(normal, Vec2.s_t0), timeOfImpact_s_axisB);

          indexA[0] = -1;
          indexB[0] = this.proxyB.getSupport(axisB);

          const localPointB: Vec2 = this.proxyB.getVertex(indexB[0]);
          const pointB: Vec2 = Transform.mulXV(xfB, localPointB, timeOfImpact_s_pointB);

          const separation: number = Vec2.DotVV(Vec2.SubVV(pointB, pointA, Vec2.s_t0), normal);
          return separation;
        }

        case SeparationFunctionType.FaceB: {
          const normal: Vec2 = Rot.mulRV(xfB.q, this.axis, timeOfImpact_s_normal);
          const pointB: Vec2 = Transform.mulXV(xfB, this.localPoint, timeOfImpact_s_pointB);

          const axisA: Vec2 = Rot.mulTRV(xfA.q, Vec2.NegV(normal, Vec2.s_t0), timeOfImpact_s_axisA);

          indexB[0] = -1;
          indexA[0] = this.proxyA.getSupport(axisA);

          const localPointA: Vec2 = this.proxyA.getVertex(indexA[0]);
          const pointA: Vec2 = Transform.mulXV(xfA, localPointA, timeOfImpact_s_pointA);

          const separation: number = Vec2.DotVV(Vec2.SubVV(pointA, pointB, Vec2.s_t0), normal);
          return separation;
        }

        default:
          // DEBUG: Assert(false);
          indexA[0] = -1;
          indexB[0] = -1;
          return 0;
      }
    }

    public evaluate(indexA: number, indexB: number, t: number): number {
      const xfA: Transform = timeOfImpact_s_xfA;
      const xfB: Transform = timeOfImpact_s_xfB;
      this.sweepA.getTransform(xfA, t);
      this.sweepB.getTransform(xfB, t);

      switch (this.type) {
        case SeparationFunctionType.Points: {
          const localPointA: Vec2 = this.proxyA.getVertex(indexA);
          const localPointB: Vec2 = this.proxyB.getVertex(indexB);

          const pointA: Vec2 = Transform.mulXV(xfA, localPointA, timeOfImpact_s_pointA);
          const pointB: Vec2 = Transform.mulXV(xfB, localPointB, timeOfImpact_s_pointB);
          const separation: number = Vec2.DotVV(Vec2.SubVV(pointB, pointA, Vec2.s_t0), this.axis);

          return separation;
        }

        case SeparationFunctionType.FaceA: {
          const normal: Vec2 = Rot.mulRV(xfA.q, this.axis, timeOfImpact_s_normal);
          const pointA: Vec2 = Transform.mulXV(xfA, this.localPoint, timeOfImpact_s_pointA);

          const localPointB: Vec2 = this.proxyB.getVertex(indexB);
          const pointB: Vec2 = Transform.mulXV(xfB, localPointB, timeOfImpact_s_pointB);

          const separation: number = Vec2.DotVV(Vec2.SubVV(pointB, pointA, Vec2.s_t0), normal);
          return separation;
        }

        case SeparationFunctionType.FaceB: {
          const normal: Vec2 = Rot.mulRV(xfB.q, this.axis, timeOfImpact_s_normal);
          const pointB: Vec2 = Transform.mulXV(xfB, this.localPoint, timeOfImpact_s_pointB);

          const localPointA: Vec2 = this.proxyA.getVertex(indexA);
          const pointA: Vec2 = Transform.mulXV(xfA, localPointA, timeOfImpact_s_pointA);

          const separation: number = Vec2.DotVV(Vec2.SubVV(pointA, pointB, Vec2.s_t0), normal);
          return separation;
        }

        default:
          // DEBUG: Assert(false);
          return 0;
      }
    }
  }

  const timeOfImpact_s_timer: Timer = new Timer();
  const timeOfImpact_s_cache: SimplexCache = new SimplexCache();
  const timeOfImpact_s_distanceInput: DistanceInput = new DistanceInput();
  const timeOfImpact_s_distanceOutput: DistanceOutput = new DistanceOutput();
  const timeOfImpact_s_fcn: SeparationFunction = new SeparationFunction();
  const timeOfImpact_s_indexA: [number] = [ 0 ];
  const timeOfImpact_s_indexB: [number] = [ 0 ];
  const timeOfImpact_s_sweepA: Sweep = new Sweep();
  const timeOfImpact_s_sweepB: Sweep = new Sweep();
  export function timeOfImpact(output: TOIOutput, input: TOIInput): void {
    const timer: Timer = timeOfImpact_s_timer.reset();

    ++toiCalls;

    output.state = TOIOutputState.Unknown;
    output.t = input.tMax;

    const proxyA: DistanceProxy = input.proxyA;
    const proxyB: DistanceProxy = input.proxyB;
    const maxVertices: number = Max(maxPolygonVertices, Max(proxyA.count, proxyB.count));

    const sweepA: Sweep = timeOfImpact_s_sweepA.copy(input.sweepA);
    const sweepB: Sweep = timeOfImpact_s_sweepB.copy(input.sweepB);

    // Large rotations can make the root finder fail, so we normalize the
    // sweep angles.
    sweepA.normalize();
    sweepB.normalize();

    const tMax: number = input.tMax;

    const totalRadius: number = proxyA.radius + proxyB.radius;
    const target: number = Max(linearSlop, totalRadius - 3 * linearSlop);
    const tolerance: number = 0.25 * linearSlop;
    // DEBUG: Assert(target > tolerance);

    let t1: number = 0;
    const k_maxIterations: number = 20; // TODO_ERIN Settings
    let iter: number = 0;

    // Prepare input for distance query.
    const cache: SimplexCache = timeOfImpact_s_cache;
    cache.count = 0;
    const distanceInput: DistanceInput = timeOfImpact_s_distanceInput;
    distanceInput.proxyA.copy(input.proxyA);
    distanceInput.proxyB.copy(input.proxyB);
    distanceInput.useRadii = false;

    // The outer loop progressively attempts to compute new separating axes.
    // This loop terminates when an axis is repeated (no progress is made).
    for (; ; ) {
      const xfA: Transform = timeOfImpact_s_xfA;
      const xfB: Transform = timeOfImpact_s_xfB;
      sweepA.getTransform(xfA, t1);
      sweepB.getTransform(xfB, t1);

      // Get the distance between shapes. We can also use the results
      // to get a separating axis.
      distanceInput.transformA.copy(xfA);
      distanceInput.transformB.copy(xfB);
      const distanceOutput: DistanceOutput = timeOfImpact_s_distanceOutput;
      distance(distanceOutput, cache, distanceInput);

      // If the shapes are overlapped, we give up on continuous collision.
      if (distanceOutput.distance <= 0) {
        // Failure!
        output.state = TOIOutputState.Overlapped;
        output.t = 0;
        break;
      }

      if (distanceOutput.distance < target + tolerance) {
        // Victory!
        output.state = TOIOutputState.Touching;
        output.t = t1;
        break;
      }

      // Initialize the separating axis.
      const fcn: SeparationFunction = timeOfImpact_s_fcn;
      fcn.initialize(cache, proxyA, sweepA, proxyB, sweepB, t1);
      /*
      #if 0
          // Dump the curve seen by the root finder {
            const int32 N = 100;
            float32 dx = 1.0f / N;
            float32 xs[N+1];
            float32 fs[N+1];

            float32 x = 0.0f;

            for (int32 i = 0; i <= N; ++i) {
              sweepA.GetTransform(&xfA, x);
              sweepB.GetTransform(&xfB, x);
              float32 f = fcn.Evaluate(xfA, xfB) - target;

              printf("%g %g\n", x, f);

              xs[i] = x;
              fs[i] = f;

              x += dx;
            }
          }
      #endif
      */

      // Compute the TOI on the separating axis. We do this by successively
      // resolving the deepest point. This loop is bounded by the number of vertices.
      let done: boolean = false;
      let t2: number = tMax;
      let pushBackIter: number = 0;
      for (; ; ) {
        // Find the deepest point at t2. Store the witness point indices.
        const indexA: [number] = timeOfImpact_s_indexA;
        const indexB: [number] = timeOfImpact_s_indexB;
        let s2: number = fcn.findMinSeparation(indexA, indexB, t2);

        // Is the final configuration separated?
        if (s2 > (target + tolerance)) {
          // Victory!
          output.state = TOIOutputState.Separated;
          output.t = tMax;
          done = true;
          break;
        }

        // Has the separation reached tolerance?
        if (s2 > (target - tolerance)) {
          // Advance the sweeps
          t1 = t2;
          break;
        }

        // Compute the initial separation of the witness points.
        let s1: number = fcn.evaluate(indexA[0], indexB[0], t1);

        // Check for initial overlap. This might happen if the root finder
        // runs out of iterations.
        if (s1 < (target - tolerance)) {
          output.state = TOIOutputState.Failed;
          output.t = t1;
          done = true;
          break;
        }

        // Check for touching
        if (s1 <= (target + tolerance)) {
          // Victory! t1 should hold the TOI (could be 0.0).
          output.state = TOIOutputState.Touching;
          output.t = t1;
          done = true;
          break;
        }

        // Compute 1D root of: f(x) - target = 0
        let rootIterCount: number = 0;
        let a1: number = t1;
        let a2: number = t2;
        for (; ; ) {
          // Use a mix of the secant rule and bisection.
          let t: number = 0;
          if (rootIterCount & 1) {
            // Secant rule to improve convergence.
            t = a1 + (target - s1) * (a2 - a1) / (s2 - s1);
          } else {
            // Bisection to guarantee progress.
            t = 0.5 * (a1 + a2);
          }

          ++rootIterCount;
          ++toiRootIters;

          const s: number = fcn.evaluate(indexA[0], indexB[0], t);

          if (Abs(s - target) < tolerance) {
            // t2 holds a tentative value for t1
            t2 = t;
            break;
          }

          // Ensure we continue to bracket the root.
          if (s > target) {
            a1 = t;
            s1 = s;
          } else {
            a2 = t;
            s2 = s;
          }

          if (rootIterCount === 50) {
            break;
          }
        }

        toiMaxRootIters = Max(toiMaxRootIters, rootIterCount);

        ++pushBackIter;

        if (pushBackIter === maxVertices) {
          break;
        }
      }

      ++iter;
      ++toiIters;

      if (done) {
        break;
      }

      if (iter === k_maxIterations) {
        // Root finder got stuck. Semi-victory.
        output.state = TOIOutputState.Failed;
        output.t = t1;
        break;
      }
    }

    toiMaxIters = Max(toiMaxIters, iter);

    const time: number = timer.getMilliseconds();
    toiMaxTime = Max(toiMaxTime, time);
    toiTime += time;
  }

}
