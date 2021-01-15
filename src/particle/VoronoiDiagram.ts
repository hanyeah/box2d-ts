/*
 * Copyright (c) 2013 Google, Inc.
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
/**
 * A field representing the nearest generator from each point.
 */
namespace b2 {
  export class VoronoiDiagram {
    public generatorBuffer: VoronoiDiagraGenerator[];
    public generatorCapacity = 0;
    public generatorCount = 0;
    public countX = 0;
    public countY = 0;
    public diagram: VoronoiDiagraGenerator[] = [];

    constructor(generatorCapacity: number) {
      this.generatorBuffer = MakeArray(generatorCapacity, (index) => new VoronoiDiagraGenerator());
      this.generatorCapacity = generatorCapacity;
    }

    /**
     * Add a generator.
     *
     * @param center the position of the generator.
     * @param tag a tag used to identify the generator in callback functions.
     * @param necessary whether to callback for nodes associated with the generator.
     */
    public AddGenerator(center: Vec2, tag: number, necessary: boolean): void {
      // DEBUG: Assert(this.generatorCount < this.generatorCapacity);
      const g = this.generatorBuffer[this.generatorCount++];
      g.center.Copy(center);
      g.tag = tag;
      g.necessary = necessary;
    }

    /**
     * Generate the Voronoi diagram. It is rasterized with a given
     * interval in the same range as the necessary generators exist.
     *
     * @param radius the interval of the diagram.
     * @param margin margin for which the range of the diagram is extended.
     */
    public Generate(radius: number, margin: number): void {
      const inverseRadius = 1 / radius;
      const lower = new Vec2(+maxFloat, +maxFloat);
      const upper = new Vec2(-maxFloat, -maxFloat);
      let necessary_count = 0;
      for (let k = 0; k < this.generatorCount; k++) {
        const g = this.generatorBuffer[k];
        if (g.necessary) {
          Vec2.MinV(lower, g.center, lower);
          Vec2.MaxV(upper, g.center, upper);
          ++necessary_count;
        }
      }
      if (necessary_count === 0) {
        ///debugger;
        this.countX = 0;
        this.countY = 0;
        return;
      }
      lower.x -= margin;
      lower.y -= margin;
      upper.x += margin;
      upper.y += margin;
      this.countX = 1 + Math.floor(inverseRadius * (upper.x - lower.x));
      this.countY = 1 + Math.floor(inverseRadius * (upper.y - lower.y));
      this.diagram = []; // MakeArray(this.countX * this.countY, (index) => null);

      // (4 * countX * countY) is the queue capacity that is experimentally
      // known to be necessary and sufficient for general particle distributions.
      const queue = new StackQueue<VoronoiDiagraTask>(4 * this.countX * this.countY);
      for (let k = 0; k < this.generatorCount; k++) {
        const g = this.generatorBuffer[k];
        ///  g.center = inverseRadius * (g.center - lower);
        g.center.SelfSub(lower).SelfMul(inverseRadius);
        const x = Math.floor(g.center.x);
        const y = Math.floor(g.center.y);
        if (x >= 0 && y >= 0 && x < this.countX && y < this.countY) {
          queue.Push(new VoronoiDiagraTask(x, y, x + y * this.countX, g));
        }
      }
      while (!queue.Empty()) {
        const task = queue.Front();
        const x = task.x;
        const y = task.y;
        const i = task.i;
        const g = task.generator;
        queue.Pop();
        if (!this.diagram[i]) {
          this.diagram[i] = g;
          if (x > 0) {
            queue.Push(new VoronoiDiagraTask(x - 1, y, i - 1, g));
          }
          if (y > 0) {
            queue.Push(new VoronoiDiagraTask(x, y - 1, i - this.countX, g));
          }
          if (x < this.countX - 1) {
            queue.Push(new VoronoiDiagraTask(x + 1, y, i + 1, g));
          }
          if (y < this.countY - 1) {
            queue.Push(new VoronoiDiagraTask(x, y + 1, i + this.countX, g));
          }
        }
      }
      for (let y = 0; y < this.countY; y++) {
        for (let x = 0; x < this.countX - 1; x++) {
          const i = x + y * this.countX;
          const a = this.diagram[i];
          const b = this.diagram[i + 1];
          if (a !== b) {
            queue.Push(new VoronoiDiagraTask(x, y, i, b));
            queue.Push(new VoronoiDiagraTask(x + 1, y, i + 1, a));
          }
        }
      }
      for (let y = 0; y < this.countY - 1; y++) {
        for (let x = 0; x < this.countX; x++) {
          const i = x + y * this.countX;
          const a = this.diagram[i];
          const b = this.diagram[i + this.countX];
          if (a !== b) {
            queue.Push(new VoronoiDiagraTask(x, y, i, b));
            queue.Push(new VoronoiDiagraTask(x, y + 1, i + this.countX, a));
          }
        }
      }
      while (!queue.Empty()) {
        const task = queue.Front();
        const x = task.x;
        const y = task.y;
        const i = task.i;
        const k = task.generator;
        queue.Pop();
        const a = this.diagram[i];
        const b = k;
        if (a !== b) {
          const ax = a.center.x - x;
          const ay = a.center.y - y;
          const bx = b.center.x - x;
          const by = b.center.y - y;
          const a2 = ax * ax + ay * ay;
          const b2 = bx * bx + by * by;
          if (a2 > b2) {
            this.diagram[i] = b;
            if (x > 0) {
              queue.Push(new VoronoiDiagraTask(x - 1, y, i - 1, b));
            }
            if (y > 0) {
              queue.Push(new VoronoiDiagraTask(x, y - 1, i - this.countX, b));
            }
            if (x < this.countX - 1) {
              queue.Push(new VoronoiDiagraTask(x + 1, y, i + 1, b));
            }
            if (y < this.countY - 1) {
              queue.Push(new VoronoiDiagraTask(x, y + 1, i + this.countX, b));
            }
          }
        }
      }
    }

    /**
     * Enumerate all nodes that contain at least one necessary
     * generator.
     */
    public GetNodes(callback: VoronoiDiagraNodeCallback): void {
      for (let y = 0; y < this.countY - 1; y++) {
        for (let x = 0; x < this.countX - 1; x++) {
          const i = x + y * this.countX;
          const a = this.diagram[i];
          const b = this.diagram[i + 1];
          const c = this.diagram[i + this.countX];
          const d = this.diagram[i + 1 + this.countX];
          if (b !== c) {
            if (a !== b && a !== c &&
              (a.necessary || b.necessary || c.necessary)) {
              callback(a.tag, b.tag, c.tag);
            }
            if (d !== b && d !== c &&
              (a.necessary || b.necessary || c.necessary)) {
              callback(b.tag, d.tag, c.tag);
            }
          }
        }
      }
    }
  }

  /**
   * Callback used by GetNodes().
   *
   * Receive tags for generators associated with a node.
   */
  export type VoronoiDiagraNodeCallback = (a: number, b: number, c: number) => void;

  export class VoronoiDiagraGenerator {
    public center: Vec2 = new Vec2();
    public tag: number = 0;
    public necessary: boolean = false;
  }

  export class VoronoiDiagraTask {
    public x: number;
    public y: number;
    public i: number;
    public generator: VoronoiDiagraGenerator;
    constructor(x: number, y: number, i: number, g: VoronoiDiagraGenerator) {
      this.x = x;
      this.y = y;
      this.i = i;
      this.generator = g;
    }
  }

}

