/*
* Copyright (c) 2011 Erin Catto http://box2d.org
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
  /// Timer for profiling. This has platform specific code and may
/// not work on every platform.
  export class Timer {
    public start: number = Date.now();

    /// Reset the timer.
    public reset(): Timer {
      this.start = Date.now();
      return this;
    }

    /// Get the time since construction or the last reset.
    public getMilliseconds(): number {
      return Date.now() - this.start;
    }
  }

  export class Counter {
    public count: number = 0;
    public minCount: number = 0;
    public maxCount: number = 0;

    public increment(): void {
      this.count++;
      if (this.maxCount < this.count) {
        this.maxCount = this.count;
      }
    }

    public decrement(): void {
      this.count--;
      if (this.minCount > this.count) {
        this.minCount = this.count;
      }
    }
  }

}
