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
    public Reset(): Timer {
      this.start = Date.now();
      return this;
    }

    /// Get the time since construction or the last reset.
    public GetMilliseconds(): number {
      return Date.now() - this.start;
    }
  }

  export class Counter {
    public count: number = 0;
    public min_count: number = 0;
    public max_count: number = 0;

    public GetCount(): number {
      return this.count;
    }

    public GetMinCount(): number {
      return this.min_count;
    }

    public GetMaxCount(): number {
      return this.max_count;
    }

    public ResetCount(): number {
      const count: number = this.count;
      this.count = 0;
      return count;
    }

    public ResetMinCount(): void {
      this.min_count = 0;
    }

    public ResetMaxCount(): void {
      this.max_count = 0;
    }

    public Increment(): void {
      this.count++;

      if (this.max_count < this.count) {
        this.max_count = this.count;
      }
    }

    public Decrement(): void {
      this.count--;

      if (this.min_count > this.count) {
        this.min_count = this.count;
      }
    }
  }

}
