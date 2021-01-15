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
namespace b2 {
  export class StackQueue<T> {
    public readonly buffer: Array<T> = [];
    public front: number = 0;
    public back: number = 0;
    public get capacity(): number { return this.buffer.length; }
    constructor(capacity: number) {
      for(let i: number = 0; i < capacity; i++) {
        this.buffer[i] = null;
      }
    }
    public push(item: T): void {
      if (this.back >= this.capacity) {
        for (let i = this.front; i < this.back; i++) {
          this.buffer[i - this.front] = this.buffer[i];
        }
        this.back -= this.front;
        this.front = 0;
      }
      this.buffer[this.back] = item;
      this.back++;
    }
    public pop(): void {
      // DEBUG: Assert(this.front < this.back);
      this.buffer[this.front] = null;
      this.front++;
    }
    public empty(): boolean {
      // DEBUG: Assert(this.front <= this.back);
      return this.front === this.back;
    }
    public getFront(): T {
      const item = this.buffer[this.front];
      if (!item) { throw new Error(); }
      return item;
    }
  }

}

