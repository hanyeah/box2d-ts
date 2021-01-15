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

  export const pi_over_180: number = pi / 180;
  export const _180_over_pi: number = 180 / pi;
  export const two_pi: number = 2 * pi;

  export function clamp(a: number, lo: number, hi: number): number {
    return (a < lo) ? (lo) : ((a > hi) ? (hi) : (a));
  }

  export function swap<T>(a: T[], b: T[]): void {
    // DEBUG: Assert(false);
    const tmp: T = a[0];
    a[0] = b[0];
    b[0] = tmp;
  }

/// This function is used to ensure that a floating point number is
/// not a NaN or infinity.
  export const isValid = isFinite;

  export function sqrt(n: number): number {
    return n * n;
  }

/// This is a approximate yet fast inverse square-root.
  export function invSqrt(n: number): number {
    return 1 / Math.sqrt(n);
  }

  export function degToRad(degrees: number): number {
    return degrees * pi_over_180;
  }

  export function radToDeg(radians: number): number {
    return radians * _180_over_pi;
  }

  export const Pow = Math.pow;
  export const Sqrt = Math.sqrt;
  export const Abs = Math.abs;
  export const Cos = Math.cos;
  export const Sin = Math.sin;
  export const Acos = Math.acos;
  export const Asin = Math.asin;
  export const Atan2 = Math.atan2;
  export const Random = Math.random;
  export const Min = Math.min;
  export const Max = Math.max;

  export function nextPowerOfTwo(x: number): number {
    x |= (x >> 1) & 0x7FFFFFFF;
    x |= (x >> 2) & 0x3FFFFFFF;
    x |= (x >> 4) & 0x0FFFFFFF;
    x |= (x >> 8) & 0x00FFFFFF;
    x |= (x >> 16) & 0x0000FFFF;
    return x + 1;
  }

  export function isPowerOfTwo(x: number): boolean {
    return x > 0 && (x & (x - 1)) === 0;
  }

  export function randomRange(lo: number, hi: number): number {
    return (hi - lo) * Math.random() + lo;
  }

  export interface XY {
    x: number;
    y: number;
  }

/// A 2D column vector.
  export class Vec2 implements XY {
    public static readonly ZERO: Vec2 = new Vec2(0, 0);
    public static readonly UNITX: Vec2 = new Vec2(1, 0);
    public static readonly UNITY: Vec2 = new Vec2(0, 1);

    public static readonly s_t0: Vec2 = new Vec2();
    public static readonly s_t1: Vec2 = new Vec2();
    public static readonly s_t2: Vec2 = new Vec2();
    public static readonly s_t3: Vec2 = new Vec2();

    public constructor(public x: number = 0, public y: number = 0) {}

    public clone(): Vec2 {
      return new Vec2(this.x, this.y);
    }

    public setZero(): this {
      this.x = 0;
      this.y = 0;
      return this;
    }

    public set(x: number, y: number): this {
      this.x = x;
      this.y = y;
      return this;
    }

    public copy(other: XY): this {
      this.x = other.x;
      this.y = other.y;
      return this;
    }

    public selfAdd(v: XY): this {
      this.x += v.x;
      this.y += v.y;
      return this;
    }

    public selfAddXY(x: number, y: number): this {
      this.x += x;
      this.y += y;
      return this;
    }

    public selfSub(v: XY): this {
      this.x -= v.x;
      this.y -= v.y;
      return this;
    }

    public selfSubXY(x: number, y: number): this {
      this.x -= x;
      this.y -= y;
      return this;
    }

    public selfMul(s: number): this {
      this.x *= s;
      this.y *= s;
      return this;
    }

    public selfMulAdd(s: number, v: XY): this {
      this.x += s * v.x;
      this.y += s * v.y;
      return this;
    }

    public selfMulSub(s: number, v: XY): this {
      this.x -= s * v.x;
      this.y -= s * v.y;
      return this;
    }

    public dot(v: XY): number {
      return this.x * v.x + this.y * v.y;
    }

    public cross(v: XY): number {
      return this.x * v.y - this.y * v.x;
    }

    public length(): number {
      const x: number = this.x, y: number = this.y;
      return Math.sqrt(x * x + y * y);
    }

    public lengthSquared(): number {
      const x: number = this.x, y: number = this.y;
      return (x * x + y * y);
    }

    public normalize(): number {
      const length: number = this.length();
      if (length >= epsilon) {
        const inv_length: number = 1 / length;
        this.x *= inv_length;
        this.y *= inv_length;
      }
      return length;
    }

    public selfNormalize(): this {
      this.normalize();
      return this;
    }

    public selfRotate(radians: number): this {
      const c: number = Math.cos(radians);
      const s: number = Math.sin(radians);
      const x: number = this.x;
      this.x = c * x - s * this.y;
      this.y = s * x + c * this.y;
      return this;
    }

    public selfRotateCosSin(c: number, s: number): this {
      const x: number = this.x;
      this.x = c * x - s * this.y;
      this.y = s * x + c * this.y;
      return this;
    }

    public isValid(): boolean {
      return isFinite(this.x) && isFinite(this.y);
    }

    public selfCrossVS(s: number): this {
      const x: number = this.x;
      this.x =  s * this.y;
      this.y = -s * x;
      return this;
    }

    public selfCrossSV(s: number): this {
      const x: number = this.x;
      this.x = -s * this.y;
      this.y =  s * x;
      return this;
    }

    public selfMinV(v: XY): this {
      this.x = Min(this.x, v.x);
      this.y = Min(this.y, v.y);
      return this;
    }

    public selfMaxV(v: XY): this {
      this.x = Max(this.x, v.x);
      this.y = Max(this.y, v.y);
      return this;
    }

    public selfAbs(): this {
      this.x = Abs(this.x);
      this.y = Abs(this.y);
      return this;
    }

    public selfNeg(): this {
      this.x = (-this.x);
      this.y = (-this.y);
      return this;
    }

    public selfSkew(): this {
      const x: number = this.x;
      this.x = -this.y;
      this.y = x;
      return this;
    }

    public static MakeArray(length: number): Vec2[] {
      return MakeArray(length, (i: number): Vec2 => new Vec2());
    }

    public static AbsV<T extends XY>(v: XY, out: T): T {
      out.x = Abs(v.x);
      out.y = Abs(v.y);
      return out;
    }

    public static MinV<T extends XY>(a: XY, b: XY, out: T): T {
      out.x = Min(a.x, b.x);
      out.y = Min(a.y, b.y);
      return out;
    }

    public static MaxV<T extends XY>(a: XY, b: XY, out: T): T {
      out.x = Max(a.x, b.x);
      out.y = Max(a.y, b.y);
      return out;
    }

    public static ClampV<T extends XY>(v: XY, lo: XY, hi: XY, out: T): T {
      out.x = clamp(v.x, lo.x, hi.x);
      out.y = clamp(v.y, lo.y, hi.y);
      return out;
    }

    public static RotateV<T extends XY>(v: XY, radians: number, out: T): T {
      const v_x: number = v.x, v_y: number = v.y;
      const c: number = Math.cos(radians);
      const s: number = Math.sin(radians);
      out.x = c * v_x - s * v_y;
      out.y = s * v_x + c * v_y;
      return out;
    }

    public static DotVV(a: XY, b: XY): number {
      return a.x * b.x + a.y * b.y;
    }

    public static CrossVV(a: XY, b: XY): number {
      return a.x * b.y - a.y * b.x;
    }

    public static CrossVS<T extends XY>(v: XY, s: number, out: T): T {
      const v_x: number = v.x;
      out.x =  s * v.y;
      out.y = -s * v_x;
      return out;
    }

    public static CrossVOne<T extends XY>(v: XY, out: T): T {
      const v_x: number = v.x;
      out.x =  v.y;
      out.y = -v_x;
      return out;
    }

    public static CrossSV<T extends XY>(s: number, v: XY, out: T): T {
      const v_x: number = v.x;
      out.x = -s * v.y;
      out.y =  s * v_x;
      return out;
    }

    public static CrossOneV<T extends XY>(v: XY, out: T): T {
      const v_x: number = v.x;
      out.x = -v.y;
      out.y =  v_x;
      return out;
    }

    public static AddVV<T extends XY>(a: XY, b: XY, out: T): T { out.x = a.x + b.x; out.y = a.y + b.y; return out; }

    public static SubVV<T extends XY>(a: XY, b: XY, out: T): T { out.x = a.x - b.x; out.y = a.y - b.y; return out; }

    public static MulSV<T extends XY>(s: number, v: XY, out: T): T { out.x = v.x * s; out.y = v.y * s; return out; }
    public static MulVS<T extends XY>(v: XY, s: number, out: T): T { out.x = v.x * s; out.y = v.y * s; return out; }

    public static AddVMulSV<T extends XY>(a: XY, s: number, b: XY, out: T): T { out.x = a.x + (s * b.x); out.y = a.y + (s * b.y); return out; }
    public static SubVMulSV<T extends XY>(a: XY, s: number, b: XY, out: T): T { out.x = a.x - (s * b.x); out.y = a.y - (s * b.y); return out; }

    public static AddVCrossSV<T extends XY>(a: XY, s: number, v: XY, out: T): T {
      const v_x: number = v.x;
      out.x = a.x - (s * v.y);
      out.y = a.y + (s * v_x);
      return out;
    }

    public static MidVV<T extends XY>(a: XY, b: XY, out: T): T { out.x = (a.x + b.x) * 0.5; out.y = (a.y + b.y) * 0.5; return out; }

    public static ExtVV<T extends XY>(a: XY, b: XY, out: T): T { out.x = (b.x - a.x) * 0.5; out.y = (b.y - a.y) * 0.5; return out; }

    public static IsEqualToV(a: XY, b: XY): boolean {
      return a.x === b.x && a.y === b.y;
    }

    public static DistanceVV(a: XY, b: XY): number {
      const c_x: number = a.x - b.x;
      const c_y: number = a.y - b.y;
      return Math.sqrt(c_x * c_x + c_y * c_y);
    }

    public static DistanceSquaredVV(a: XY, b: XY): number {
      const c_x: number = a.x - b.x;
      const c_y: number = a.y - b.y;
      return (c_x * c_x + c_y * c_y);
    }

    public static NegV<T extends XY>(v: XY, out: T): T { out.x = -v.x; out.y = -v.y; return out; }

  }

  export const Vec2_zero: Vec2 = new Vec2(0, 0);

  export class TypedVec2 implements Vec2 {
    public readonly data: Float32Array;
    public get x(): number { return this.data[0]; } public set x(value: number) { this.data[0] = value; }
    public get y(): number { return this.data[1]; } public set y(value: number) { this.data[1] = value; }

    constructor();
    constructor(data: Float32Array);
    constructor(x: number, y: number);
    constructor(...args: any[]) {
      if (args[0] instanceof Float32Array) {
        if (args[0].length !== 2) { throw new Error(); }
        this.data = args[0];
      } else {
        const x: number = typeof args[0] === "number" ? args[0] : 0;
        const y: number = typeof args[1] === "number" ? args[1] : 0;
        this.data = new Float32Array([ x, y ]);
      }
    }

    public clone(): TypedVec2 {
      return new TypedVec2(new Float32Array(this.data));
    }

    public setZero(): this {
      this.x = 0;
      this.y = 0;
      return this;
    }

    public set(x: number, y: number): this {
      this.x = x;
      this.y = y;
      return this;
    }

    public copy(other: XY): this {
      if (other instanceof TypedVec2) {
        this.data.set(other.data);
      }
      else {
        this.x = other.x;
        this.y = other.y;
      }
      return this;
    }

    public selfAdd(v: XY): this {
      this.x += v.x;
      this.y += v.y;
      return this;
    }

    public selfAddXY(x: number, y: number): this {
      this.x += x;
      this.y += y;
      return this;
    }

    public selfSub(v: XY): this {
      this.x -= v.x;
      this.y -= v.y;
      return this;
    }

    public selfSubXY(x: number, y: number): this {
      this.x -= x;
      this.y -= y;
      return this;
    }

    public selfMul(s: number): this {
      this.x *= s;
      this.y *= s;
      return this;
    }

    public selfMulAdd(s: number, v: XY): this {
      this.x += s * v.x;
      this.y += s * v.y;
      return this;
    }

    public selfMulSub(s: number, v: XY): this {
      this.x -= s * v.x;
      this.y -= s * v.y;
      return this;
    }

    public dot(v: XY): number {
      return this.x * v.x + this.y * v.y;
    }

    public cross(v: XY): number {
      return this.x * v.y - this.y * v.x;
    }

    public length(): number {
      const x: number = this.x, y: number = this.y;
      return Math.sqrt(x * x + y * y);
    }

    public lengthSquared(): number {
      const x: number = this.x, y: number = this.y;
      return (x * x + y * y);
    }

    public normalize(): number {
      const length: number = this.length();
      if (length >= epsilon) {
        const inv_length: number = 1 / length;
        this.x *= inv_length;
        this.y *= inv_length;
      }
      return length;
    }

    public selfNormalize(): this {
      const length: number = this.length();
      if (length >= epsilon) {
        const inv_length: number = 1 / length;
        this.x *= inv_length;
        this.y *= inv_length;
      }
      return this;
    }

    public selfRotate(radians: number): this {
      const c: number = Math.cos(radians);
      const s: number = Math.sin(radians);
      const x: number = this.x;
      this.x = c * x - s * this.y;
      this.y = s * x + c * this.y;
      return this;
    }

    public selfRotateCosSin(c: number, s: number): this {
      const x: number = this.x;
      this.x = c * x - s * this.y;
      this.y = s * x + c * this.y;
      return this;
    }

    public isValid(): boolean {
      return isFinite(this.x) && isFinite(this.y);
    }

    public selfCrossVS(s: number): this {
      const x: number = this.x;
      this.x =  s * this.y;
      this.y = -s * x;
      return this;
    }

    public selfCrossSV(s: number): this {
      const x: number = this.x;
      this.x = -s * this.y;
      this.y =  s * x;
      return this;
    }

    public selfMinV(v: XY): this {
      this.x = Min(this.x, v.x);
      this.y = Min(this.y, v.y);
      return this;
    }

    public selfMaxV(v: XY): this {
      this.x = Max(this.x, v.x);
      this.y = Max(this.y, v.y);
      return this;
    }

    public selfAbs(): this {
      this.x = Abs(this.x);
      this.y = Abs(this.y);
      return this;
    }

    public selfNeg(): this {
      this.x = (-this.x);
      this.y = (-this.y);
      return this;
    }

    public selfSkew(): this {
      const x: number = this.x;
      this.x = -this.y;
      this.y = x;
      return this;
    }
  }

  export interface XYZ extends XY {
    z: number;
  }

/// A 2D column vector with 3 elements.
  export class Vec3 implements XYZ {
    public static readonly ZERO: Vec3 = new Vec3(0, 0, 0);

    public static readonly s_t0: Vec3 = new Vec3();

    public readonly data: Float32Array;
    public get x(): number { return this.data[0]; } public set x(value: number) { this.data[0] = value; }
    public get y(): number { return this.data[1]; } public set y(value: number) { this.data[1] = value; }
    public get z(): number { return this.data[2]; } public set z(value: number) { this.data[2] = value; }

    constructor();
    constructor(data: Float32Array);
    constructor(x: number, y: number, z: number);
    constructor(...args: any[]) {
      if (args[0] instanceof Float32Array) {
        if (args[0].length !== 3) { throw new Error(); }
        this.data = args[0];
      } else {
        const x: number = typeof args[0] === "number" ? args[0] : 0;
        const y: number = typeof args[1] === "number" ? args[1] : 0;
        const z: number = typeof args[2] === "number" ? args[2] : 0;
        this.data = new Float32Array([ x, y, z ]);
      }
    }

    public clone(): Vec3 {
      return new Vec3(this.x, this.y, this.z);
    }

    public setZero(): this {
      this.x = 0;
      this.y = 0;
      this.z = 0;
      return this;
    }

    public setXYZ(x: number, y: number, z: number): this {
      this.x = x;
      this.y = y;
      this.z = z;
      return this;
    }

    public copy(other: XYZ): this {
      this.x = other.x;
      this.y = other.y;
      this.z = other.z;
      return this;
    }

    public selfNeg(): this {
      this.x = (-this.x);
      this.y = (-this.y);
      this.z = (-this.z);
      return this;
    }

    public selfAdd(v: XYZ): this {
      this.x += v.x;
      this.y += v.y;
      this.z += v.z;
      return this;
    }

    public selfAddXYZ(x: number, y: number, z: number): this {
      this.x += x;
      this.y += y;
      this.z += z;
      return this;
    }

    public selfSub(v: XYZ): this {
      this.x -= v.x;
      this.y -= v.y;
      this.z -= v.z;
      return this;
    }

    public selfSubXYZ(x: number, y: number, z: number): this {
      this.x -= x;
      this.y -= y;
      this.z -= z;
      return this;
    }

    public selfMul(s: number): this {
      this.x *= s;
      this.y *= s;
      this.z *= s;
      return this;
    }

    public static dotV3V3(a: XYZ, b: XYZ): number {
      return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    public static crossV3V3<T extends XYZ>(a: XYZ, b: XYZ, out: T): T {
      const a_x: number = a.x, a_y = a.y, a_z = a.z;
      const b_x: number = b.x, b_y = b.y, b_z = b.z;
      out.x = a_y * b_z - a_z * b_y;
      out.y = a_z * b_x - a_x * b_z;
      out.z = a_x * b_y - a_y * b_x;
      return out;
    }
  }

/// A 2-by-2 matrix. Stored in column-major order.
  export class Mat22 {
    public static readonly IDENTITY: Mat22 = new Mat22();

    // public readonly data: Float32Array = new Float32Array([ 1, 0, 0, 1 ]);
    // public readonly ex: Vec2 = new Vec2(this.data.subarray(0, 2));
    // public readonly ey: Vec2 = new Vec2(this.data.subarray(2, 4));

    public readonly ex: Vec2 = new Vec2(1, 0);
    public readonly ey: Vec2 = new Vec2(0, 1);

    public clone(): Mat22 {
      return new Mat22().copy(this);
    }

    public static fromVV(c1: XY, c2: XY): Mat22 {
      return new Mat22().setVV(c1, c2);
    }

    public static fromSSSS(r1c1: number, r1c2: number, r2c1: number, r2c2: number): Mat22 {
      return new Mat22().setSSSS(r1c1, r1c2, r2c1, r2c2);
    }

    public static fromAngle(radians: number): Mat22 {
      return new Mat22().setAngle(radians);
    }

    public setSSSS(r1c1: number, r1c2: number, r2c1: number, r2c2: number): this {
      this.ex.set(r1c1, r2c1);
      this.ey.set(r1c2, r2c2);
      return this;
    }

    public setVV(c1: XY, c2: XY): this {
      this.ex.copy(c1);
      this.ey.copy(c2);
      return this;
    }

    public setAngle(radians: number): this {
      const c: number = Math.cos(radians);
      const s: number = Math.sin(radians);
      this.ex.set( c, s);
      this.ey.set(-s, c);
      return this;
    }

    public copy(other: Mat22): this {
      this.ex.copy(other.ex);
      this.ey.copy(other.ey);
      return this;
    }

    public setIdentity(): this {
      this.ex.set(1, 0);
      this.ey.set(0, 1);
      return this;
    }

    public setZero(): this {
      this.ex.setZero();
      this.ey.setZero();
      return this;
    }

    public getAngle(): number {
      return Math.atan2(this.ex.y, this.ex.x);
    }

    public getInverse(out: Mat22): Mat22 {
      const a: number = this.ex.x;
      const b: number = this.ey.x;
      const c: number = this.ex.y;
      const d: number = this.ey.y;
      let det: number = a * d - b * c;
      if (det !== 0) {
        det = 1 / det;
      }
      out.ex.x =   det * d;
      out.ey.x = (-det * b);
      out.ex.y = (-det * c);
      out.ey.y =   det * a;
      return out;
    }

    public solve<T extends XY>(b_x: number, b_y: number, out: T): T {
      const a11: number = this.ex.x, a12 = this.ey.x;
      const a21: number = this.ex.y, a22 = this.ey.y;
      let det: number = a11 * a22 - a12 * a21;
      if (det !== 0) {
        det = 1 / det;
      }
      out.x = det * (a22 * b_x - a12 * b_y);
      out.y = det * (a11 * b_y - a21 * b_x);
      return out;
    }

    public selfAbs(): this {
      this.ex.selfAbs();
      this.ey.selfAbs();
      return this;
    }

    public selfInv(): this {
      this.getInverse(this);
      return this;
    }

    public selfAddM(M: Mat22): this {
      this.ex.selfAdd(M.ex);
      this.ey.selfAdd(M.ey);
      return this;
    }

    public selfSubM(M: Mat22): this {
      this.ex.selfSub(M.ex);
      this.ey.selfSub(M.ey);
      return this;
    }

    public static absM(M: Mat22, out: Mat22): Mat22 {
      const ex: Vec2 = M.ex, ey: Vec2 = M.ey;
      out.ex.x = Abs(ex.x);
      out.ex.y = Abs(ex.y);
      out.ey.x = Abs(ey.x);
      out.ey.y = Abs(ey.y);
      return out;
    }

    public static mulMV<T extends XY>(M: Mat22, v: XY, out: T): T {
      const ex: Vec2 = M.ex, ey: Vec2 = M.ey;
      const v_x: number = v.x, v_y: number = v.y;
      out.x = ex.x * v_x + ey.x * v_y;
      out.y = ex.y * v_x + ey.y * v_y;
      return out;
    }

    public static mulTMV<T extends XY>(M: Mat22, v: XY, out: T): T {
      const ex: Vec2 = M.ex, ey: Vec2 = M.ey;
      const v_x: number = v.x, v_y: number = v.y;
      out.x = ex.x * v_x + ex.y * v_y;
      out.y = ey.x * v_x + ey.y * v_y;
      return out;
    }

    public static addMM(A: Mat22, B: Mat22, out: Mat22): Mat22 {
      const A_ex: Vec2 = A.ex, A_ey: Vec2 = A.ey;
      const B_ex: Vec2 = B.ex, B_ey: Vec2 = B.ey;
      out.ex.x = A_ex.x + B_ex.x;
      out.ex.y = A_ex.y + B_ex.y;
      out.ey.x = A_ey.x + B_ey.x;
      out.ey.y = A_ey.y + B_ey.y;
      return out;
    }

    public static mulMM(A: Mat22, B: Mat22, out: Mat22): Mat22 {
      const A_ex_x: number = A.ex.x, A_ex_y: number = A.ex.y;
      const A_ey_x: number = A.ey.x, A_ey_y: number = A.ey.y;
      const B_ex_x: number = B.ex.x, B_ex_y: number = B.ex.y;
      const B_ey_x: number = B.ey.x, B_ey_y: number = B.ey.y;
      out.ex.x = A_ex_x * B_ex_x + A_ey_x * B_ex_y;
      out.ex.y = A_ex_y * B_ex_x + A_ey_y * B_ex_y;
      out.ey.x = A_ex_x * B_ey_x + A_ey_x * B_ey_y;
      out.ey.y = A_ex_y * B_ey_x + A_ey_y * B_ey_y;
      return out;
    }

    public static mulTMM(A: Mat22, B: Mat22, out: Mat22): Mat22 {
      const A_ex_x: number = A.ex.x, A_ex_y: number = A.ex.y;
      const A_ey_x: number = A.ey.x, A_ey_y: number = A.ey.y;
      const B_ex_x: number = B.ex.x, B_ex_y: number = B.ex.y;
      const B_ey_x: number = B.ey.x, B_ey_y: number = B.ey.y;
      out.ex.x = A_ex_x * B_ex_x + A_ex_y * B_ex_y;
      out.ex.y = A_ey_x * B_ex_x + A_ey_y * B_ex_y;
      out.ey.x = A_ex_x * B_ey_x + A_ex_y * B_ey_y;
      out.ey.y = A_ey_x * B_ey_x + A_ey_y * B_ey_y;
      return out;
    }
  }

/// A 3-by-3 matrix. Stored in column-major order.
  export class Mat33 {
    public static readonly IDENTITY: Mat33 = new Mat33();

    public readonly data: Float32Array = new Float32Array([ 1, 0, 0, 0, 1, 0, 0, 0, 1 ]);
    public readonly ex: Vec3 = new Vec3(this.data.subarray(0, 3));
    public readonly ey: Vec3 = new Vec3(this.data.subarray(3, 6));
    public readonly ez: Vec3 = new Vec3(this.data.subarray(6, 9));

    public clone(): Mat33 {
      return new Mat33().copy(this);
    }

    public setVVV(c1: XYZ, c2: XYZ, c3: XYZ): this {
      this.ex.copy(c1);
      this.ey.copy(c2);
      this.ez.copy(c3);
      return this;
    }

    public copy(other: Mat33): this {
      this.ex.copy(other.ex);
      this.ey.copy(other.ey);
      this.ez.copy(other.ez);
      return this;
    }

    public setIdentity(): this {
      this.ex.setXYZ(1, 0, 0);
      this.ey.setXYZ(0, 1, 0);
      this.ez.setXYZ(0, 0, 1);
      return this;
    }

    public setZero(): this {
      this.ex.setZero();
      this.ey.setZero();
      this.ez.setZero();
      return this;
    }

    public selfAddM(M: Mat33): this {
      this.ex.selfAdd(M.ex);
      this.ey.selfAdd(M.ey);
      this.ez.selfAdd(M.ez);
      return this;
    }

    public solve33<T extends XYZ>(b_x: number, b_y: number, b_z: number, out: T): T {
      const a11: number = this.ex.x, a21: number = this.ex.y, a31: number = this.ex.z;
      const a12: number = this.ey.x, a22: number = this.ey.y, a32: number = this.ey.z;
      const a13: number = this.ez.x, a23: number = this.ez.y, a33: number = this.ez.z;
      let det: number = a11 * (a22 * a33 - a32 * a23) + a21 * (a32 * a13 - a12 * a33) + a31 * (a12 * a23 - a22 * a13);
      if (det !== 0) {
        det = 1 / det;
      }
      out.x = det * (b_x * (a22 * a33 - a32 * a23) + b_y * (a32 * a13 - a12 * a33) + b_z * (a12 * a23 - a22 * a13));
      out.y = det * (a11 * (b_y * a33 - b_z * a23) + a21 * (b_z * a13 - b_x * a33) + a31 * (b_x * a23 - b_y * a13));
      out.z = det * (a11 * (a22 * b_z - a32 * b_y) + a21 * (a32 * b_x - a12 * b_z) + a31 * (a12 * b_y - a22 * b_x));
      return out;
    }

    public solve22<T extends XY>(b_x: number, b_y: number, out: T): T {
      const a11: number = this.ex.x, a12: number = this.ey.x;
      const a21: number = this.ex.y, a22: number = this.ey.y;
      let det: number = a11 * a22 - a12 * a21;
      if (det !== 0) {
        det = 1 / det;
      }
      out.x = det * (a22 * b_x - a12 * b_y);
      out.y = det * (a11 * b_y - a21 * b_x);
      return out;
    }

    public getInverse22(M: Mat33): void {
      const a: number = this.ex.x, b: number = this.ey.x, c: number = this.ex.y, d: number = this.ey.y;
      let det: number = a * d - b * c;
      if (det !== 0) {
        det = 1 / det;
      }

      M.ex.x =  det * d; M.ey.x = -det * b; M.ex.z = 0;
      M.ex.y = -det * c; M.ey.y =  det * a; M.ey.z = 0;
      M.ez.x =        0; M.ez.y =        0; M.ez.z = 0;
    }

    public getSymInverse33(M: Mat33): void {
      let det: number = Vec3.dotV3V3(this.ex, Vec3.crossV3V3(this.ey, this.ez, Vec3.s_t0));
      if (det !== 0) {
        det = 1 / det;
      }

      const a11: number = this.ex.x, a12: number = this.ey.x, a13: number = this.ez.x;
      const a22: number = this.ey.y, a23: number = this.ez.y;
      const a33: number = this.ez.z;

      M.ex.x = det * (a22 * a33 - a23 * a23);
      M.ex.y = det * (a13 * a23 - a12 * a33);
      M.ex.z = det * (a12 * a23 - a13 * a22);

      M.ey.x = M.ex.y;
      M.ey.y = det * (a11 * a33 - a13 * a13);
      M.ey.z = det * (a13 * a12 - a11 * a23);

      M.ez.x = M.ex.z;
      M.ez.y = M.ey.z;
      M.ez.z = det * (a11 * a22 - a12 * a12);
    }

    public static mulM33V3<T extends XYZ>(A: Mat33, v: XYZ, out: T): T {
      const v_x: number = v.x, v_y: number = v.y, v_z: number = v.z;
      out.x = A.ex.x * v_x + A.ey.x * v_y + A.ez.x * v_z;
      out.y = A.ex.y * v_x + A.ey.y * v_y + A.ez.y * v_z;
      out.z = A.ex.z * v_x + A.ey.z * v_y + A.ez.z * v_z;
      return out;
    }
    public static mulM33XYZ<T extends XYZ>(A: Mat33, x: number, y: number, z: number, out: T): T {
      out.x = A.ex.x * x + A.ey.x * y + A.ez.x * z;
      out.y = A.ex.y * x + A.ey.y * y + A.ez.y * z;
      out.z = A.ex.z * x + A.ey.z * y + A.ez.z * z;
      return out;
    }
    public static mulM33V2<T extends XY>(A: Mat33, v: XY, out: T): T {
      const v_x: number = v.x, v_y: number = v.y;
      out.x = A.ex.x * v_x + A.ey.x * v_y;
      out.y = A.ex.y * v_x + A.ey.y * v_y;
      return out;
    }
    public static mulM33XY<T extends XY>(A: Mat33, x: number, y: number, out: T): T {
      out.x = A.ex.x * x + A.ey.x * y;
      out.y = A.ex.y * x + A.ey.y * y;
      return out;
    }
  }

/// Rotation
  export class Rot {
    public static readonly IDENTITY: Rot = new Rot();

    public s: number = 0;
    public c: number = 1;

    constructor(angle: number = 0) {
      if (angle) {
        this.s = Math.sin(angle);
        this.c = Math.cos(angle);
      }
    }

    public clone(): Rot {
      return new Rot().copy(this);
    }

    public copy(other: Rot): this {
      this.s = other.s;
      this.c = other.c;
      return this;
    }

    public setAngle(angle: number): this {
      this.s = Math.sin(angle);
      this.c = Math.cos(angle);
      return this;
    }

    public setIdentity(): this {
      this.s = 0;
      this.c = 1;
      return this;
    }

    public getAngle(): number {
      return Math.atan2(this.s, this.c);
    }

    public getXAxis<T extends XY>(out: T): T {
      out.x = this.c;
      out.y = this.s;
      return out;
    }

    public getYAxis<T extends XY>(out: T): T {
      out.x = -this.s;
      out.y = this.c;
      return out;
    }

    public static mulRR(q: Rot, r: Rot, out: Rot): Rot {
      // [qc -qs] * [rc -rs] = [qc*rc-qs*rs -qc*rs-qs*rc]
      // [qs  qc]   [rs  rc]   [qs*rc+qc*rs -qs*rs+qc*rc]
      // s = qs * rc + qc * rs
      // c = qc * rc - qs * rs
      const q_c: number = q.c, q_s: number = q.s;
      const r_c: number = r.c, r_s: number = r.s;
      out.s = q_s * r_c + q_c * r_s;
      out.c = q_c * r_c - q_s * r_s;
      return out;
    }

    public static mulTRR(q: Rot, r: Rot, out: Rot): Rot {
      // [ qc qs] * [rc -rs] = [qc*rc+qs*rs -qc*rs+qs*rc]
      // [-qs qc]   [rs  rc]   [-qs*rc+qc*rs qs*rs+qc*rc]
      // s = qc * rs - qs * rc
      // c = qc * rc + qs * rs
      const q_c: number = q.c, q_s: number = q.s;
      const r_c: number = r.c, r_s: number = r.s;
      out.s = q_c * r_s - q_s * r_c;
      out.c = q_c * r_c + q_s * r_s;
      return out;
    }

    public static mulRV<T extends XY>(q: Rot, v: XY, out: T): T {
      const q_c: number = q.c, q_s: number = q.s;
      const v_x: number = v.x, v_y: number = v.y;
      out.x = q_c * v_x - q_s * v_y;
      out.y = q_s * v_x + q_c * v_y;
      return out;
    }

    public static mulTRV<T extends XY>(q: Rot, v: XY, out: T): T {
      const q_c: number = q.c, q_s: number = q.s;
      const v_x: number = v.x, v_y: number = v.y;
      out.x =  q_c * v_x + q_s * v_y;
      out.y = -q_s * v_x + q_c * v_y;
      return out;
    }
  }

/// A transform contains translation and rotation. It is used to represent
/// the position and orientation of rigid frames.
  export class Transform {
    public static readonly IDENTITY: Transform = new Transform();

    public readonly p: Vec2 = new Vec2();
    public readonly q: Rot = new Rot();

    public clone(): Transform {
      return new Transform().copy(this);
    }

    public copy(other: Transform): this {
      this.p.copy(other.p);
      this.q.copy(other.q);
      return this;
    }

    public setIdentity(): this {
      this.p.setZero();
      this.q.setIdentity();
      return this;
    }

    public setPositionRotation(position: XY, q: Rot): this {
      this.p.copy(position);
      this.q.copy(q);
      return this;
    }

    public setPositionAngle(pos: XY, a: number): this {
      this.p.copy(pos);
      this.q.setAngle(a);
      return this;
    }

    public setPosition(position: XY): this {
      this.p.copy(position);
      return this;
    }

    public setPositionXY(x: number, y: number): this {
      this.p.set(x, y);
      return this;
    }

    public setRotation(rotation: Rot): this {
      this.q.copy(rotation);
      return this;
    }

    public setRotationAngle(radians: number): this {
      this.q.setAngle(radians);
      return this;
    }

    public getPosition(): Vec2 {
      return this.p;
    }

    public getRotation(): Rot {
      return this.q;
    }

    public getRotationAngle(): number {
      return this.q.getAngle();
    }

    public getAngle(): number {
      return this.q.getAngle();
    }

    public static mulXV<T extends XY>(T: Transform, v: XY, out: T): T {
      // float32 x = (T.q.c * v.x - T.q.s * v.y) + T.p.x;
      // float32 y = (T.q.s * v.x + T.q.c * v.y) + T.p.y;
      // return Vec2(x, y);
      const T_q_c: number = T.q.c, T_q_s: number = T.q.s;
      const v_x: number = v.x, v_y: number = v.y;
      out.x = (T_q_c * v_x - T_q_s * v_y) + T.p.x;
      out.y = (T_q_s * v_x + T_q_c * v_y) + T.p.y;
      return out;
    }

    public static mulTXV<T extends XY>(T: Transform, v: XY, out: T): T {
      // float32 px = v.x - T.p.x;
      // float32 py = v.y - T.p.y;
      // float32 x = (T.q.c * px + T.q.s * py);
      // float32 y = (-T.q.s * px + T.q.c * py);
      // return Vec2(x, y);
      const T_q_c: number = T.q.c, T_q_s: number = T.q.s;
      const p_x: number = v.x - T.p.x;
      const p_y: number = v.y - T.p.y;
      out.x = ( T_q_c * p_x + T_q_s * p_y);
      out.y = (-T_q_s * p_x + T_q_c * p_y);
      return out;
    }

    public static mulXX(A: Transform, B: Transform, out: Transform): Transform {
      Rot.mulRR(A.q, B.q, out.q);
      Vec2.AddVV(Rot.mulRV(A.q, B.p, out.p), A.p, out.p);
      return out;
    }

    public static mulTXX(A: Transform, B: Transform, out: Transform): Transform {
      Rot.mulTRR(A.q, B.q, out.q);
      Rot.mulTRV(A.q, Vec2.SubVV(B.p, A.p, out.p), out.p);
      return out;
    }

  }

/// This describes the motion of a body/shape for TOI computation.
/// Shapes are defined with respect to the body origin, which may
/// no coincide with the center of mass. However, to support dynamics
/// we must interpolate the center of mass position.
  export class Sweep {
    public readonly localCenter: Vec2 = new Vec2();
    public readonly c0: Vec2 = new Vec2();
    public readonly c: Vec2 = new Vec2();
    public a0: number = 0;
    public a: number = 0;
    public alpha0: number = 0;

    public clone(): Sweep {
      return new Sweep().copy(this);
    }

    public copy(other: Sweep): this {
      this.localCenter.copy(other.localCenter);
      this.c0.copy(other.c0);
      this.c.copy(other.c);
      this.a0 = other.a0;
      this.a = other.a;
      this.alpha0 = other.alpha0;
      return this;
    }

    // https://fgiesen.wordpress.com/2012/08/15/linear-interpolation-past-present-and-future/
    public getTransform(xf: Transform, beta: number): Transform {
      xf.p.x = (1.0 - beta) * this.c0.x + beta * this.c.x;
      xf.p.y = (1.0 - beta) * this.c0.y + beta * this.c.y;
      const angle: number = (1.0 - beta) * this.a0 + beta * this.a;
      xf.q.setAngle(angle);

      xf.p.selfSub(Rot.mulRV(xf.q, this.localCenter, Vec2.s_t0));
      return xf;
    }

    public advance(alpha: number): void {
      // DEBUG: Assert(this.alpha0 < 1);
      const beta: number = (alpha - this.alpha0) / (1 - this.alpha0);
      const one_minus_beta: number = (1 - beta);
      this.c0.x = one_minus_beta * this.c0.x + beta * this.c.x;
      this.c0.y = one_minus_beta * this.c0.y + beta * this.c.y;
      this.a0 = one_minus_beta * this.a0 + beta * this.a;
      this.alpha0 = alpha;
    }

    public normalize(): void {
      const d: number = two_pi * Math.floor(this.a0 / two_pi);
      this.a0 -= d;
      this.a -= d;
    }
  }

}
