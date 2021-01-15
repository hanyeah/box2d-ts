declare namespace b2 {
    function Assert(condition: boolean, ...args: any[]): void;
    function Maybe<T>(value: T | undefined, def: T): T;
    const maxFloat: number;
    const epsilon: number;
    const epsilon_sq: number;
    const pi: number;
    const lengthUnitsPerMeter: number;
    const maxPolygonVertices: number;
    const maxManifoldPoints: number;
    const aabbExtension: number;
    const aabbMultiplier: number;
    const linearSlop: number;
    const angularSlop: number;
    const polygonRadius: number;
    const maxSubSteps: number;
    const maxTOIContacts: number;
    const maxLinearCorrection: number;
    const maxAngularCorrection: number;
    const maxTranslation: number;
    const maxTranslationSquared: number;
    const maxRotation: number;
    const maxRotationSquared: number;
    const baumgarte: number;
    const toiBaumgarte: number;
    const invalidParticleIndex: number;
    const maxParticleIndex: number;
    const particleStride: number;
    const minParticleWeight: number;
    const maxParticlePressure: number;
    const maxParticleForce: number;
    const maxTriadDistance: number;
    const maxTriadDistanceSquared: number;
    const minParticleSystemBufferCapacity: number;
    const barrierCollisionTime: number;
    const timeToSleep: number;
    const linearSleepTolerance: number;
    const angularSleepTolerance: number;
    class Version {
        major: number;
        minor: number;
        revision: number;
        constructor(major?: number, minor?: number, revision?: number);
        toString(): string;
    }
    const version: Version;
    const branch: string;
    const commit: string;
    function ParseInt(v: string): number;
    function ParseUInt(v: string): number;
    function MakeArray<T>(length: number, init: (i: number) => T): T[];
    function MakeNullArray<T>(length: number): Array<T | null>;
    function MakeNumberArray(length: number, init?: number): number[];
}
declare namespace b2 {
    const pi_over_180: number;
    const _180_over_pi: number;
    const two_pi: number;
    const Abs: (x: number) => number;
    function Min(a: number, b: number): number;
    function Max(a: number, b: number): number;
    function Clamp(a: number, lo: number, hi: number): number;
    function Swap<T>(a: T[], b: T[]): void;
    const IsValid: typeof isFinite;
    function Sq(n: number): number;
    function InvSqrt(n: number): number;
    const Sqrt: (x: number) => number;
    const Pow: (x: number, y: number) => number;
    function DegToRad(degrees: number): number;
    function RadToDeg(radians: number): number;
    const Cos: (x: number) => number;
    const Sin: (x: number) => number;
    const Acos: (x: number) => number;
    const Asin: (x: number) => number;
    const Atan2: (y: number, x: number) => number;
    function NextPowerOfTwo(x: number): number;
    function IsPowerOfTwo(x: number): boolean;
    function Random(): number;
    function RandomRange(lo: number, hi: number): number;
    interface XY {
        x: number;
        y: number;
    }
    class Vec2 implements XY {
        x: number;
        y: number;
        static readonly ZERO: Vec2;
        static readonly UNITX: Vec2;
        static readonly UNITY: Vec2;
        static readonly s_t0: Vec2;
        static readonly s_t1: Vec2;
        static readonly s_t2: Vec2;
        static readonly s_t3: Vec2;
        constructor(x?: number, y?: number);
        Clone(): Vec2;
        SetZero(): this;
        Set(x: number, y: number): this;
        Copy(other: XY): this;
        SelfAdd(v: XY): this;
        SelfAddXY(x: number, y: number): this;
        SelfSub(v: XY): this;
        SelfSubXY(x: number, y: number): this;
        SelfMul(s: number): this;
        SelfMulAdd(s: number, v: XY): this;
        SelfMulSub(s: number, v: XY): this;
        Dot(v: XY): number;
        Cross(v: XY): number;
        Length(): number;
        LengthSquared(): number;
        Normalize(): number;
        SelfNormalize(): this;
        SelfRotate(radians: number): this;
        SelfRotateCosSin(c: number, s: number): this;
        IsValid(): boolean;
        SelfCrossVS(s: number): this;
        SelfCrossSV(s: number): this;
        SelfMinV(v: XY): this;
        SelfMaxV(v: XY): this;
        SelfAbs(): this;
        SelfNeg(): this;
        SelfSkew(): this;
        static MakeArray(length: number): Vec2[];
        static AbsV<T extends XY>(v: XY, out: T): T;
        static MinV<T extends XY>(a: XY, b: XY, out: T): T;
        static MaxV<T extends XY>(a: XY, b: XY, out: T): T;
        static ClampV<T extends XY>(v: XY, lo: XY, hi: XY, out: T): T;
        static RotateV<T extends XY>(v: XY, radians: number, out: T): T;
        static DotVV(a: XY, b: XY): number;
        static CrossVV(a: XY, b: XY): number;
        static CrossVS<T extends XY>(v: XY, s: number, out: T): T;
        static CrossVOne<T extends XY>(v: XY, out: T): T;
        static CrossSV<T extends XY>(s: number, v: XY, out: T): T;
        static CrossOneV<T extends XY>(v: XY, out: T): T;
        static AddVV<T extends XY>(a: XY, b: XY, out: T): T;
        static SubVV<T extends XY>(a: XY, b: XY, out: T): T;
        static MulSV<T extends XY>(s: number, v: XY, out: T): T;
        static MulVS<T extends XY>(v: XY, s: number, out: T): T;
        static AddVMulSV<T extends XY>(a: XY, s: number, b: XY, out: T): T;
        static SubVMulSV<T extends XY>(a: XY, s: number, b: XY, out: T): T;
        static AddVCrossSV<T extends XY>(a: XY, s: number, v: XY, out: T): T;
        static MidVV<T extends XY>(a: XY, b: XY, out: T): T;
        static ExtVV<T extends XY>(a: XY, b: XY, out: T): T;
        static IsEqualToV(a: XY, b: XY): boolean;
        static DistanceVV(a: XY, b: XY): number;
        static DistanceSquaredVV(a: XY, b: XY): number;
        static NegV<T extends XY>(v: XY, out: T): T;
    }
    const Vec2_zero: Vec2;
    class TypedVec2 implements Vec2 {
        readonly data: Float32Array;
        x: number;
        y: number;
        constructor();
        constructor(data: Float32Array);
        constructor(x: number, y: number);
        Clone(): TypedVec2;
        SetZero(): this;
        Set(x: number, y: number): this;
        Copy(other: XY): this;
        SelfAdd(v: XY): this;
        SelfAddXY(x: number, y: number): this;
        SelfSub(v: XY): this;
        SelfSubXY(x: number, y: number): this;
        SelfMul(s: number): this;
        SelfMulAdd(s: number, v: XY): this;
        SelfMulSub(s: number, v: XY): this;
        Dot(v: XY): number;
        Cross(v: XY): number;
        Length(): number;
        LengthSquared(): number;
        Normalize(): number;
        SelfNormalize(): this;
        SelfRotate(radians: number): this;
        SelfRotateCosSin(c: number, s: number): this;
        IsValid(): boolean;
        SelfCrossVS(s: number): this;
        SelfCrossSV(s: number): this;
        SelfMinV(v: XY): this;
        SelfMaxV(v: XY): this;
        SelfAbs(): this;
        SelfNeg(): this;
        SelfSkew(): this;
    }
    interface XYZ extends XY {
        z: number;
    }
    class Vec3 implements XYZ {
        static readonly ZERO: Vec3;
        static readonly s_t0: Vec3;
        readonly data: Float32Array;
        x: number;
        y: number;
        z: number;
        constructor();
        constructor(data: Float32Array);
        constructor(x: number, y: number, z: number);
        Clone(): Vec3;
        SetZero(): this;
        SetXYZ(x: number, y: number, z: number): this;
        Copy(other: XYZ): this;
        SelfNeg(): this;
        SelfAdd(v: XYZ): this;
        SelfAddXYZ(x: number, y: number, z: number): this;
        SelfSub(v: XYZ): this;
        SelfSubXYZ(x: number, y: number, z: number): this;
        SelfMul(s: number): this;
        static DotV3V3(a: XYZ, b: XYZ): number;
        static CrossV3V3<T extends XYZ>(a: XYZ, b: XYZ, out: T): T;
    }
    class Mat22 {
        static readonly IDENTITY: Mat22;
        readonly ex: Vec2;
        readonly ey: Vec2;
        Clone(): Mat22;
        static FromVV(c1: XY, c2: XY): Mat22;
        static FromSSSS(r1c1: number, r1c2: number, r2c1: number, r2c2: number): Mat22;
        static FromAngle(radians: number): Mat22;
        SetSSSS(r1c1: number, r1c2: number, r2c1: number, r2c2: number): this;
        SetVV(c1: XY, c2: XY): this;
        SetAngle(radians: number): this;
        Copy(other: Mat22): this;
        SetIdentity(): this;
        SetZero(): this;
        GetAngle(): number;
        GetInverse(out: Mat22): Mat22;
        Solve<T extends XY>(b_x: number, b_y: number, out: T): T;
        SelfAbs(): this;
        SelfInv(): this;
        SelfAddM(M: Mat22): this;
        SelfSubM(M: Mat22): this;
        static AbsM(M: Mat22, out: Mat22): Mat22;
        static MulMV<T extends XY>(M: Mat22, v: XY, out: T): T;
        static MulTMV<T extends XY>(M: Mat22, v: XY, out: T): T;
        static AddMM(A: Mat22, B: Mat22, out: Mat22): Mat22;
        static MulMM(A: Mat22, B: Mat22, out: Mat22): Mat22;
        static MulTMM(A: Mat22, B: Mat22, out: Mat22): Mat22;
    }
    class Mat33 {
        static readonly IDENTITY: Mat33;
        readonly data: Float32Array;
        readonly ex: Vec3;
        readonly ey: Vec3;
        readonly ez: Vec3;
        Clone(): Mat33;
        SetVVV(c1: XYZ, c2: XYZ, c3: XYZ): this;
        Copy(other: Mat33): this;
        SetIdentity(): this;
        SetZero(): this;
        SelfAddM(M: Mat33): this;
        Solve33<T extends XYZ>(b_x: number, b_y: number, b_z: number, out: T): T;
        Solve22<T extends XY>(b_x: number, b_y: number, out: T): T;
        GetInverse22(M: Mat33): void;
        GetSymInverse33(M: Mat33): void;
        static MulM33V3<T extends XYZ>(A: Mat33, v: XYZ, out: T): T;
        static MulM33XYZ<T extends XYZ>(A: Mat33, x: number, y: number, z: number, out: T): T;
        static MulM33V2<T extends XY>(A: Mat33, v: XY, out: T): T;
        static MulM33XY<T extends XY>(A: Mat33, x: number, y: number, out: T): T;
    }
    class Rot {
        static readonly IDENTITY: Rot;
        s: number;
        c: number;
        constructor(angle?: number);
        Clone(): Rot;
        Copy(other: Rot): this;
        SetAngle(angle: number): this;
        SetIdentity(): this;
        GetAngle(): number;
        GetXAxis<T extends XY>(out: T): T;
        GetYAxis<T extends XY>(out: T): T;
        static MulRR(q: Rot, r: Rot, out: Rot): Rot;
        static MulTRR(q: Rot, r: Rot, out: Rot): Rot;
        static MulRV<T extends XY>(q: Rot, v: XY, out: T): T;
        static MulTRV<T extends XY>(q: Rot, v: XY, out: T): T;
    }
    class Transform {
        static readonly IDENTITY: Transform;
        readonly p: Vec2;
        readonly q: Rot;
        Clone(): Transform;
        Copy(other: Transform): this;
        SetIdentity(): this;
        SetPositionRotation(position: XY, q: Rot): this;
        SetPositionAngle(pos: XY, a: number): this;
        SetPosition(position: XY): this;
        SetPositionXY(x: number, y: number): this;
        SetRotation(rotation: Rot): this;
        SetRotationAngle(radians: number): this;
        GetPosition(): Vec2;
        GetRotation(): Rot;
        GetRotationAngle(): number;
        GetAngle(): number;
        static MulXV<T extends XY>(T: Transform, v: XY, out: T): T;
        static MulTXV<T extends XY>(T: Transform, v: XY, out: T): T;
        static MulXX(A: Transform, B: Transform, out: Transform): Transform;
        static MulTXX(A: Transform, B: Transform, out: Transform): Transform;
    }
    class Sweep {
        readonly localCenter: Vec2;
        readonly c0: Vec2;
        readonly c: Vec2;
        a0: number;
        a: number;
        alpha0: number;
        Clone(): Sweep;
        Copy(other: Sweep): this;
        GetTransform(xf: Transform, beta: number): Transform;
        Advance(alpha: number): void;
        Normalize(): void;
    }
}
declare namespace b2 {
    class DistanceProxy {
        readonly buffer: Vec2[];
        vertices: Vec2[];
        count: number;
        radius: number;
        Copy(other: DistanceProxy): this;
        Reset(): DistanceProxy;
        SetShape(shape: Shape, index: number): void;
        SetVerticesRadius(vertices: Vec2[], count: number, radius: number): void;
        GetSupport(d: Vec2): number;
        GetSupportVertex(d: Vec2): Vec2;
        GetVertexCount(): number;
        GetVertex(index: number): Vec2;
    }
    class SimplexCache {
        metric: number;
        count: number;
        readonly indexA: [number, number, number];
        readonly indexB: [number, number, number];
        Reset(): SimplexCache;
    }
    class DistanceInput {
        readonly proxyA: DistanceProxy;
        readonly proxyB: DistanceProxy;
        readonly transformA: Transform;
        readonly transformB: Transform;
        useRadii: boolean;
        Reset(): DistanceInput;
    }
    class DistanceOutput {
        readonly pointA: Vec2;
        readonly pointB: Vec2;
        distance: number;
        iterations: number;
        Reset(): DistanceOutput;
    }
    class ShapeCastInput {
        readonly proxyA: DistanceProxy;
        readonly proxyB: DistanceProxy;
        readonly transformA: Transform;
        readonly transformB: Transform;
        readonly translationB: Vec2;
    }
    class ShapeCastOutput {
        readonly point: Vec2;
        readonly normal: Vec2;
        lambda: number;
        iterations: number;
    }
    let gjkCalls: number;
    let gjkIters: number;
    let gjkMaxIters: number;
    function gjk_reset(): void;
    class SimplexVertex {
        readonly wA: Vec2;
        readonly wB: Vec2;
        readonly w: Vec2;
        a: number;
        indexA: number;
        indexB: number;
        Copy(other: SimplexVertex): SimplexVertex;
    }
    class Simplex {
        readonly v1: SimplexVertex;
        readonly v2: SimplexVertex;
        readonly v3: SimplexVertex;
        readonly vertices: SimplexVertex[];
        count: number;
        constructor();
        ReadCache(cache: SimplexCache, proxyA: DistanceProxy, transformA: Transform, proxyB: DistanceProxy, transformB: Transform): void;
        WriteCache(cache: SimplexCache): void;
        GetSearchDirection(out: Vec2): Vec2;
        GetClosestPoint(out: Vec2): Vec2;
        GetWitnessPoints(pA: Vec2, pB: Vec2): void;
        GetMetric(): number;
        Solve2(): void;
        Solve3(): void;
        private static s_e12;
        private static s_e13;
        private static s_e23;
    }
    function Distance(output: DistanceOutput, cache: SimplexCache, input: DistanceInput): void;
    function ShapeCast(output: ShapeCastOutput, input: ShapeCastInput): boolean;
}
declare namespace b2 {
    class Timer {
        start: number;
        Reset(): Timer;
        GetMilliseconds(): number;
    }
    class Counter {
        count: number;
        min_count: number;
        max_count: number;
        GetCount(): number;
        GetMinCount(): number;
        GetMaxCount(): number;
        ResetCount(): number;
        ResetMinCount(): void;
        ResetMaxCount(): void;
        Increment(): void;
        Decrement(): void;
    }
}
declare namespace b2 {
    enum ContactFeatureType {
        e_vertex = 0,
        e_face = 1
    }
    class ContactFeature {
        private _key;
        private _key_invalid;
        private _indexA;
        private _indexB;
        private _typeA;
        private _typeB;
        key: number;
        indexA: number;
        indexB: number;
        typeA: number;
        typeB: number;
    }
    class ContactID {
        readonly cf: ContactFeature;
        Copy(o: ContactID): ContactID;
        Clone(): ContactID;
        key: number;
    }
    class ManifoldPoint {
        readonly localPoint: Vec2;
        normalImpulse: number;
        tangentImpulse: number;
        readonly id: ContactID;
        static MakeArray(length: number): ManifoldPoint[];
        Reset(): void;
        Copy(o: ManifoldPoint): ManifoldPoint;
    }
    enum ManifoldType {
        e_unknown = -1,
        e_circles = 0,
        e_faceA = 1,
        e_faceB = 2
    }
    class Manifold {
        readonly points: ManifoldPoint[];
        readonly localNormal: Vec2;
        readonly localPoint: Vec2;
        type: ManifoldType;
        pointCount: number;
        Reset(): void;
        Copy(o: Manifold): Manifold;
        Clone(): Manifold;
    }
    class WorldManifold {
        readonly normal: Vec2;
        readonly points: Vec2[];
        readonly separations: number[];
        private static Initialize_s_pointA;
        private static Initialize_s_pointB;
        private static Initialize_s_cA;
        private static Initialize_s_cB;
        private static Initialize_s_planePoint;
        private static Initialize_s_clipPoint;
        Initialize(manifold: Manifold, xfA: Transform, radiusA: number, xfB: Transform, radiusB: number): void;
    }
    enum PointState {
        nullState = 0,
        addState = 1,
        persistState = 2,
        removeState = 3
    }
    function GetPointStates(state1: PointState[], state2: PointState[], manifold1: Manifold, manifold2: Manifold): void;
    class ClipVertex {
        readonly v: Vec2;
        readonly id: ContactID;
        static MakeArray(length: number): ClipVertex[];
        Copy(other: ClipVertex): ClipVertex;
    }
    class RayCastInput {
        readonly p1: Vec2;
        readonly p2: Vec2;
        maxFraction: number;
        Copy(o: RayCastInput): RayCastInput;
    }
    class RayCastOutput {
        readonly normal: Vec2;
        fraction: number;
        Copy(o: RayCastOutput): RayCastOutput;
    }
    class AABB {
        readonly lowerBound: Vec2;
        readonly upperBound: Vec2;
        private readonly cache_center;
        private readonly cache_extent;
        Copy(o: AABB): AABB;
        IsValid(): boolean;
        GetCenter(): Vec2;
        GetExtents(): Vec2;
        GetPerimeter(): number;
        Combine1(aabb: AABB): AABB;
        Combine2(aabb1: AABB, aab: AABB): AABB;
        static Combine(aabb1: AABB, aab: AABB, out: AABB): AABB;
        Contains(aabb: AABB): boolean;
        RayCast(output: RayCastOutput, input: RayCastInput): boolean;
        TestContain(point: XY): boolean;
        TestOverlap(other: AABB): boolean;
    }
    function TestOverlapAABB(a: AABB, b: AABB): boolean;
    function ClipSegmentToLine(vOut: [ClipVertex, ClipVertex], vIn: [ClipVertex, ClipVertex], normal: Vec2, offset: number, vertexIndexA: number): number;
    function TestOverlapShape(shapeA: Shape, indexA: number, shapeB: Shape, indexB: number, xfA: Transform, xfB: Transform): boolean;
}
declare namespace b2 {
    class MassData {
        mass: number;
        readonly center: Vec2;
        I: number;
    }
    enum ShapeType {
        e_unknown = -1,
        e_circleShape = 0,
        e_edgeShape = 1,
        e_polygonShape = 2,
        e_chainShape = 3,
        e_shapeTypeCount = 4
    }
    abstract class Shape {
        readonly type: ShapeType;
        radius: number;
        constructor(type: ShapeType, radius: number);
        abstract Clone(): Shape;
        Copy(other: Shape): Shape;
        GetType(): ShapeType;
        abstract GetChildCount(): number;
        abstract TestPoint(xf: Transform, p: XY): boolean;
        abstract ComputeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number;
        abstract RayCast(output: RayCastOutput, input: RayCastInput, transform: Transform, childIndex: number): boolean;
        abstract ComputeAABB(aabb: AABB, xf: Transform, childIndex: number): void;
        abstract ComputeMass(massData: MassData, density: number): void;
        abstract SetupDistanceProxy(proxy: DistanceProxy, index: number): void;
        abstract ComputeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number;
        abstract Dump(log: (format: string, ...args: any[]) => void): void;
    }
}
declare namespace b2 {
    let toiTime: number;
    let toiMaxTime: number;
    let toiCalls: number;
    let toiIters: number;
    let toiMaxIters: number;
    let toiRootIters: number;
    let toiMaxRootIters: number;
    function toi_reset(): void;
    class TOIInput {
        readonly proxyA: DistanceProxy;
        readonly proxyB: DistanceProxy;
        readonly sweepA: Sweep;
        readonly sweepB: Sweep;
        tMax: number;
    }
    enum TOIOutputState {
        e_unknown = 0,
        e_failed = 1,
        e_overlapped = 2,
        e_touching = 3,
        e_separated = 4
    }
    class TOIOutput {
        state: TOIOutputState;
        t: number;
    }
    enum SeparationFunctionType {
        e_unknown = -1,
        e_points = 0,
        e_faceA = 1,
        e_faceB = 2
    }
    class SeparationFunction {
        proxyA: DistanceProxy;
        proxyB: DistanceProxy;
        readonly sweepA: Sweep;
        readonly sweepB: Sweep;
        type: SeparationFunctionType;
        readonly localPoint: Vec2;
        readonly axis: Vec2;
        Initialize(cache: SimplexCache, proxyA: DistanceProxy, sweepA: Sweep, proxyB: DistanceProxy, sweepB: Sweep, t1: number): number;
        FindMinSeparation(indexA: [number], indexB: [number], t: number): number;
        Evaluate(indexA: number, indexB: number, t: number): number;
    }
    function TimeOfImpact(output: TOIOutput, input: TOIInput): void;
}
declare namespace b2 {
    interface RGB {
        r: number;
        g: number;
        b: number;
    }
    interface RGBA extends RGB {
        a: number;
    }
    class Color implements RGBA {
        r: number;
        g: number;
        b: number;
        a: number;
        static readonly ZERO: Color;
        static readonly RED: Color;
        static readonly GREEN: Color;
        static readonly BLUE: Color;
        constructor(r?: number, g?: number, b?: number, a?: number);
        Clone(): Color;
        Copy(other: RGBA): this;
        IsEqual(color: RGBA): boolean;
        IsZero(): boolean;
        Set(r: number, g: number, b: number, a?: number): void;
        SetByteRGB(r: number, g: number, b: number): this;
        SetByteRGBA(r: number, g: number, b: number, a: number): this;
        SetRGB(rr: number, gg: number, bb: number): this;
        SetRGBA(rr: number, gg: number, bb: number, aa: number): this;
        SelfAdd(color: RGBA): this;
        Add<T extends RGBA>(color: RGBA, out: T): T;
        SelfSub(color: RGBA): this;
        Sub<T extends RGBA>(color: RGBA, out: T): T;
        SelfMul(s: number): this;
        Mul<T extends RGBA>(s: number, out: T): T;
        Mix(mixColor: RGBA, strength: number): void;
        static MixColors(colorA: RGBA, colorB: RGBA, strength: number): void;
        MakeStyleString(alpha?: number): string;
        static MakeStyleString(r: number, g: number, b: number, a?: number): string;
    }
    class TypedColor implements Color {
        readonly data: Float32Array;
        r: number;
        g: number;
        b: number;
        a: number;
        constructor();
        constructor(data: Float32Array);
        constructor(rr: number, gg: number, bb: number);
        constructor(rr: number, gg: number, bb: number, aa: number);
        Clone(): TypedColor;
        Copy(other: RGBA): this;
        IsEqual(color: RGBA): boolean;
        IsZero(): boolean;
        Set(r: number, g: number, b: number, a?: number): void;
        SetByteRGB(r: number, g: number, b: number): this;
        SetByteRGBA(r: number, g: number, b: number, a: number): this;
        SetRGB(rr: number, gg: number, bb: number): this;
        SetRGBA(rr: number, gg: number, bb: number, aa: number): this;
        SelfAdd(color: RGBA): this;
        Add<T extends RGBA>(color: RGBA, out: T): T;
        SelfSub(color: RGBA): this;
        Sub<T extends RGBA>(color: RGBA, out: T): T;
        SelfMul(s: number): this;
        Mul<T extends RGBA>(s: number, out: T): T;
        Mix(mixColor: RGBA, strength: number): void;
        MakeStyleString(alpha?: number): string;
    }
    enum DrawFlags {
        e_none = 0,
        e_shapeBit = 1,
        e_jointBit = 2,
        e_aabbBit = 4,
        e_pairBit = 8,
        e_centerOfMassBit = 16,
        e_particleBit = 32,
        e_controllerBit = 64,
        e_all = 63
    }
    class Draw {
        drawFlags: DrawFlags;
        SetFlags(flags: DrawFlags): void;
        GetFlags(): DrawFlags;
        AppendFlags(flags: DrawFlags): void;
        ClearFlags(flags: DrawFlags): void;
        PushTransform(xf: Transform): void;
        PopTransform(xf: Transform): void;
        DrawPolygon(vertices: XY[], vertexCount: number, color: RGBA): void;
        DrawSolidPolygon(vertices: XY[], vertexCount: number, color: RGBA): void;
        DrawCircle(center: XY, radius: number, color: RGBA): void;
        DrawSolidCircle(center: XY, radius: number, axis: XY, color: RGBA): void;
        DrawParticles(centers: XY[], radius: number, colors: RGBA[] | null, count: number): void;
        DrawSegment(p1: XY, p2: XY, color: RGBA): void;
        DrawTransform(xf: Transform): void;
        DrawPoint(p: XY, size: number, color: RGBA): void;
    }
}
declare namespace b2 {
    class Profile {
        step: number;
        collide: number;
        solve: number;
        solveInit: number;
        solveVelocity: number;
        solvePosition: number;
        broadphase: number;
        solveTOI: number;
        Reset(): this;
    }
    class TimeStep {
        dt: number;
        inv_dt: number;
        dtRatio: number;
        velocityIterations: number;
        positionIterations: number;
        particleIterations: number;
        warmStarting: boolean;
        Copy(step: TimeStep): TimeStep;
    }
    class Position {
        readonly c: Vec2;
        a: number;
        static MakeArray(length: number): Position[];
    }
    class Velocity {
        readonly v: Vec2;
        w: number;
        static MakeArray(length: number): Velocity[];
    }
    class SolverData {
        readonly step: TimeStep;
        positions: Position[];
        velocities: Velocity[];
    }
}
declare namespace b2 {
    class EdgeShape extends Shape {
        readonly vertex1: Vec2;
        readonly vertex2: Vec2;
        readonly vertex0: Vec2;
        readonly vertex3: Vec2;
        oneSided: boolean;
        constructor();
        SetOneSided(v0: XY, v1: XY, v2: XY, v3: XY): EdgeShape;
        SetTwoSided(v1: XY, v2: XY): EdgeShape;
        Clone(): EdgeShape;
        Copy(other: EdgeShape): EdgeShape;
        GetChildCount(): number;
        TestPoint(xf: Transform, p: XY): boolean;
        private static ComputeDistance_s_v1;
        private static ComputeDistance_s_v2;
        private static ComputeDistance_s_d;
        private static ComputeDistance_s_s;
        ComputeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number;
        private static RayCast_s_p1;
        private static RayCast_s_p2;
        private static RayCast_s_d;
        private static RayCast_s_e;
        private static RayCast_s_q;
        private static RayCast_s_r;
        RayCast(output: RayCastOutput, input: RayCastInput, xf: Transform, childIndex: number): boolean;
        private static ComputeAABB_s_v1;
        private static ComputeAABB_s_v2;
        ComputeAABB(aabb: AABB, xf: Transform, childIndex: number): void;
        ComputeMass(massData: MassData, density: number): void;
        SetupDistanceProxy(proxy: DistanceProxy, index: number): void;
        ComputeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number;
        Dump(log: (format: string, ...args: any[]) => void): void;
    }
}
declare namespace b2 {
    /**
     * A controller edge is used to connect bodies and controllers
     * together in a bipartite graph.
     */
    class ControllerEdge {
        readonly controller: Controller;
        readonly body: Body;
        prevBody: ControllerEdge | null;
        nextBody: ControllerEdge | null;
        prevController: ControllerEdge | null;
        nextController: ControllerEdge | null;
        constructor(controller: Controller, body: Body);
    }
    /**
     * Base class for controllers. Controllers are a convience for
     * encapsulating common per-step functionality.
     */
    abstract class Controller {
        bodyList: ControllerEdge | null;
        bodyCount: number;
        prev: Controller | null;
        next: Controller | null;
        /**
         * Controllers override this to implement per-step functionality.
         */
        abstract Step(step: TimeStep): void;
        /**
         * Controllers override this to provide debug drawing.
         */
        abstract Draw(debugDraw: Draw): void;
        /**
         * Get the next controller in the world's body list.
         */
        GetNext(): Controller | null;
        /**
         * Get the previous controller in the world's body list.
         */
        GetPrev(): Controller | null;
        /**
         * Get the parent world of this body.
         */
        /**
         * Get the attached body list
         */
        GetBodyList(): ControllerEdge | null;
        /**
         * Adds a body to the controller list.
         */
        AddBody(body: Body): void;
        /**
         * Removes a body from the controller list.
         */
        RemoveBody(body: Body): void;
        /**
         * Removes all bodies from the controller list.
         */
        Clear(): void;
    }
}
declare namespace b2 {
    function MixFriction(friction1: number, friction2: number): number;
    function MixRestitution(restitution1: number, restitution2: number): number;
    function MixRestitutionThreshold(threshold1: number, threshold2: number): number;
    class ContactEdge {
        private _other;
        other: Body;
        readonly contact: Contact;
        prev: ContactEdge | null;
        next: ContactEdge | null;
        constructor(contact: Contact);
        Reset(): void;
    }
    abstract class Contact<A extends Shape = Shape, B extends Shape = Shape> {
        islandFlag: boolean;
        touchingFlag: boolean;
        enabledFlag: boolean;
        filterFlag: boolean;
        bulletHitFlag: boolean;
        toiFlag: boolean;
        prev: Contact | null;
        next: Contact | null;
        readonly nodeA: ContactEdge;
        readonly nodeB: ContactEdge;
        fixtureA: Fixture;
        fixtureB: Fixture;
        indexA: number;
        indexB: number;
        manifold: Manifold;
        toiCount: number;
        toi: number;
        friction: number;
        restitution: number;
        restitutionThreshold: number;
        tangentSpeed: number;
        oldManifold: Manifold;
        GetManifold(): Manifold;
        GetWorldManifold(worldManifold: WorldManifold): void;
        IsTouching(): boolean;
        SetEnabled(flag: boolean): void;
        IsEnabled(): boolean;
        GetNext(): Contact | null;
        GetFixtureA(): Fixture;
        GetChildIndexA(): number;
        GetShapeA(): A;
        GetFixtureB(): Fixture;
        GetChildIndexB(): number;
        GetShapeB(): B;
        abstract Evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
        FlagForFiltering(): void;
        SetFriction(friction: number): void;
        GetFriction(): number;
        ResetFriction(): void;
        SetRestitution(restitution: number): void;
        GetRestitution(): number;
        ResetRestitution(): void;
        SetRestitutionThreshold(threshold: number): void;
        GetRestitutionThreshold(): number;
        ResetRestitutionThreshold(): void;
        SetTangentSpeed(speed: number): void;
        GetTangentSpeed(): number;
        Reset(fixtureA: Fixture, indexA: number, fixtureB: Fixture, indexB: number): void;
        Update(listener: ContactListener): void;
        private static ComputeTOI_s_input;
        private static ComputeTOI_s_output;
        ComputeTOI(sweepA: Sweep, sweepB: Sweep): number;
    }
}
declare namespace b2 {
    let g_blockSolve: boolean;
    function get_g_blockSolve(): boolean;
    function set_g_blockSolve(value: boolean): void;
    class VelocityConstraintPoint {
        readonly rA: Vec2;
        readonly rB: Vec2;
        normalImpulse: number;
        tangentImpulse: number;
        normalMass: number;
        tangentMass: number;
        velocityBias: number;
        static MakeArray(length: number): VelocityConstraintPoint[];
    }
    class ContactVelocityConstraint {
        readonly points: VelocityConstraintPoint[];
        readonly normal: Vec2;
        readonly tangent: Vec2;
        readonly normalMass: Mat22;
        readonly K: Mat22;
        indexA: number;
        indexB: number;
        invMassA: number;
        invMassB: number;
        invIA: number;
        invIB: number;
        friction: number;
        restitution: number;
        threshold: number;
        tangentSpeed: number;
        pointCount: number;
        contactIndex: number;
        static MakeArray(length: number): ContactVelocityConstraint[];
    }
    class ContactPositionConstraint {
        readonly localPoints: Vec2[];
        readonly localNormal: Vec2;
        readonly localPoint: Vec2;
        indexA: number;
        indexB: number;
        invMassA: number;
        invMassB: number;
        readonly localCenterA: Vec2;
        readonly localCenterB: Vec2;
        invIA: number;
        invIB: number;
        type: ManifoldType;
        radiusA: number;
        radiusB: number;
        pointCount: number;
        static MakeArray(length: number): ContactPositionConstraint[];
    }
    class ContactSolverDef {
        readonly step: TimeStep;
        contacts: Contact[];
        count: number;
        positions: Position[];
        velocities: Velocity[];
    }
    class PositionSolverManifold {
        readonly normal: Vec2;
        readonly point: Vec2;
        separation: number;
        private static Initialize_s_pointA;
        private static Initialize_s_pointB;
        private static Initialize_s_planePoint;
        private static Initialize_s_clipPoint;
        Initialize(pc: ContactPositionConstraint, xfA: Transform, xfB: Transform, index: number): void;
    }
    class ContactSolver {
        readonly step: TimeStep;
        positions: Position[];
        velocities: Velocity[];
        readonly positionConstraints: ContactPositionConstraint[];
        readonly velocityConstraints: ContactVelocityConstraint[];
        contacts: Contact[];
        count: number;
        Initialize(def: ContactSolverDef): ContactSolver;
        private static InitializeVelocityConstraints_s_xfA;
        private static InitializeVelocityConstraints_s_xfB;
        private static InitializeVelocityConstraints_s_worldManifold;
        InitializeVelocityConstraints(): void;
        private static WarmStart_s_P;
        WarmStart(): void;
        private static SolveVelocityConstraints_s_dv;
        private static SolveVelocityConstraints_s_dv1;
        private static SolveVelocityConstraints_s_dv2;
        private static SolveVelocityConstraints_s_P;
        private static SolveVelocityConstraints_s_a;
        private static SolveVelocityConstraints_s_b;
        private static SolveVelocityConstraints_s_x;
        private static SolveVelocityConstraints_s_d;
        private static SolveVelocityConstraints_s_P1;
        private static SolveVelocityConstraints_s_P2;
        private static SolveVelocityConstraints_s_P1P2;
        SolveVelocityConstraints(): void;
        StoreImpulses(): void;
        private static SolvePositionConstraints_s_xfA;
        private static SolvePositionConstraints_s_xfB;
        private static SolvePositionConstraints_s_psm;
        private static SolvePositionConstraints_s_rA;
        private static SolvePositionConstraints_s_rB;
        private static SolvePositionConstraints_s_P;
        SolvePositionConstraints(): boolean;
        private static SolveTOIPositionConstraints_s_xfA;
        private static SolveTOIPositionConstraints_s_xfB;
        private static SolveTOIPositionConstraints_s_psm;
        private static SolveTOIPositionConstraints_s_rA;
        private static SolveTOIPositionConstraints_s_rB;
        private static SolveTOIPositionConstraints_s_P;
        SolveTOIPositionConstraints(toiIndexA: number, toiIndexB: number): boolean;
    }
}
declare namespace b2 {
    interface IFilter {
        categoryBits: number;
        maskBits: number;
        groupIndex?: number;
    }
    class Filter implements IFilter {
        static readonly DEFAULT: Filter;
        categoryBits: number;
        maskBits: number;
        groupIndex: number;
        Clone(): Filter;
        Copy(other: IFilter): this;
    }
    interface IFixtureDef {
        shape: Shape;
        userData?: any;
        friction?: number;
        restitution?: number;
        restitutionThreshold?: number;
        density?: number;
        isSensor?: boolean;
        filter?: IFilter;
    }
    class FixtureDef implements IFixtureDef {
        shape: Shape;
        userData: any;
        friction: number;
        restitution: number;
        restitutionThreshold: number;
        density: number;
        isSensor: boolean;
        readonly filter: Filter;
    }
    class FixtureProxy {
        readonly aabb: AABB;
        readonly fixture: Fixture;
        readonly childIndex: number;
        treeNode: TreeNode<FixtureProxy>;
        constructor(fixture: Fixture, childIndex: number);
        Reset(): void;
        Touch(): void;
        private static Synchronize_s_aabb1;
        private static Synchronize_s_aab;
        private static Synchronize_s_displacement;
        Synchronize(transform1: Transform, transform2: Transform): void;
    }
    class Fixture {
        density: number;
        next: Fixture | null;
        readonly body: Body;
        readonly shape: Shape;
        friction: number;
        restitution: number;
        restitutionThreshold: number;
        readonly proxies: FixtureProxy[];
        readonly proxyCount: number;
        readonly filter: Filter;
        isSensor: boolean;
        userData: any;
        constructor(body: Body, def: IFixtureDef);
        Reset(): void;
        GetType(): ShapeType;
        GetShape(): Shape;
        SetSensor(sensor: boolean): void;
        IsSensor(): boolean;
        SetFilterData(filter: Filter): void;
        GetFilterData(): Filter;
        Refilter(): void;
        GetBody(): Body;
        GetNext(): Fixture | null;
        GetUserData(): any;
        SetUserData(data: any): void;
        TestPoint(p: XY): boolean;
        ComputeDistance(p: Vec2, normal: Vec2, childIndex: number): number;
        RayCast(output: RayCastOutput, input: RayCastInput, childIndex: number): boolean;
        GetMassData(massData?: MassData): MassData;
        SetDensity(density: number): void;
        GetDensity(): number;
        GetFriction(): number;
        SetFriction(friction: number): void;
        GetRestitution(): number;
        SetRestitution(restitution: number): void;
        GetRestitutionThreshold(): number;
        SetRestitutionThreshold(threshold: number): void;
        GetAABB(childIndex: number): AABB;
        Dump(log: (format: string, ...args: any[]) => void, bodyIndex: number): void;
        CreateProxies(): void;
        DestroyProxies(): void;
        TouchProxies(): void;
        SynchronizeProxies(transform1: Transform, transform2: Transform): void;
    }
}
declare namespace b2 {
    enum JointType {
        e_unknownJoint = 0,
        e_revoluteJoint = 1,
        e_prismaticJoint = 2,
        e_distanceJoint = 3,
        e_pulleyJoint = 4,
        e_mouseJoint = 5,
        e_gearJoint = 6,
        e_wheelJoint = 7,
        e_weldJoint = 8,
        e_frictionJoint = 9,
        e_ropeJoint = 10,
        e_motorJoint = 11,
        e_areaJoint = 12
    }
    class Jacobian {
        readonly linear: Vec2;
        angularA: number;
        angularB: number;
        SetZero(): Jacobian;
        Set(x: XY, a1: number, a2: number): Jacobian;
    }
    class JointEdge {
        private _other;
        other: Body;
        readonly joint: Joint;
        prev: JointEdge | null;
        next: JointEdge | null;
        constructor(joint: Joint);
        Reset(): void;
    }
    interface IJointDef {
        type: JointType;
        userData?: any;
        bodyA: Body;
        bodyB: Body;
        collideConnected?: boolean;
    }
    abstract class JointDef implements IJointDef {
        readonly type: JointType;
        userData: any;
        bodyA: Body;
        bodyB: Body;
        collideConnected: boolean;
        constructor(type: JointType);
    }
    function LinearStiffness(def: {
        stiffness: number;
        damping: number;
    }, frequencyHertz: number, dampingRatio: number, bodyA: Body, bodyB: Body): void;
    function AngularStiffness(def: {
        stiffness: number;
        damping: number;
    }, frequencyHertz: number, dampingRatio: number, bodyA: Body, bodyB: Body): void;
    abstract class Joint {
        readonly type: JointType;
        prev: Joint | null;
        next: Joint | null;
        readonly edgeA: JointEdge;
        readonly edgeB: JointEdge;
        bodyA: Body;
        bodyB: Body;
        index: number;
        islandFlag: boolean;
        collideConnected: boolean;
        userData: any;
        constructor(def: IJointDef);
        GetType(): JointType;
        GetBodyA(): Body;
        GetBodyB(): Body;
        abstract GetAnchorA<T extends XY>(out: T): T;
        abstract GetAnchorB<T extends XY>(out: T): T;
        abstract GetReactionForce<T extends XY>(inv_dt: number, out: T): T;
        abstract GetReactionTorque(inv_dt: number): number;
        GetNext(): Joint | null;
        GetUserData(): any;
        SetUserData(data: any): void;
        IsEnabled(): boolean;
        GetCollideConnected(): boolean;
        Dump(log: (format: string, ...args: any[]) => void): void;
        ShiftOrigin(newOrigin: XY): void;
        private static Draw_s_p1;
        private static Draw_s_p2;
        private static Draw_s_color;
        private static Draw_s_c;
        Draw(draw: Draw): void;
        abstract InitVelocityConstraints(data: SolverData): void;
        abstract SolveVelocityConstraints(data: SolverData): void;
        abstract SolvePositionConstraints(data: SolverData): boolean;
    }
}
declare namespace b2 {
    class DestructionListener {
        SayGoodbyeJoint(joint: Joint): void;
        SayGoodbyeFixture(fixture: Fixture): void;
        SayGoodbyeParticleGroup(group: ParticleGroup): void;
        SayGoodbyeParticle(system: ParticleSystem, index: number): void;
    }
    class ContactFilter {
        ShouldCollide(fixtureA: Fixture, fixtureB: Fixture): boolean;
        ShouldCollideFixtureParticle(fixture: Fixture, system: ParticleSystem, index: number): boolean;
        ShouldCollideParticleParticle(system: ParticleSystem, indexA: number, indexB: number): boolean;
        static readonly defaultFilter: ContactFilter;
    }
    class ContactImpulse {
        normalImpulses: number[];
        tangentImpulses: number[];
        count: number;
    }
    class ContactListener {
        BeginContact(contact: Contact): void;
        EndContact(contact: Contact): void;
        BeginContactFixtureParticle(system: ParticleSystem, contact: ParticleBodyContact): void;
        EndContactFixtureParticle(system: ParticleSystem, contact: ParticleBodyContact): void;
        BeginContactParticleParticle(system: ParticleSystem, contact: ParticleContact): void;
        EndContactParticleParticle(system: ParticleSystem, contact: ParticleContact): void;
        PreSolve(contact: Contact, oldManifold: Manifold): void;
        PostSolve(contact: Contact, impulse: ContactImpulse): void;
        static readonly defaultListener: ContactListener;
    }
    class QueryCallback {
        ReportFixture(fixture: Fixture): boolean;
        ReportParticle(system: ParticleSystem, index: number): boolean;
        ShouldQueryParticleSystem(system: ParticleSystem): boolean;
    }
    type QueryCallbackFunction = (fixture: Fixture) => boolean;
    class RayCastCallback {
        ReportFixture(fixture: Fixture, point: Vec2, normal: Vec2, fraction: number): number;
        ReportParticle(system: ParticleSystem, index: number, point: Vec2, normal: Vec2, fraction: number): number;
        ShouldQueryParticleSystem(system: ParticleSystem): boolean;
    }
    type RayCastCallbackFunction = (fixture: Fixture, point: Vec2, normal: Vec2, fraction: number) => number;
}
declare namespace b2 {
    /**
     * The particle type. Can be combined with the | operator.
     */
    enum ParticleFlag {
        waterParticle = 0,
        zombieParticle = 2,
        wallParticle = 4,
        springParticle = 8,
        elasticParticle = 16,
        viscousParticle = 32,
        powderParticle = 64,
        tensileParticle = 128,
        colorMixingParticle = 256,
        destructionListenerParticle = 512,
        barrierParticle = 1024,
        staticPressureParticle = 2048,
        reactiveParticle = 4096,
        repulsiveParticle = 8192,
        fixtureContactListenerParticle = 16384,
        particleContactListenerParticle = 32768,
        fixtureContactFilterParticle = 65536,
        particleContactFilterParticle = 131072
    }
    interface IParticleDef {
        flags?: ParticleFlag;
        position?: XY;
        velocity?: XY;
        color?: RGBA;
        lifetime?: number;
        userData?: any;
        group?: ParticleGroup | null;
    }
    class ParticleDef implements IParticleDef {
        flags: ParticleFlag;
        readonly position: Vec2;
        readonly velocity: Vec2;
        readonly color: Color;
        lifetime: number;
        userData: any;
        group: ParticleGroup | null;
    }
    function CalculateParticleIterations(gravity: number, radius: number, timeStep: number): number;
    class ParticleHandle {
        index: number;
        GetIndex(): number;
        SetIndex(index: number): void;
    }
}
declare namespace b2 {
    class Pair<T> {
        proxyA: TreeNode<T>;
        proxyB: TreeNode<T>;
        constructor(proxyA: TreeNode<T>, proxyB: TreeNode<T>);
    }
    class BroadPhase<T> {
        readonly tree: DynamicTree<T>;
        proxyCount: number;
        moveCount: number;
        readonly moveBuffer: Array<TreeNode<T> | null>;
        pairCount: number;
        readonly pairBuffer: Array<Pair<T>>;
        CreateProxy(aabb: AABB, userData: T): TreeNode<T>;
        DestroyProxy(proxy: TreeNode<T>): void;
        MoveProxy(proxy: TreeNode<T>, aabb: AABB, displacement: Vec2): void;
        TouchProxy(proxy: TreeNode<T>): void;
        GetProxyCount(): number;
        UpdatePairs(callback: (a: T, b: T) => void): void;
        Query(aabb: AABB, callback: (node: TreeNode<T>) => boolean): void;
        QueryPoint(point: XY, callback: (node: TreeNode<T>) => boolean): void;
        RayCast(input: RayCastInput, callback: (input: RayCastInput, node: TreeNode<T>) => number): void;
        GetTreeHeight(): number;
        GetTreeBalance(): number;
        GetTreeQuality(): number;
        ShiftOrigin(newOrigin: XY): void;
        BufferMove(proxy: TreeNode<T>): void;
        UnBufferMove(proxy: TreeNode<T>): void;
    }
}
declare namespace b2 {
    class ChainShape extends Shape {
        vertices: Vec2[];
        count: number;
        readonly prevVertex: Vec2;
        readonly nextVertex: Vec2;
        constructor();
        CreateLoop(vertices: XY[]): ChainShape;
        CreateLoop(vertices: XY[], count: number): ChainShape;
        CreateLoop(vertices: number[]): ChainShape;
        private _CreateLoop;
        CreateChain(vertices: XY[], prevVertex: XY, nextVertex: XY): ChainShape;
        CreateChain(vertices: XY[], count: number, prevVertex: XY, nextVertex: XY): ChainShape;
        CreateChain(vertices: number[], prevVertex: XY, nextVertex: XY): ChainShape;
        private _CreateChain;
        Clone(): ChainShape;
        Copy(other: ChainShape): ChainShape;
        GetChildCount(): number;
        GetChildEdge(edge: EdgeShape, index: number): void;
        TestPoint(xf: Transform, p: XY): boolean;
        private static ComputeDistance_s_edgeShape;
        ComputeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number;
        private static RayCast_s_edgeShape;
        RayCast(output: RayCastOutput, input: RayCastInput, xf: Transform, childIndex: number): boolean;
        private static ComputeAABB_s_v1;
        private static ComputeAABB_s_v2;
        private static ComputeAABB_s_lower;
        private static ComputeAABB_s_upper;
        ComputeAABB(aabb: AABB, xf: Transform, childIndex: number): void;
        ComputeMass(massData: MassData, density: number): void;
        SetupDistanceProxy(proxy: DistanceProxy, index: number): void;
        ComputeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number;
        Dump(log: (format: string, ...args: any[]) => void): void;
    }
}
declare namespace b2 {
    class CircleShape extends Shape {
        readonly p: Vec2;
        constructor(radius?: number);
        Set(position: XY, radius?: number): this;
        Clone(): CircleShape;
        Copy(other: CircleShape): CircleShape;
        GetChildCount(): number;
        private static TestPoint_s_center;
        private static TestPoint_s_d;
        TestPoint(transform: Transform, p: XY): boolean;
        private static ComputeDistance_s_center;
        ComputeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number;
        private static RayCast_s_position;
        private static RayCast_s_s;
        private static RayCast_s_r;
        RayCast(output: RayCastOutput, input: RayCastInput, transform: Transform, childIndex: number): boolean;
        private static ComputeAABB_s_p;
        ComputeAABB(aabb: AABB, transform: Transform, childIndex: number): void;
        ComputeMass(massData: MassData, density: number): void;
        SetupDistanceProxy(proxy: DistanceProxy, index: number): void;
        ComputeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number;
        Dump(log: (format: string, ...args: any[]) => void): void;
    }
}
declare namespace b2 {
    function CollideCircles(manifold: Manifold, circleA: CircleShape, xfA: Transform, circleB: CircleShape, xfB: Transform): void;
    function CollidePolygonAndCircle(manifold: Manifold, polygonA: PolygonShape, xfA: Transform, circleB: CircleShape, xfB: Transform): void;
}
declare namespace b2 {
    function CollideEdgeAndCircle(manifold: Manifold, edgeA: EdgeShape, xfA: Transform, circleB: CircleShape, xfB: Transform): void;
    function CollideEdgeAndPolygon(manifold: Manifold, edgeA: EdgeShape, xfA: Transform, polygonB: PolygonShape, xfB: Transform): void;
}
declare namespace b2 {
    function CollidePolygons(manifold: Manifold, polyA: PolygonShape, xfA: Transform, polyB: PolygonShape, xfB: Transform): void;
}
declare namespace b2 {
    class TreeNode<T> {
        readonly id: number;
        readonly aabb: AABB;
        private _userData;
        userData: T;
        parent: TreeNode<T> | null;
        child1: TreeNode<T> | null;
        child2: TreeNode<T> | null;
        height: number;
        moved: boolean;
        constructor(id?: number);
        Reset(): void;
        IsLeaf(): boolean;
    }
    class DynamicTree<T> {
        root: TreeNode<T> | null;
        freeList: TreeNode<T> | null;
        insertionCount: number;
        readonly stack: GrowableStack<TreeNode<T>>;
        static readonly s_r: Vec2;
        static readonly s_v: Vec2;
        static readonly s_abs_v: Vec2;
        static readonly s_segmentAABB: AABB;
        static readonly s_subInput: RayCastInput;
        static readonly s_combinedAABB: AABB;
        static readonly s_aabb: AABB;
        Query(aabb: AABB, callback: (node: TreeNode<T>) => boolean): void;
        QueryPoint(point: XY, callback: (node: TreeNode<T>) => boolean): void;
        RayCast(input: RayCastInput, callback: (input: RayCastInput, node: TreeNode<T>) => number): void;
        static s_node_id: number;
        AllocateNode(): TreeNode<T>;
        FreeNode(node: TreeNode<T>): void;
        CreateProxy(aabb: AABB, userData: T): TreeNode<T>;
        DestroyProxy(node: TreeNode<T>): void;
        private static MoveProxy_s_fatAABB;
        private static MoveProxy_s_hugeAABB;
        MoveProxy(node: TreeNode<T>, aabb: AABB, displacement: Vec2): boolean;
        InsertLeaf(leaf: TreeNode<T>): void;
        RemoveLeaf(leaf: TreeNode<T>): void;
        Balance(A: TreeNode<T>): TreeNode<T>;
        GetHeight(): number;
        private static GetAreaNode;
        GetAreaRatio(): number;
        static ComputeHeightNode<T>(node: TreeNode<T> | null): number;
        ComputeHeight(): number;
        ValidateStructure(node: TreeNode<T> | null): void;
        ValidateMetrics(node: TreeNode<T> | null): void;
        Validate(): void;
        private static GetMaxBalanceNode;
        GetMaxBalance(): number;
        RebuildBottomUp(): void;
        private static ShiftOriginNode;
        ShiftOrigin(newOrigin: XY): void;
    }
}
declare namespace b2 {
    class PolygonShape extends Shape {
        readonly centroid: Vec2;
        vertices: Vec2[];
        normals: Vec2[];
        count: number;
        constructor();
        Clone(): PolygonShape;
        Copy(other: PolygonShape): PolygonShape;
        GetChildCount(): number;
        private static Set_s_r;
        private static Set_s_v;
        Set(vertices: XY[]): PolygonShape;
        Set(vertices: XY[], count: number): PolygonShape;
        Set(vertices: number[]): PolygonShape;
        _Set(vertices: (index: number) => XY, count: number): PolygonShape;
        SetAsBox(hx: number, hy: number, center?: XY, angle?: number): PolygonShape;
        private static TestPoint_s_pLocal;
        TestPoint(xf: Transform, p: XY): boolean;
        private static ComputeDistance_s_pLocal;
        private static ComputeDistance_s_normalForMaxDistance;
        private static ComputeDistance_s_minDistance;
        private static ComputeDistance_s_distance;
        ComputeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number;
        private static RayCast_s_p1;
        private static RayCast_s_p2;
        private static RayCast_s_d;
        RayCast(output: RayCastOutput, input: RayCastInput, xf: Transform, childIndex: number): boolean;
        private static ComputeAABB_s_v;
        ComputeAABB(aabb: AABB, xf: Transform, childIndex: number): void;
        private static ComputeMass_s_center;
        private static ComputeMass_s_s;
        private static ComputeMass_s_e1;
        private static ComputeMass_s_e2;
        ComputeMass(massData: MassData, density: number): void;
        private static Validate_s_e;
        private static Validate_s_v;
        Validate(): boolean;
        SetupDistanceProxy(proxy: DistanceProxy, index: number): void;
        private static ComputeSubmergedArea_s_normalL;
        private static ComputeSubmergedArea_s_md;
        private static ComputeSubmergedArea_s_intoVec;
        private static ComputeSubmergedArea_s_outoVec;
        private static ComputeSubmergedArea_s_center;
        ComputeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number;
        Dump(log: (format: string, ...args: any[]) => void): void;
        private static ComputeCentroid_s_s;
        private static ComputeCentroid_s_p1;
        private static ComputeCentroid_s_p2;
        private static ComputeCentroid_s_p3;
        private static ComputeCentroid_s_e1;
        private static ComputeCentroid_s_e2;
        static ComputeCentroid(vs: Vec2[], count: number, out: Vec2): Vec2;
    }
}
declare namespace b2 {
    class BlockAllocator {
    }
}
declare namespace b2 {
    class GrowableStack<T> {
        stack: Array<T | null>;
        count: number;
        constructor(N: number);
        Reset(): this;
        Push(element: T): void;
        Pop(): T | null;
        GetCount(): number;
    }
}
declare namespace b2 {
    function Alloc(size: number): any;
    function Free(mem: any): void;
    function Log(message: string, ...args: any[]): void;
}
declare namespace b2 {
    class StackAllocator {
    }
}
declare namespace b2 {
    /**
     * Calculates buoyancy forces for fluids in the form of a half
     * plane.
     */
    class BuoyancyController extends Controller {
        /**
         * The outer surface normal
         */
        readonly normal: Vec2;
        /**
         * The height of the fluid surface along the normal
         */
        offset: number;
        /**
         * The fluid density
         */
        density: number;
        /**
         * Fluid velocity, for drag calculations
         */
        readonly velocity: Vec2;
        /**
         * Linear drag co-efficient
         */
        linearDrag: number;
        /**
         * Angular drag co-efficient
         */
        angularDrag: number;
        /**
         * If false, bodies are assumed to be uniformly dense, otherwise
         * use the shapes densities
         */
        useDensity: boolean;
        /**
         * If true, gravity is taken from the world instead of the
         */
        useWorldGravity: boolean;
        /**
         * Gravity vector, if the world's gravity is not used
         */
        readonly gravity: Vec2;
        Step(step: TimeStep): void;
        Draw(debugDraw: Draw): void;
    }
}
declare namespace b2 {
    /**
     * Applies a force every frame
     */
    class ConstantAccelController extends Controller {
        /**
         * The acceleration to apply
         */
        readonly A: Vec2;
        Step(step: TimeStep): void;
        private static Step_s_dtA;
        Draw(draw: Draw): void;
    }
}
declare namespace b2 {
    /**
     * Applies a force every frame
     */
    class ConstantForceController extends Controller {
        /**
         * The force to apply
         */
        readonly F: Vec2;
        Step(step: TimeStep): void;
        Draw(draw: Draw): void;
    }
}
declare namespace b2 {
    /**
     * Applies simplified gravity between every pair of bodies
     */
    class GravityController extends Controller {
        /**
         * Specifies the strength of the gravitiation force
         */
        G: number;
        /**
         * If true, gravity is proportional to r^-2, otherwise r^-1
         */
        invSqr: boolean;
        /**
         * @see Controller::Step
         */
        Step(step: TimeStep): void;
        private static Step_s_f;
        Draw(draw: Draw): void;
    }
}
declare namespace b2 {
    /**
     * Applies top down linear damping to the controlled bodies
     * The damping is calculated by multiplying velocity by a matrix
     * in local co-ordinates.
     */
    class TensorDampingController extends Controller {
        readonly T: Mat22;
        maxTimestep: number;
        /**
         * @see Controller::Step
         */
        Step(step: TimeStep): void;
        private static Step_s_damping;
        Draw(draw: Draw): void;
        /**
         * Sets damping independantly along the x and y axes
         */
        SetAxisAligned(xDamping: number, yDamping: number): void;
    }
}
declare namespace b2 {
    interface IAreaJointDef extends IJointDef {
        bodies: Body[];
        stiffness?: number;
        damping?: number;
    }
    class AreaJointDef extends JointDef implements IAreaJointDef {
        bodies: Body[];
        stiffness: number;
        damping: number;
        constructor();
        AddBody(body: Body): void;
    }
    class AreaJoint extends Joint {
        bodies: Body[];
        stiffness: number;
        damping: number;
        impulse: number;
        readonly targetLengths: number[];
        targetArea: number;
        readonly normals: Vec2[];
        readonly joints: DistanceJoint[];
        readonly deltas: Vec2[];
        readonly delta: Vec2;
        constructor(def: IAreaJointDef);
        GetAnchorA<T extends XY>(out: T): T;
        GetAnchorB<T extends XY>(out: T): T;
        GetReactionForce<T extends XY>(inv_dt: number, out: T): T;
        GetReactionTorque(inv_dt: number): number;
        SetStiffness(stiffness: number): void;
        GetStiffness(): number;
        SetDamping(damping: number): void;
        GetDamping(): number;
        Dump(log: (format: string, ...args: any[]) => void): void;
        InitVelocityConstraints(data: SolverData): void;
        SolveVelocityConstraints(data: SolverData): void;
        SolvePositionConstraints(data: SolverData): boolean;
    }
}
declare namespace b2 {
    enum BodyType {
        unknown = -1,
        staticBody = 0,
        kinematicBody = 1,
        dynamicBody = 2
    }
    interface IBodyDef {
        type?: BodyType;
        position?: XY;
        angle?: number;
        linearVelocity?: XY;
        angularVelocity?: number;
        linearDamping?: number;
        angularDamping?: number;
        allowSleep?: boolean;
        awake?: boolean;
        fixedRotation?: boolean;
        bullet?: boolean;
        enabled?: boolean;
        userData?: any;
        gravityScale?: number;
    }
    class BodyDef implements IBodyDef {
        type: BodyType;
        readonly position: Vec2;
        angle: number;
        readonly linearVelocity: Vec2;
        angularVelocity: number;
        linearDamping: number;
        angularDamping: number;
        allowSleep: boolean;
        awake: boolean;
        fixedRotation: boolean;
        bullet: boolean;
        enabled: boolean;
        userData: any;
        gravityScale: number;
    }
    class Body {
        type: BodyType;
        islandFlag: boolean;
        awakeFlag: boolean;
        autoSleepFlag: boolean;
        bulletFlag: boolean;
        fixedRotationFlag: boolean;
        enabledFlag: boolean;
        toiFlag: boolean;
        islandIndex: number;
        readonly xf: Transform;
        readonly xf0: Transform;
        readonly sweep: Sweep;
        readonly linearVelocity: Vec2;
        angularVelocity: number;
        readonly force: Vec2;
        torque: number;
        world: World;
        prev: Body | null;
        next: Body | null;
        fixtureList: Fixture | null;
        fixtureCount: number;
        jointList: JointEdge | null;
        contactList: ContactEdge | null;
        mass: number;
        invMass: number;
        I: number;
        invI: number;
        linearDamping: number;
        angularDamping: number;
        gravityScale: number;
        sleepTime: number;
        userData: any;
        controllerList: ControllerEdge | null;
        controllerCount: number;
        constructor(bd: IBodyDef, world: World);
        CreateFixture(def: IFixtureDef): Fixture;
        CreateFixture(shape: Shape): Fixture;
        CreateFixture(shape: Shape, density: number): Fixture;
        CreateFixtureDef(def: IFixtureDef): Fixture;
        private static CreateFixtureShapeDensity_s_def;
        CreateFixtureShapeDensity(shape: Shape, density?: number): Fixture;
        DestroyFixture(fixture: Fixture): void;
        SetTransformVec(position: XY, angle: number): void;
        SetTransformXY(x: number, y: number, angle: number): void;
        SetTransform(xf: Transform): void;
        GetTransform(): Transform;
        GetPosition(): Vec2;
        SetPosition(position: XY): void;
        SetPositionXY(x: number, y: number): void;
        GetAngle(): number;
        SetAngle(angle: number): void;
        GetWorldCenter(): Vec2;
        GetLocalCenter(): Vec2;
        SetLinearVelocity(v: XY): void;
        GetLinearVelocity(): Vec2;
        SetAngularVelocity(w: number): void;
        GetAngularVelocity(): number;
        GetDefinition(bd: BodyDef): BodyDef;
        ApplyForce(force: XY, point: XY, wake?: boolean): void;
        ApplyForceToCenter(force: XY, wake?: boolean): void;
        ApplyTorque(torque: number, wake?: boolean): void;
        ApplyLinearImpulse(impulse: XY, point: XY, wake?: boolean): void;
        ApplyLinearImpulseToCenter(impulse: XY, wake?: boolean): void;
        ApplyAngularImpulse(impulse: number, wake?: boolean): void;
        GetMass(): number;
        GetInertia(): number;
        GetMassData(data: MassData): MassData;
        private static SetMassData_s_oldCenter;
        SetMassData(massData: MassData): void;
        private static ResetMassData_s_localCenter;
        private static ResetMassData_s_oldCenter;
        private static ResetMassData_s_massData;
        ResetMassData(): void;
        GetWorldPoint<T extends XY>(localPoint: XY, out: T): T;
        GetWorldVector<T extends XY>(localVector: XY, out: T): T;
        GetLocalPoint<T extends XY>(worldPoint: XY, out: T): T;
        GetLocalVector<T extends XY>(worldVector: XY, out: T): T;
        GetLinearVelocityFromWorldPoint<T extends XY>(worldPoint: XY, out: T): T;
        GetLinearVelocityFromLocalPoint<T extends XY>(localPoint: XY, out: T): T;
        GetLinearDamping(): number;
        SetLinearDamping(linearDamping: number): void;
        GetAngularDamping(): number;
        SetAngularDamping(angularDamping: number): void;
        GetGravityScale(): number;
        SetGravityScale(scale: number): void;
        SetType(type: BodyType): void;
        GetType(): BodyType;
        SetBullet(flag: boolean): void;
        IsBullet(): boolean;
        SetSleepingAllowed(flag: boolean): void;
        IsSleepingAllowed(): boolean;
        SetAwake(flag: boolean): void;
        IsAwake(): boolean;
        SetEnabled(flag: boolean): void;
        IsEnabled(): boolean;
        SetFixedRotation(flag: boolean): void;
        IsFixedRotation(): boolean;
        GetFixtureList(): Fixture | null;
        GetJointList(): JointEdge | null;
        GetContactList(): ContactEdge | null;
        GetNext(): Body | null;
        GetUserData(): any;
        SetUserData(data: any): void;
        GetWorld(): World;
        Dump(log: (format: string, ...args: any[]) => void): void;
        private static SynchronizeFixtures_s_xf1;
        SynchronizeFixtures(): void;
        SynchronizeTransform(): void;
        ShouldCollide(other: Body): boolean;
        ShouldCollideConnected(other: Body): boolean;
        Advance(alpha: number): void;
        GetControllerList(): ControllerEdge | null;
        GetControllerCount(): number;
    }
}
declare namespace b2 {
    class ChainAndCircleContact extends Contact<ChainShape, CircleShape> {
        static Create(): Contact;
        static Destroy(contact: Contact): void;
        private static Evaluate_s_edge;
        Evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
    }
}
declare namespace b2 {
    class ChainAndPolygonContact extends Contact<ChainShape, PolygonShape> {
        static Create(): Contact;
        static Destroy(contact: Contact): void;
        private static Evaluate_s_edge;
        Evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
    }
}
declare namespace b2 {
    class CircleContact extends Contact<CircleShape, CircleShape> {
        static Create(): Contact;
        static Destroy(contact: Contact): void;
        Evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
    }
}
declare namespace b2 {
    class ContactRegister {
        pool: Contact[];
        createFcn: (() => Contact) | null;
        destroyFcn: ((contact: Contact) => void) | null;
        primary: boolean;
    }
    class ContactFactory {
        readonly registers: ContactRegister[][];
        constructor();
        private AddType;
        private InitializeRegisters;
        Create(fixtureA: Fixture, indexA: number, fixtureB: Fixture, indexB: number): Contact | null;
        Destroy(contact: Contact): void;
    }
}
declare namespace b2 {
    class ContactManager {
        readonly broadPhase: BroadPhase<FixtureProxy>;
        contactList: Contact | null;
        contactCount: number;
        contactFilter: ContactFilter;
        contactListener: ContactListener;
        readonly contactFactory: ContactFactory;
        AddPair(proxyA: FixtureProxy, proxyB: FixtureProxy): void;
        FindNewContacts(): void;
        Destroy(c: Contact): void;
        Collide(): void;
    }
}
declare namespace b2 {
    interface IDistanceJointDef extends IJointDef {
        localAnchorA?: XY;
        localAnchorB?: XY;
        length?: number;
        minLength?: number;
        maxLength?: number;
        stiffness?: number;
        damping?: number;
    }
    class DistanceJointDef extends JointDef implements IDistanceJointDef {
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        length: number;
        minLength: number;
        maxLength: number;
        stiffness: number;
        damping: number;
        constructor();
        Initialize(b1: Body, b2: Body, anchor1: XY, anchor2: XY): void;
    }
    class DistanceJoint extends Joint {
        stiffness: number;
        damping: number;
        bias: number;
        length: number;
        minLength: number;
        maxLength: number;
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        gamma: number;
        impulse: number;
        lowerImpulse: number;
        upperImpulse: number;
        indexA: number;
        indexB: number;
        readonly u: Vec2;
        readonly rA: Vec2;
        readonly rB: Vec2;
        readonly localCenterA: Vec2;
        readonly localCenterB: Vec2;
        currentLength: number;
        invMassA: number;
        invMassB: number;
        invIA: number;
        invIB: number;
        softMass: number;
        mass: number;
        readonly qA: Rot;
        readonly qB: Rot;
        readonly lalcA: Vec2;
        readonly lalcB: Vec2;
        constructor(def: IDistanceJointDef);
        GetAnchorA<T extends XY>(out: T): T;
        GetAnchorB<T extends XY>(out: T): T;
        GetReactionForce<T extends XY>(inv_dt: number, out: T): T;
        GetReactionTorque(inv_dt: number): number;
        GetLocalAnchorA(): Vec2;
        GetLocalAnchorB(): Vec2;
        SetLength(length: number): number;
        GetLength(): number;
        SetMinLength(minLength: number): number;
        SetMaxLength(maxLength: number): number;
        GetCurrentLength(): number;
        SetStiffness(stiffness: number): void;
        GetStiffness(): number;
        SetDamping(damping: number): void;
        GetDamping(): number;
        Dump(log: (format: string, ...args: any[]) => void): void;
        private static InitVelocityConstraints_s_P;
        InitVelocityConstraints(data: SolverData): void;
        private static SolveVelocityConstraints_s_vpA;
        private static SolveVelocityConstraints_s_vpB;
        private static SolveVelocityConstraints_s_P;
        SolveVelocityConstraints(data: SolverData): void;
        private static SolvePositionConstraints_s_P;
        SolvePositionConstraints(data: SolverData): boolean;
        private static Draw_s_pA;
        private static Draw_s_pB;
        private static Draw_s_axis;
        private static Draw_s_c1;
        private static Draw_s_c2;
        private static Draw_s_c3;
        private static Draw_s_c4;
        private static Draw_s_pRest;
        private static Draw_s_pMin;
        private static Draw_s_pMax;
        Draw(draw: Draw): void;
    }
}
declare namespace b2 {
    class EdgeAndCircleContact extends Contact<EdgeShape, CircleShape> {
        static Create(): Contact;
        static Destroy(contact: Contact): void;
        Evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
    }
}
declare namespace b2 {
    class EdgeAndPolygonContact extends Contact<EdgeShape, PolygonShape> {
        static Create(): Contact;
        static Destroy(contact: Contact): void;
        Evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
    }
}
declare namespace b2 {
    interface IFrictionJointDef extends IJointDef {
        localAnchorA?: XY;
        localAnchorB?: XY;
        maxForce?: number;
        maxTorque?: number;
    }
    class FrictionJointDef extends JointDef implements IFrictionJointDef {
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        maxForce: number;
        maxTorque: number;
        constructor();
        Initialize(bA: Body, bB: Body, anchor: Vec2): void;
    }
    class FrictionJoint extends Joint {
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        readonly linearImpulse: Vec2;
        angularImpulse: number;
        maxForce: number;
        maxTorque: number;
        indexA: number;
        indexB: number;
        readonly rA: Vec2;
        readonly rB: Vec2;
        readonly localCenterA: Vec2;
        readonly localCenterB: Vec2;
        invMassA: number;
        invMassB: number;
        invIA: number;
        invIB: number;
        readonly linearMass: Mat22;
        angularMass: number;
        readonly qA: Rot;
        readonly qB: Rot;
        readonly lalcA: Vec2;
        readonly lalcB: Vec2;
        readonly K: Mat22;
        constructor(def: IFrictionJointDef);
        InitVelocityConstraints(data: SolverData): void;
        private static SolveVelocityConstraints_s_Cdot_v2;
        private static SolveVelocityConstraints_s_impulseV;
        private static SolveVelocityConstraints_s_oldImpulseV;
        SolveVelocityConstraints(data: SolverData): void;
        SolvePositionConstraints(data: SolverData): boolean;
        GetAnchorA<T extends XY>(out: T): T;
        GetAnchorB<T extends XY>(out: T): T;
        GetReactionForce<T extends XY>(inv_dt: number, out: T): T;
        GetReactionTorque(inv_dt: number): number;
        GetLocalAnchorA(): Vec2;
        GetLocalAnchorB(): Vec2;
        SetMaxForce(force: number): void;
        GetMaxForce(): number;
        SetMaxTorque(torque: number): void;
        GetMaxTorque(): number;
        Dump(log: (format: string, ...args: any[]) => void): void;
    }
}
declare namespace b2 {
    interface IGearJointDef extends IJointDef {
        joint1: RevoluteJoint | PrismaticJoint;
        joint2: RevoluteJoint | PrismaticJoint;
        ratio?: number;
    }
    class GearJointDef extends JointDef implements IGearJointDef {
        joint1: RevoluteJoint | PrismaticJoint;
        joint2: RevoluteJoint | PrismaticJoint;
        ratio: number;
        constructor();
    }
    class GearJoint extends Joint {
        joint1: RevoluteJoint | PrismaticJoint;
        joint2: RevoluteJoint | PrismaticJoint;
        typeA: JointType;
        typeB: JointType;
        bodyC: Body;
        bodyD: Body;
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        readonly localAnchorC: Vec2;
        readonly localAnchorD: Vec2;
        readonly localAxisC: Vec2;
        readonly localAxisD: Vec2;
        referenceAngleA: number;
        referenceAngleB: number;
        constant: number;
        ratio: number;
        impulse: number;
        indexA: number;
        indexB: number;
        indexC: number;
        indexD: number;
        readonly lcA: Vec2;
        readonly lcB: Vec2;
        readonly lcC: Vec2;
        readonly lcD: Vec2;
        mA: number;
        mB: number;
        mC: number;
        mD: number;
        iA: number;
        iB: number;
        iC: number;
        iD: number;
        readonly JvAC: Vec2;
        readonly JvBD: Vec2;
        JwA: number;
        JwB: number;
        JwC: number;
        JwD: number;
        mass: number;
        readonly qA: Rot;
        readonly qB: Rot;
        readonly qC: Rot;
        readonly qD: Rot;
        readonly lalcA: Vec2;
        readonly lalcB: Vec2;
        readonly lalcC: Vec2;
        readonly lalcD: Vec2;
        constructor(def: IGearJointDef);
        private static InitVelocityConstraints_s_u;
        private static InitVelocityConstraints_s_rA;
        private static InitVelocityConstraints_s_rB;
        private static InitVelocityConstraints_s_rC;
        private static InitVelocityConstraints_s_rD;
        InitVelocityConstraints(data: SolverData): void;
        SolveVelocityConstraints(data: SolverData): void;
        private static SolvePositionConstraints_s_u;
        private static SolvePositionConstraints_s_rA;
        private static SolvePositionConstraints_s_rB;
        private static SolvePositionConstraints_s_rC;
        private static SolvePositionConstraints_s_rD;
        SolvePositionConstraints(data: SolverData): boolean;
        GetAnchorA<T extends XY>(out: T): T;
        GetAnchorB<T extends XY>(out: T): T;
        GetReactionForce<T extends XY>(inv_dt: number, out: T): T;
        GetReactionTorque(inv_dt: number): number;
        GetJoint1(): RevoluteJoint | PrismaticJoint;
        GetJoint2(): RevoluteJoint | PrismaticJoint;
        GetRatio(): number;
        SetRatio(ratio: number): void;
        Dump(log: (format: string, ...args: any[]) => void): void;
    }
}
declare namespace b2 {
    class Island {
        listener: ContactListener;
        readonly bodies: Body[];
        readonly contacts: Contact[];
        readonly joints: Joint[];
        readonly positions: Position[];
        readonly velocities: Velocity[];
        bodyCount: number;
        jointCount: number;
        contactCount: number;
        bodyCapacity: number;
        contactCapacity: number;
        jointCapacity: number;
        Initialize(bodyCapacity: number, contactCapacity: number, jointCapacity: number, listener: ContactListener): void;
        Clear(): void;
        AddBody(body: Body): void;
        AddContact(contact: Contact): void;
        AddJoint(joint: Joint): void;
        private static s_timer;
        private static s_solverData;
        private static s_contactSolverDef;
        private static s_contactSolver;
        private static s_translation;
        Solve(profile: Profile, step: TimeStep, gravity: Vec2, allowSleep: boolean): void;
        SolveTOI(subStep: TimeStep, toiIndexA: number, toiIndexB: number): void;
        private static s_impulse;
        Report(constraints: ContactVelocityConstraint[]): void;
    }
}
declare namespace b2 {
    interface IMotorJointDef extends IJointDef {
        linearOffset?: XY;
        angularOffset?: number;
        maxForce?: number;
        maxTorque?: number;
        correctionFactor?: number;
    }
    class MotorJointDef extends JointDef implements IMotorJointDef {
        readonly linearOffset: Vec2;
        angularOffset: number;
        maxForce: number;
        maxTorque: number;
        correctionFactor: number;
        constructor();
        Initialize(bA: Body, bB: Body): void;
    }
    class MotorJoint extends Joint {
        readonly linearOffset: Vec2;
        angularOffset: number;
        readonly linearImpulse: Vec2;
        angularImpulse: number;
        maxForce: number;
        maxTorque: number;
        correctionFactor: number;
        indexA: number;
        indexB: number;
        readonly rA: Vec2;
        readonly rB: Vec2;
        readonly localCenterA: Vec2;
        readonly localCenterB: Vec2;
        readonly linearError: Vec2;
        angularError: number;
        invMassA: number;
        invMassB: number;
        invIA: number;
        invIB: number;
        readonly linearMass: Mat22;
        angularMass: number;
        readonly qA: Rot;
        readonly qB: Rot;
        readonly K: Mat22;
        constructor(def: IMotorJointDef);
        GetAnchorA<T extends XY>(out: T): T;
        GetAnchorB<T extends XY>(out: T): T;
        GetReactionForce<T extends XY>(inv_dt: number, out: T): T;
        GetReactionTorque(inv_dt: number): number;
        SetLinearOffset(linearOffset: Vec2): void;
        GetLinearOffset(): Vec2;
        SetAngularOffset(angularOffset: number): void;
        GetAngularOffset(): number;
        SetMaxForce(force: number): void;
        GetMaxForce(): number;
        SetMaxTorque(torque: number): void;
        GetMaxTorque(): number;
        InitVelocityConstraints(data: SolverData): void;
        private static SolveVelocityConstraints_s_Cdot_v2;
        private static SolveVelocityConstraints_s_impulse_v2;
        private static SolveVelocityConstraints_s_oldImpulse_v2;
        SolveVelocityConstraints(data: SolverData): void;
        SolvePositionConstraints(data: SolverData): boolean;
        Dump(log: (format: string, ...args: any[]) => void): void;
    }
}
declare namespace b2 {
    interface IMouseJointDef extends IJointDef {
        target?: XY;
        maxForce?: number;
        stiffness?: number;
        damping?: number;
    }
    class MouseJointDef extends JointDef implements IMouseJointDef {
        readonly target: Vec2;
        maxForce: number;
        stiffness: number;
        damping: number;
        constructor();
    }
    class MouseJoint extends Joint {
        readonly localAnchorB: Vec2;
        readonly targetA: Vec2;
        stiffness: number;
        damping: number;
        beta: number;
        readonly impulse: Vec2;
        maxForce: number;
        gamma: number;
        indexA: number;
        indexB: number;
        readonly rB: Vec2;
        readonly localCenterB: Vec2;
        invMassB: number;
        invIB: number;
        readonly mass: Mat22;
        readonly C: Vec2;
        readonly qB: Rot;
        readonly lalcB: Vec2;
        readonly K: Mat22;
        constructor(def: IMouseJointDef);
        SetTarget(target: Vec2): void;
        GetTarget(): Vec2;
        SetMaxForce(maxForce: number): void;
        GetMaxForce(): number;
        SetStiffness(stiffness: number): void;
        GetStiffness(): number;
        SetDamping(damping: number): void;
        GetDamping(): number;
        InitVelocityConstraints(data: SolverData): void;
        private static SolveVelocityConstraints_s_Cdot;
        private static SolveVelocityConstraints_s_impulse;
        private static SolveVelocityConstraints_s_oldImpulse;
        SolveVelocityConstraints(data: SolverData): void;
        SolvePositionConstraints(data: SolverData): boolean;
        GetAnchorA<T extends XY>(out: T): T;
        GetAnchorB<T extends XY>(out: T): T;
        GetReactionForce<T extends XY>(inv_dt: number, out: T): T;
        GetReactionTorque(inv_dt: number): number;
        Dump(log: (format: string, ...args: any[]) => void): void;
        ShiftOrigin(newOrigin: Vec2): void;
    }
}
declare namespace b2 {
    class PolygonAndCircleContact extends Contact<PolygonShape, CircleShape> {
        static Create(): Contact;
        static Destroy(contact: Contact): void;
        Evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
    }
}
declare namespace b2 {
    class PolygonContact extends Contact<PolygonShape, PolygonShape> {
        static Create(): Contact;
        static Destroy(contact: Contact): void;
        Evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
    }
}
declare namespace b2 {
    interface IPrismaticJointDef extends IJointDef {
        localAnchorA?: XY;
        localAnchorB?: XY;
        localAxisA?: XY;
        referenceAngle?: number;
        enableLimit?: boolean;
        lowerTranslation?: number;
        upperTranslation?: number;
        enableMotor?: boolean;
        maxMotorForce?: number;
        motorSpeed?: number;
    }
    class PrismaticJointDef extends JointDef implements IPrismaticJointDef {
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        readonly localAxisA: Vec2;
        referenceAngle: number;
        enableLimit: boolean;
        lowerTranslation: number;
        upperTranslation: number;
        enableMotor: boolean;
        maxMotorForce: number;
        motorSpeed: number;
        constructor();
        Initialize(bA: Body, bB: Body, anchor: Vec2, axis: Vec2): void;
    }
    class PrismaticJoint extends Joint {
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        readonly localXAxisA: Vec2;
        readonly localYAxisA: Vec2;
        referenceAngle: number;
        readonly impulse: Vec2;
        motorImpulse: number;
        lowerImpulse: number;
        upperImpulse: number;
        lowerTranslation: number;
        upperTranslation: number;
        maxMotorForce: number;
        motorSpeed: number;
        enableLimit: boolean;
        enableMotor: boolean;
        indexA: number;
        indexB: number;
        readonly localCenterA: Vec2;
        readonly localCenterB: Vec2;
        invMassA: number;
        invMassB: number;
        invIA: number;
        invIB: number;
        readonly axis: Vec2;
        readonly perp: Vec2;
        s1: number;
        s2: number;
        a1: number;
        a2: number;
        readonly K: Mat22;
        readonly K3: Mat33;
        readonly K2: Mat22;
        translation: number;
        axialMass: number;
        readonly qA: Rot;
        readonly qB: Rot;
        readonly lalcA: Vec2;
        readonly lalcB: Vec2;
        readonly rA: Vec2;
        readonly rB: Vec2;
        constructor(def: IPrismaticJointDef);
        private static InitVelocityConstraints_s_d;
        private static InitVelocityConstraints_s_P;
        InitVelocityConstraints(data: SolverData): void;
        private static SolveVelocityConstraints_s_P;
        private static SolveVelocityConstraints_s_df;
        SolveVelocityConstraints(data: SolverData): void;
        private static SolvePositionConstraints_s_d;
        private static SolvePositionConstraints_s_impulse;
        private static SolvePositionConstraints_s_impulse1;
        private static SolvePositionConstraints_s_P;
        SolvePositionConstraints(data: SolverData): boolean;
        GetAnchorA<T extends XY>(out: T): T;
        GetAnchorB<T extends XY>(out: T): T;
        GetReactionForce<T extends XY>(inv_dt: number, out: T): T;
        GetReactionTorque(inv_dt: number): number;
        GetLocalAnchorA(): Vec2;
        GetLocalAnchorB(): Vec2;
        GetLocalAxisA(): Vec2;
        GetReferenceAngle(): number;
        private static GetJointTranslation_s_pA;
        private static GetJointTranslation_s_pB;
        private static GetJointTranslation_s_d;
        private static GetJointTranslation_s_axis;
        GetJointTranslation(): number;
        GetJointSpeed(): number;
        IsLimitEnabled(): boolean;
        EnableLimit(flag: boolean): void;
        GetLowerLimit(): number;
        GetUpperLimit(): number;
        SetLimits(lower: number, upper: number): void;
        IsMotorEnabled(): boolean;
        EnableMotor(flag: boolean): void;
        SetMotorSpeed(speed: number): void;
        GetMotorSpeed(): number;
        SetMaxMotorForce(force: number): void;
        GetMaxMotorForce(): number;
        GetMotorForce(inv_dt: number): number;
        Dump(log: (format: string, ...args: any[]) => void): void;
        private static Draw_s_pA;
        private static Draw_s_pB;
        private static Draw_s_axis;
        private static Draw_s_c1;
        private static Draw_s_c2;
        private static Draw_s_c3;
        private static Draw_s_c4;
        private static Draw_s_c5;
        private static Draw_s_lower;
        private static Draw_s_upper;
        private static Draw_s_perp;
        Draw(draw: Draw): void;
    }
}
declare namespace b2 {
    const minPulleyLength: number;
    interface IPulleyJointDef extends IJointDef {
        groundAnchorA?: XY;
        groundAnchorB?: XY;
        localAnchorA?: XY;
        localAnchorB?: XY;
        lengthA?: number;
        lengthB?: number;
        ratio?: number;
    }
    class PulleyJointDef extends JointDef implements IPulleyJointDef {
        readonly groundAnchorA: Vec2;
        readonly groundAnchorB: Vec2;
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        lengthA: number;
        lengthB: number;
        ratio: number;
        constructor();
        Initialize(bA: Body, bB: Body, groundA: Vec2, groundB: Vec2, anchorA: Vec2, anchorB: Vec2, r: number): void;
    }
    class PulleyJoint extends Joint {
        readonly groundAnchorA: Vec2;
        readonly groundAnchorB: Vec2;
        lengthA: number;
        lengthB: number;
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        constant: number;
        ratio: number;
        impulse: number;
        indexA: number;
        indexB: number;
        readonly uA: Vec2;
        readonly uB: Vec2;
        readonly rA: Vec2;
        readonly rB: Vec2;
        readonly localCenterA: Vec2;
        readonly localCenterB: Vec2;
        invMassA: number;
        invMassB: number;
        invIA: number;
        invIB: number;
        mass: number;
        readonly qA: Rot;
        readonly qB: Rot;
        readonly lalcA: Vec2;
        readonly lalcB: Vec2;
        constructor(def: IPulleyJointDef);
        private static InitVelocityConstraints_s_PA;
        private static InitVelocityConstraints_s_PB;
        InitVelocityConstraints(data: SolverData): void;
        private static SolveVelocityConstraints_s_vpA;
        private static SolveVelocityConstraints_s_vpB;
        private static SolveVelocityConstraints_s_PA;
        private static SolveVelocityConstraints_s_PB;
        SolveVelocityConstraints(data: SolverData): void;
        private static SolvePositionConstraints_s_PA;
        private static SolvePositionConstraints_s_PB;
        SolvePositionConstraints(data: SolverData): boolean;
        GetAnchorA<T extends XY>(out: T): T;
        GetAnchorB<T extends XY>(out: T): T;
        GetReactionForce<T extends XY>(inv_dt: number, out: T): T;
        GetReactionTorque(inv_dt: number): number;
        GetGroundAnchorA(): Vec2;
        GetGroundAnchorB(): Vec2;
        GetLengthA(): number;
        GetLengthB(): number;
        GetRatio(): number;
        private static GetCurrentLengthA_s_p;
        GetCurrentLengthA(): number;
        private static GetCurrentLengthB_s_p;
        GetCurrentLengthB(): number;
        Dump(log: (format: string, ...args: any[]) => void): void;
        ShiftOrigin(newOrigin: Vec2): void;
    }
}
declare namespace b2 {
    interface IRevoluteJointDef extends IJointDef {
        localAnchorA?: XY;
        localAnchorB?: XY;
        referenceAngle?: number;
        enableLimit?: boolean;
        lowerAngle?: number;
        upperAngle?: number;
        enableMotor?: boolean;
        motorSpeed?: number;
        maxMotorTorque?: number;
    }
    class RevoluteJointDef extends JointDef implements IRevoluteJointDef {
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        referenceAngle: number;
        enableLimit: boolean;
        lowerAngle: number;
        upperAngle: number;
        enableMotor: boolean;
        motorSpeed: number;
        maxMotorTorque: number;
        constructor();
        Initialize(bA: Body, bB: Body, anchor: XY): void;
    }
    class RevoluteJoint extends Joint {
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        readonly impulse: Vec2;
        motorImpulse: number;
        lowerImpulse: number;
        upperImpulse: number;
        enableMotor: boolean;
        maxMotorTorque: number;
        motorSpeed: number;
        enableLimit: boolean;
        referenceAngle: number;
        lowerAngle: number;
        upperAngle: number;
        indexA: number;
        indexB: number;
        readonly rA: Vec2;
        readonly rB: Vec2;
        readonly localCenterA: Vec2;
        readonly localCenterB: Vec2;
        invMassA: number;
        invMassB: number;
        invIA: number;
        invIB: number;
        readonly K: Mat22;
        angle: number;
        axialMass: number;
        readonly qA: Rot;
        readonly qB: Rot;
        readonly lalcA: Vec2;
        readonly lalcB: Vec2;
        constructor(def: IRevoluteJointDef);
        private static InitVelocityConstraints_s_P;
        InitVelocityConstraints(data: SolverData): void;
        private static SolveVelocityConstraints_s_Cdot_v2;
        private static SolveVelocityConstraints_s_impulse_v2;
        SolveVelocityConstraints(data: SolverData): void;
        private static SolvePositionConstraints_s_C_v2;
        private static SolvePositionConstraints_s_impulse;
        SolvePositionConstraints(data: SolverData): boolean;
        GetAnchorA<T extends XY>(out: T): T;
        GetAnchorB<T extends XY>(out: T): T;
        GetReactionForce<T extends XY>(inv_dt: number, out: T): T;
        GetReactionTorque(inv_dt: number): number;
        GetLocalAnchorA(): Vec2;
        GetLocalAnchorB(): Vec2;
        GetReferenceAngle(): number;
        GetJointAngle(): number;
        GetJointSpeed(): number;
        IsMotorEnabled(): boolean;
        EnableMotor(flag: boolean): void;
        GetMotorTorque(inv_dt: number): number;
        GetMotorSpeed(): number;
        SetMaxMotorTorque(torque: number): void;
        GetMaxMotorTorque(): number;
        IsLimitEnabled(): boolean;
        EnableLimit(flag: boolean): void;
        GetLowerLimit(): number;
        GetUpperLimit(): number;
        SetLimits(lower: number, upper: number): void;
        SetMotorSpeed(speed: number): void;
        Dump(log: (format: string, ...args: any[]) => void): void;
        private static Draw_s_pA;
        private static Draw_s_pB;
        private static Draw_s_c1;
        private static Draw_s_c2;
        private static Draw_s_c3;
        private static Draw_s_c4;
        private static Draw_s_c5;
        private static Draw_s_color_;
        private static Draw_s_r;
        private static Draw_s_rlo;
        private static Draw_s_rhi;
        Draw(draw: Draw): void;
    }
}
declare namespace b2 {
    interface IWeldJointDef extends IJointDef {
        localAnchorA?: XY;
        localAnchorB?: XY;
        referenceAngle?: number;
        stiffness?: number;
        damping?: number;
    }
    class WeldJointDef extends JointDef implements IWeldJointDef {
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        referenceAngle: number;
        stiffness: number;
        damping: number;
        constructor();
        Initialize(bA: Body, bB: Body, anchor: Vec2): void;
    }
    class WeldJoint extends Joint {
        stiffness: number;
        damping: number;
        bias: number;
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        referenceAngle: number;
        gamma: number;
        readonly impulse: Vec3;
        indexA: number;
        indexB: number;
        readonly rA: Vec2;
        readonly rB: Vec2;
        readonly localCenterA: Vec2;
        readonly localCenterB: Vec2;
        invMassA: number;
        invMassB: number;
        invIA: number;
        invIB: number;
        readonly mass: Mat33;
        readonly qA: Rot;
        readonly qB: Rot;
        readonly lalcA: Vec2;
        readonly lalcB: Vec2;
        readonly K: Mat33;
        constructor(def: IWeldJointDef);
        private static InitVelocityConstraints_s_P;
        InitVelocityConstraints(data: SolverData): void;
        private static SolveVelocityConstraints_s_Cdot1;
        private static SolveVelocityConstraints_s_impulse1;
        private static SolveVelocityConstraints_s_impulse;
        private static SolveVelocityConstraints_s_P;
        SolveVelocityConstraints(data: SolverData): void;
        private static SolvePositionConstraints_s_C1;
        private static SolvePositionConstraints_s_P;
        private static SolvePositionConstraints_s_impulse;
        SolvePositionConstraints(data: SolverData): boolean;
        GetAnchorA<T extends XY>(out: T): T;
        GetAnchorB<T extends XY>(out: T): T;
        GetReactionForce<T extends XY>(inv_dt: number, out: T): T;
        GetReactionTorque(inv_dt: number): number;
        GetLocalAnchorA(): Vec2;
        GetLocalAnchorB(): Vec2;
        GetReferenceAngle(): number;
        SetStiffness(stiffness: number): void;
        GetStiffness(): number;
        SetDamping(damping: number): void;
        GetDamping(): number;
        Dump(log: (format: string, ...args: any[]) => void): void;
    }
}
declare namespace b2 {
    interface IWheelJointDef extends IJointDef {
        localAnchorA?: XY;
        localAnchorB?: XY;
        localAxisA?: XY;
        enableLimit?: boolean;
        lowerTranslation?: number;
        upperTranslation?: number;
        enableMotor?: boolean;
        maxMotorTorque?: number;
        motorSpeed?: number;
        stiffness?: number;
        damping?: number;
    }
    class WheelJointDef extends JointDef implements IWheelJointDef {
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        readonly localAxisA: Vec2;
        enableLimit: boolean;
        lowerTranslation: number;
        upperTranslation: number;
        enableMotor: boolean;
        maxMotorTorque: number;
        motorSpeed: number;
        stiffness: number;
        damping: number;
        constructor();
        Initialize(bA: Body, bB: Body, anchor: Vec2, axis: Vec2): void;
    }
    class WheelJoint extends Joint {
        readonly localAnchorA: Vec2;
        readonly localAnchorB: Vec2;
        readonly localXAxisA: Vec2;
        readonly localYAxisA: Vec2;
        impulse: number;
        motorImpulse: number;
        springImpulse: number;
        lowerImpulse: number;
        upperImpulse: number;
        translation: number;
        lowerTranslation: number;
        upperTranslation: number;
        maxMotorTorque: number;
        motorSpeed: number;
        enableLimit: boolean;
        enableMotor: boolean;
        stiffness: number;
        damping: number;
        indexA: number;
        indexB: number;
        readonly localCenterA: Vec2;
        readonly localCenterB: Vec2;
        invMassA: number;
        invMassB: number;
        invIA: number;
        invIB: number;
        readonly ax: Vec2;
        readonly ay: Vec2;
        sAx: number;
        sBx: number;
        sAy: number;
        sBy: number;
        mass: number;
        motorMass: number;
        axialMass: number;
        springMass: number;
        bias: number;
        gamma: number;
        readonly qA: Rot;
        readonly qB: Rot;
        readonly lalcA: Vec2;
        readonly lalcB: Vec2;
        readonly rA: Vec2;
        readonly rB: Vec2;
        constructor(def: IWheelJointDef);
        GetMotorSpeed(): number;
        GetMaxMotorTorque(): number;
        SetSpringFrequencyHz(hz: number): void;
        GetSpringFrequencyHz(): number;
        SetSpringDampingRatio(ratio: number): void;
        GetSpringDampingRatio(): number;
        private static InitVelocityConstraints_s_d;
        private static InitVelocityConstraints_s_P;
        InitVelocityConstraints(data: SolverData): void;
        private static SolveVelocityConstraints_s_P;
        SolveVelocityConstraints(data: SolverData): void;
        private static SolvePositionConstraints_s_d;
        private static SolvePositionConstraints_s_P;
        SolvePositionConstraints(data: SolverData): boolean;
        GetDefinition(def: WheelJointDef): WheelJointDef;
        GetAnchorA<T extends XY>(out: T): T;
        GetAnchorB<T extends XY>(out: T): T;
        GetReactionForce<T extends XY>(inv_dt: number, out: T): T;
        GetReactionTorque(inv_dt: number): number;
        GetLocalAnchorA(): Vec2;
        GetLocalAnchorB(): Vec2;
        GetLocalAxisA(): Vec2;
        GetJointTranslation(): number;
        GetJointLinearSpeed(): number;
        GetJointAngle(): number;
        GetJointAngularSpeed(): number;
        GetPrismaticJointTranslation(): number;
        GetPrismaticJointSpeed(): number;
        GetRevoluteJointAngle(): number;
        GetRevoluteJointSpeed(): number;
        IsMotorEnabled(): boolean;
        EnableMotor(flag: boolean): void;
        SetMotorSpeed(speed: number): void;
        SetMaxMotorTorque(force: number): void;
        GetMotorTorque(inv_dt: number): number;
        IsLimitEnabled(): boolean;
        EnableLimit(flag: boolean): void;
        GetLowerLimit(): number;
        GetUpperLimit(): number;
        SetLimits(lower: number, upper: number): void;
        Dump(log: (format: string, ...args: any[]) => void): void;
        private static Draw_s_pA;
        private static Draw_s_pB;
        private static Draw_s_axis;
        private static Draw_s_c1;
        private static Draw_s_c2;
        private static Draw_s_c3;
        private static Draw_s_c4;
        private static Draw_s_c5;
        private static Draw_s_lower;
        private static Draw_s_upper;
        private static Draw_s_perp;
        Draw(draw: Draw): void;
    }
}
declare namespace b2 {
    class World {
        readonly contactManager: ContactManager;
        bodyList: Body | null;
        jointList: Joint | null;
        particleSystemList: ParticleSystem | null;
        bodyCount: number;
        jointCount: number;
        readonly gravity: Vec2;
        allowSleep: boolean;
        destructionListener: DestructionListener | null;
        debugDraw: Draw;
        inv_dt0: number;
        newContacts: boolean;
        locked: boolean;
        clearForces: boolean;
        warmStarting: boolean;
        continuousPhysics: boolean;
        subStepping: boolean;
        stepComplete: boolean;
        readonly profile: Profile;
        readonly island: Island;
        readonly s_stack: Array<Body | null>;
        controllerList: Controller | null;
        controllerCount: number;
        constructor(gravity: XY);
        SetDestructionListener(listener: DestructionListener | null): void;
        SetContactFilter(filter: ContactFilter): void;
        SetContactListener(listener: ContactListener): void;
        SetDebugDraw(debugDraw: Draw): void;
        CreateBody(def?: IBodyDef): Body;
        DestroyBody(b: Body): void;
        private static _Joint_Create;
        private static _Joint_Destroy;
        CreateJoint(def: IAreaJointDef): AreaJoint;
        CreateJoint(def: IDistanceJointDef): DistanceJoint;
        CreateJoint(def: IFrictionJointDef): FrictionJoint;
        CreateJoint(def: IGearJointDef): GearJoint;
        CreateJoint(def: IMotorJointDef): MotorJoint;
        CreateJoint(def: IMouseJointDef): MouseJoint;
        CreateJoint(def: IPrismaticJointDef): PrismaticJoint;
        CreateJoint(def: IPulleyJointDef): PulleyJoint;
        CreateJoint(def: IRevoluteJointDef): RevoluteJoint;
        CreateJoint(def: IWeldJointDef): WeldJoint;
        CreateJoint(def: IWheelJointDef): WheelJoint;
        DestroyJoint(j: Joint): void;
        CreateParticleSystem(def: ParticleSystemDef): ParticleSystem;
        DestroyParticleSystem(p: ParticleSystem): void;
        CalculateReasonableParticleIterations(timeStep: number): number;
        private static Step_s_step;
        private static Step_s_stepTimer;
        private static Step_s_timer;
        Step(dt: number, velocityIterations: number, positionIterations: number, particleIterations?: number): void;
        ClearForces(): void;
        DrawParticleSystem(system: ParticleSystem): void;
        private static DebugDraw_s_color;
        private static DebugDraw_s_vs;
        private static DebugDraw_s_xf;
        DebugDraw(): void;
        QueryAABB(callback: QueryCallback, aabb: AABB): void;
        QueryAABB(aabb: AABB, fn: QueryCallbackFunction): void;
        private _QueryAABB;
        QueryAllAABB(aabb: AABB, out?: Fixture[]): Fixture[];
        QueryPointAABB(callback: QueryCallback, point: XY): void;
        QueryPointAABB(point: XY, fn: QueryCallbackFunction): void;
        private _QueryPointAABB;
        QueryAllPointAABB(point: XY, out?: Fixture[]): Fixture[];
        QueryFixtureShape(callback: QueryCallback, shape: Shape, index: number, transform: Transform): void;
        QueryFixtureShape(shape: Shape, index: number, transform: Transform, fn: QueryCallbackFunction): void;
        private static QueryFixtureShape_s_aabb;
        private _QueryFixtureShape;
        QueryAllFixtureShape(shape: Shape, index: number, transform: Transform, out?: Fixture[]): Fixture[];
        QueryFixturePoint(callback: QueryCallback, point: XY): void;
        QueryFixturePoint(point: XY, fn: QueryCallbackFunction): void;
        private _QueryFixturePoint;
        QueryAllFixturePoint(point: XY, out?: Fixture[]): Fixture[];
        RayCast(callback: RayCastCallback, point1: XY, point2: XY): void;
        RayCast(point1: XY, point2: XY, fn: RayCastCallbackFunction): void;
        private static RayCast_s_input;
        private static RayCast_s_output;
        private static RayCast_s_point;
        private _RayCast;
        RayCastOne(point1: XY, point2: XY): Fixture | null;
        RayCastAll(point1: XY, point2: XY, out?: Fixture[]): Fixture[];
        GetBodyList(): Body | null;
        GetJointList(): Joint | null;
        GetParticleSystemList(): ParticleSystem | null;
        GetContactList(): Contact | null;
        SetAllowSleeping(flag: boolean): void;
        GetAllowSleeping(): boolean;
        SetWarmStarting(flag: boolean): void;
        GetWarmStarting(): boolean;
        SetContinuousPhysics(flag: boolean): void;
        GetContinuousPhysics(): boolean;
        SetSubStepping(flag: boolean): void;
        GetSubStepping(): boolean;
        GetProxyCount(): number;
        GetBodyCount(): number;
        GetJointCount(): number;
        GetContactCount(): number;
        GetTreeHeight(): number;
        GetTreeBalance(): number;
        GetTreeQuality(): number;
        SetGravity(gravity: XY, wake?: boolean): void;
        GetGravity(): Vec2;
        IsLocked(): boolean;
        SetAutoClearForces(flag: boolean): void;
        GetAutoClearForces(): boolean;
        ShiftOrigin(newOrigin: XY): void;
        GetContactManager(): ContactManager;
        GetProfile(): Profile;
        Dump(log: (format: string, ...args: any[]) => void): void;
        DrawShape(fixture: Fixture, color: Color): void;
        Solve(step: TimeStep): void;
        private static SolveTOI_s_subStep;
        private static SolveTOI_s_backup;
        private static SolveTOI_s_backup1;
        private static SolveTOI_s_backup2;
        private static SolveTOI_s_toi_input;
        private static SolveTOI_s_toi_output;
        SolveTOI(step: TimeStep): void;
        AddController(controller: Controller): Controller;
        RemoveController(controller: Controller): Controller;
    }
}
declare namespace b2 {
    enum ParticleGroupFlag {
        solidParticleGroup = 1,
        rigidParticleGroup = 2,
        particleGroupCanBeEmpty = 4,
        particleGroupWillBeDestroyed = 8,
        particleGroupNeedsUpdateDepth = 16,
        particleGroupInternalMask = 24
    }
    interface IParticleGroupDef {
        flags?: ParticleFlag;
        groupFlags?: ParticleGroupFlag;
        position?: XY;
        angle?: number;
        linearVelocity?: XY;
        angularVelocity?: number;
        color?: RGBA;
        strength?: number;
        shape?: Shape;
        shapes?: Shape[];
        shapeCount?: number;
        stride?: number;
        particleCount?: number;
        positionData?: XY[];
        lifetime?: number;
        userData?: any;
        group?: ParticleGroup | null;
    }
    class ParticleGroupDef implements IParticleGroupDef {
        flags: ParticleFlag;
        groupFlags: ParticleGroupFlag;
        readonly position: Vec2;
        angle: number;
        readonly linearVelocity: Vec2;
        angularVelocity: number;
        readonly color: Color;
        strength: number;
        shape?: Shape;
        shapes?: Shape[];
        shapeCount: number;
        stride: number;
        particleCount: number;
        positionData?: Vec2[];
        lifetime: number;
        userData: any;
        group: ParticleGroup | null;
    }
    class ParticleGroup {
        readonly system: ParticleSystem;
        firstIndex: number;
        lastIndex: number;
        groupFlags: ParticleGroupFlag;
        strength: number;
        prev: ParticleGroup | null;
        next: ParticleGroup | null;
        timestamp: number;
        mass: number;
        inertia: number;
        readonly center: Vec2;
        readonly linearVelocity: Vec2;
        angularVelocity: number;
        readonly transform: Transform;
        userData: any;
        constructor(system: ParticleSystem);
        GetNext(): ParticleGroup | null;
        GetParticleSystem(): ParticleSystem;
        GetParticleCount(): number;
        GetBufferIndex(): number;
        ContainsParticle(index: number): boolean;
        GetAllParticleFlags(): ParticleFlag;
        GetGroupFlags(): ParticleGroupFlag;
        SetGroupFlags(flags: number): void;
        GetMass(): number;
        GetInertia(): number;
        GetCenter(): Vec2;
        GetLinearVelocity(): Vec2;
        GetAngularVelocity(): number;
        GetTransform(): Transform;
        GetPosition(): Vec2;
        GetAngle(): number;
        GetLinearVelocityFromWorldPoint<T extends XY>(worldPoint: XY, out: T): T;
        static readonly GetLinearVelocityFromWorldPoint_s_t0: Vec2;
        GetUserData(): void;
        SetUserData(data: any): void;
        ApplyForce(force: XY): void;
        ApplyLinearImpulse(impulse: XY): void;
        DestroyParticles(callDestructionListener: boolean): void;
        UpdateStatistics(): void;
    }
}
declare namespace b2 {
    class GrowableBuffer<T> {
        data: T[];
        count: number;
        capacity: number;
        allocator: () => T;
        constructor(allocator: () => T);
        Append(): number;
        Reserve(newCapacity: number): void;
        Grow(): void;
        Free(): void;
        Shorten(newEnd: number): void;
        Data(): T[];
        GetCount(): number;
        SetCount(newCount: number): void;
        GetCapacity(): number;
        RemoveIf(pred: (t: T) => boolean): void;
        Unique(pred: (a: T, b: T) => boolean): void;
    }
    type ParticleIndex = number;
    class FixtureParticleQueryCallback extends QueryCallback {
        system: ParticleSystem;
        constructor(system: ParticleSystem);
        ShouldQueryParticleSystem(system: ParticleSystem): boolean;
        ReportFixture(fixture: Fixture): boolean;
        ReportParticle(system: ParticleSystem, index: number): boolean;
        ReportFixtureAndParticle(fixture: Fixture, childIndex: number, index: number): void;
    }
    class ParticleContact {
        indexA: number;
        indexB: number;
        weight: number;
        normal: Vec2;
        flags: ParticleFlag;
        SetIndices(a: number, b: number): void;
        SetWeight(w: number): void;
        SetNormal(n: Vec2): void;
        SetFlags(f: ParticleFlag): void;
        GetIndexA(): number;
        GetIndexB(): number;
        GetWeight(): number;
        GetNormal(): Vec2;
        GetFlags(): ParticleFlag;
        IsEqual(rhs: ParticleContact): boolean;
        IsNotEqual(rhs: ParticleContact): boolean;
        ApproximatelyEqual(rhs: ParticleContact): boolean;
    }
    class ParticleBodyContact {
        index: number;
        body: Body;
        fixture: Fixture;
        weight: number;
        normal: Vec2;
        mass: number;
    }
    class ParticlePair {
        indexA: number;
        indexB: number;
        flags: ParticleFlag;
        strength: number;
        distance: number;
    }
    class ParticleTriad {
        indexA: number;
        indexB: number;
        indexC: number;
        flags: ParticleFlag;
        strength: number;
        pa: Vec2;
        pb: Vec2;
        pc: Vec2;
        ka: number;
        kb: number;
        kc: number;
        s: number;
    }
    class ParticleSystemDef {
        /**
         * Enable strict Particle/Body contact check.
         * See SetStrictContactCheck for details.
         */
        strictContactCheck: boolean;
        /**
         * Set the particle density.
         * See SetDensity for details.
         */
        density: number;
        /**
         * Change the particle gravity scale. Adjusts the effect of the
         * global gravity vector on particles. Default value is 1.0f.
         */
        gravityScale: number;
        /**
         * Particles behave as circles with this radius. In Box2D units.
         */
        radius: number;
        /**
         * Set the maximum number of particles.
         * By default, there is no maximum. The particle buffers can
         * continue to grow while World's block allocator still has
         * memory.
         * See SetMaxParticleCount for details.
         */
        maxCount: number;
        /**
         * Increases pressure in response to compression
         * Smaller values allow more compression
         */
        pressureStrength: number;
        /**
         * Reduces velocity along the collision normal
         * Smaller value reduces less
         */
        dampingStrength: number;
        /**
         * Restores shape of elastic particle groups
         * Larger values increase elastic particle velocity
         */
        elasticStrength: number;
        /**
         * Restores length of spring particle groups
         * Larger values increase spring particle velocity
         */
        springStrength: number;
        /**
         * Reduces relative velocity of viscous particles
         * Larger values slow down viscous particles more
         */
        viscousStrength: number;
        /**
         * Produces pressure on tensile particles
         * 0~0.2. Larger values increase the amount of surface tension.
         */
        surfaceTensionPressureStrength: number;
        /**
         * Smoothes outline of tensile particles
         * 0~0.2. Larger values result in rounder, smoother,
         * water-drop-like clusters of particles.
         */
        surfaceTensionNormalStrength: number;
        /**
         * Produces additional pressure on repulsive particles
         * Larger values repulse more
         * Negative values mean attraction. The range where particles
         * behave stably is about -0.2 to 2.0.
         */
        repulsiveStrength: number;
        /**
         * Produces repulsion between powder particles
         * Larger values repulse more
         */
        powderStrength: number;
        /**
         * Pushes particles out of solid particle group
         * Larger values repulse more
         */
        ejectionStrength: number;
        /**
         * Produces static pressure
         * Larger values increase the pressure on neighboring partilces
         * For a description of static pressure, see
         * http://en.wikipedia.org/wiki/Static_pressure#Static_pressure_in_fluid_dynamics
         */
        staticPressureStrength: number;
        /**
         * Reduces instability in static pressure calculation
         * Larger values make stabilize static pressure with fewer
         * iterations
         */
        staticPressureRelaxation: number;
        /**
         * Computes static pressure more precisely
         * See SetStaticPressureIterations for details
         */
        staticPressureIterations: number;
        /**
         * Determines how fast colors are mixed
         * 1.0f ==> mixed immediately
         * 0.5f ==> mixed half way each simulation step (see
         * World::Step())
         */
        colorMixingStrength: number;
        /**
         * Whether to destroy particles by age when no more particles
         * can be created.  See #ParticleSystem::SetDestructionByAge()
         * for more information.
         */
        destroyByAge: boolean;
        /**
         * Granularity of particle lifetimes in seconds.  By default
         * this is set to (1.0f / 60.0f) seconds.  ParticleSystem uses
         * a 32-bit signed value to track particle lifetimes so the
         * maximum lifetime of a particle is (2^32 - 1) / (1.0f /
         * lifetimeGranularity) seconds. With the value set to 1/60 the
         * maximum lifetime or age of a particle is 2.27 years.
         */
        lifetimeGranularity: number;
        Copy(def: ParticleSystemDef): ParticleSystemDef;
        Clone(): ParticleSystemDef;
    }
    class ParticleSystem {
        paused: boolean;
        timestamp: number;
        allParticleFlags: ParticleFlag;
        needsUpdateAllParticleFlags: boolean;
        allGroupFlags: ParticleGroupFlag;
        needsUpdateAllGroupFlags: boolean;
        hasForce: boolean;
        iterationIndex: number;
        inverseDensity: number;
        particleDiameter: number;
        inverseDiameter: number;
        squaredDiameter: number;
        count: number;
        internalAllocatedCapacity: number;
        /**
         * Allocator for ParticleHandle instances.
         */
        /**
         * Maps particle indicies to handles.
         */
        handleIndexBuffer: ParticleSysteUserOverridableBuffer<ParticleHandle | null>;
        flagsBuffer: ParticleSysteUserOverridableBuffer<ParticleFlag>;
        positionBuffer: ParticleSysteUserOverridableBuffer<Vec2>;
        velocityBuffer: ParticleSysteUserOverridableBuffer<Vec2>;
        forceBuffer: Vec2[];
        /**
         * this.weightBuffer is populated in ComputeWeight and used in
         * ComputeDepth(), SolveStaticPressure() and SolvePressure().
         */
        weightBuffer: number[];
        /**
         * When any particles have the flag staticPressureParticle,
         * this.staticPressureBuffer is first allocated and used in
         * SolveStaticPressure() and SolvePressure().  It will be
         * reallocated on subsequent CreateParticle() calls.
         */
        staticPressureBuffer: number[];
        /**
         * this.accumulationBuffer is used in many functions as a temporary
         * buffer for scalar values.
         */
        accumulationBuffer: number[];
        /**
         * When any particles have the flag tensileParticle,
         * this.accumulation2Buffer is first allocated and used in
         * SolveTensile() as a temporary buffer for vector values.  It
         * will be reallocated on subsequent CreateParticle() calls.
         */
        accumulation2Buffer: Vec2[];
        /**
         * When any particle groups have the flag solidParticleGroup,
         * this.depthBuffer is first allocated and populated in
         * ComputeDepth() and used in SolveSolid(). It will be
         * reallocated on subsequent CreateParticle() calls.
         */
        depthBuffer: number[];
        colorBuffer: ParticleSysteUserOverridableBuffer<Color>;
        groupBuffer: Array<ParticleGroup | null>;
        userDataBuffer: ParticleSysteUserOverridableBuffer<any>;
        /**
         * Stuck particle detection parameters and record keeping
         */
        stuckThreshold: number;
        lastBodyContactStepBuffer: ParticleSysteUserOverridableBuffer<number>;
        bodyContactCountBuffer: ParticleSysteUserOverridableBuffer<number>;
        consecutiveContactStepsBuffer: ParticleSysteUserOverridableBuffer<number>;
        stuckParticleBuffer: GrowableBuffer<number>;
        proxyBuffer: GrowableBuffer<ParticleSysteProxy>;
        contactBuffer: GrowableBuffer<ParticleContact>;
        bodyContactBuffer: GrowableBuffer<ParticleBodyContact>;
        pairBuffer: GrowableBuffer<ParticlePair>;
        triadBuffer: GrowableBuffer<ParticleTriad>;
        /**
         * Time each particle should be destroyed relative to the last
         * time this.timeElapsed was initialized.  Each unit of time
         * corresponds to ParticleSystemDef::lifetimeGranularity
         * seconds.
         */
        expirationTimeBuffer: ParticleSysteUserOverridableBuffer<number>;
        /**
         * List of particle indices sorted by expiration time.
         */
        indexByExpirationTimeBuffer: ParticleSysteUserOverridableBuffer<number>;
        /**
         * Time elapsed in 32:32 fixed point.  Each non-fractional unit
         * of time corresponds to
         * ParticleSystemDef::lifetimeGranularity seconds.
         */
        timeElapsed: number;
        /**
         * Whether the expiration time buffer has been modified and
         * needs to be resorted.
         */
        expirationTimeBufferRequiresSorting: boolean;
        groupCount: number;
        groupList: ParticleGroup | null;
        def: ParticleSystemDef;
        world: World;
        prev: ParticleSystem | null;
        next: ParticleSystem | null;
        static readonly xTruncBits: number;
        static readonly yTruncBits: number;
        static readonly tagBits: number;
        static readonly yOffset: number;
        static readonly yShift: number;
        static readonly xShift: number;
        static readonly xScale: number;
        static readonly xOffset: number;
        static readonly yMask: number;
        static readonly xMask: number;
        static computeTag(x: number, y: number): number;
        static computeRelativeTag(tag: number, x: number, y: number): number;
        constructor(def: ParticleSystemDef, world: World);
        Drop(): void;
        /**
         * Create a particle whose properties have been defined.
         *
         * No reference to the definition is retained.
         *
         * A simulation step must occur before it's possible to interact
         * with a newly created particle.  For example,
         * DestroyParticleInShape() will not destroy a particle until
         * World::Step() has been called.
         *
         * warning: This function is locked during callbacks.
         */
        CreateParticle(def: IParticleDef): number;
        /**
         * Retrieve a handle to the particle at the specified index.
         *
         * Please see #ParticleHandle for why you might want a handle.
         */
        GetParticleHandleFromIndex(index: number): ParticleHandle;
        /**
         * Destroy a particle.
         *
         * The particle is removed after the next simulation step (see
         * World::Step()).
         *
         * @param index Index of the particle to destroy.
         * @param callDestructionListener Whether to call the
         *      destruction listener just before the particle is
         *      destroyed.
         */
        DestroyParticle(index: number, callDestructionListener?: boolean): void;
        /**
         * Destroy the Nth oldest particle in the system.
         *
         * The particle is removed after the next World::Step().
         *
         * @param index Index of the Nth oldest particle to
         *      destroy, 0 will destroy the oldest particle in the
         *      system, 1 will destroy the next oldest particle etc.
         * @param callDestructionListener Whether to call the
         *      destruction listener just before the particle is
         *      destroyed.
         */
        DestroyOldestParticle(index: number, callDestructionListener?: boolean): void;
        /**
         * Destroy particles inside a shape.
         *
         * warning: This function is locked during callbacks.
         *
         * In addition, this function immediately destroys particles in
         * the shape in constrast to DestroyParticle() which defers the
         * destruction until the next simulation step.
         *
         * @return Number of particles destroyed.
         * @param shape Shape which encloses particles
         *      that should be destroyed.
         * @param xf Transform applied to the shape.
         * @param callDestructionListener Whether to call the
         *      world DestructionListener for each particle
         *      destroyed.
         */
        DestroyParticlesInShape(shape: Shape, xf: Transform, callDestructionListener?: boolean): number;
        static readonly DestroyParticlesInShape_s_aabb: AABB;
        /**
         * Create a particle group whose properties have been defined.
         *
         * No reference to the definition is retained.
         *
         * warning: This function is locked during callbacks.
         */
        CreateParticleGroup(groupDef: IParticleGroupDef): ParticleGroup;
        static readonly CreateParticleGroup_s_transform: Transform;
        /**
         * Join two particle groups.
         *
         * warning: This function is locked during callbacks.
         *
         * @param groupA the first group. Expands to encompass the second group.
         * @param groupB the second group. It is destroyed.
         */
        JoinParticleGroups(groupA: ParticleGroup, groupB: ParticleGroup): void;
        /**
         * Split particle group into multiple disconnected groups.
         *
         * warning: This function is locked during callbacks.
         *
         * @param group the group to be split.
         */
        SplitParticleGroup(group: ParticleGroup): void;
        /**
         * Get the world particle group list. With the returned group,
         * use ParticleGroup::GetNext to get the next group in the
         * world list.
         *
         * A null group indicates the end of the list.
         *
         * @return the head of the world particle group list.
         */
        GetParticleGroupList(): ParticleGroup | null;
        /**
         * Get the number of particle groups.
         */
        GetParticleGroupCount(): number;
        /**
         * Get the number of particles.
         */
        GetParticleCount(): number;
        /**
         * Get the maximum number of particles.
         */
        GetMaxParticleCount(): number;
        /**
         * Set the maximum number of particles.
         *
         * A value of 0 means there is no maximum. The particle buffers
         * can continue to grow while World's block allocator still
         * has memory.
         *
         * Note: If you try to CreateParticle() with more than this
         * count, invalidParticleIndex is returned unless
         * SetDestructionByAge() is used to enable the destruction of
         * the oldest particles in the system.
         */
        SetMaxParticleCount(count: number): void;
        /**
         * Get all existing particle flags.
         */
        GetAllParticleFlags(): ParticleFlag;
        /**
         * Get all existing particle group flags.
         */
        GetAllGroupFlags(): ParticleGroupFlag;
        /**
         * Pause or unpause the particle system. When paused,
         * World::Step() skips over this particle system. All
         * ParticleSystem function calls still work.
         *
         * @param paused paused is true to pause, false to un-pause.
         */
        SetPaused(paused: boolean): void;
        /**
         * Initially, true, then, the last value passed into
         * SetPaused().
         *
         * @return true if the particle system is being updated in World::Step().
         */
        GetPaused(): boolean;
        /**
         * Change the particle density.
         *
         * Particle density affects the mass of the particles, which in
         * turn affects how the particles interact with Bodies. Note
         * that the density does not affect how the particles interact
         * with each other.
         */
        SetDensity(density: number): void;
        /**
         * Get the particle density.
         */
        GetDensity(): number;
        /**
         * Change the particle gravity scale. Adjusts the effect of the
         * global gravity vector on particles.
         */
        SetGravityScale(gravityScale: number): void;
        /**
         * Get the particle gravity scale.
         */
        GetGravityScale(): number;
        /**
         * Damping is used to reduce the velocity of particles. The
         * damping parameter can be larger than 1.0f but the damping
         * effect becomes sensitive to the time step when the damping
         * parameter is large.
         */
        SetDamping(damping: number): void;
        /**
         * Get damping for particles
         */
        GetDamping(): number;
        /**
         * Change the number of iterations when calculating the static
         * pressure of particles. By default, 8 iterations. You can
         * reduce the number of iterations down to 1 in some situations,
         * but this may cause instabilities when many particles come
         * together. If you see particles popping away from each other
         * like popcorn, you may have to increase the number of
         * iterations.
         *
         * For a description of static pressure, see
         * http://en.wikipedia.org/wiki/Static_pressure#Static_pressure_in_fluid_dynamics
         */
        SetStaticPressureIterations(iterations: number): void;
        /**
         * Get the number of iterations for static pressure of
         * particles.
         */
        GetStaticPressureIterations(): number;
        /**
         * Change the particle radius.
         *
         * You should set this only once, on world start.
         * If you change the radius during execution, existing particles
         * may explode, shrink, or behave unexpectedly.
         */
        SetRadius(radius: number): void;
        /**
         * Get the particle radius.
         */
        GetRadius(): number;
        /**
         * Get the position of each particle
         *
         * Array is length GetParticleCount()
         *
         * @return the pointer to the head of the particle positions array.
         */
        GetPositionBuffer(): Vec2[];
        /**
         * Get the velocity of each particle
         *
         * Array is length GetParticleCount()
         *
         * @return the pointer to the head of the particle velocities array.
         */
        GetVelocityBuffer(): Vec2[];
        /**
         * Get the color of each particle
         *
         * Array is length GetParticleCount()
         *
         * @return the pointer to the head of the particle colors array.
         */
        GetColorBuffer(): Color[];
        /**
         * Get the particle-group of each particle.
         *
         * Array is length GetParticleCount()
         *
         * @return the pointer to the head of the particle group array.
         */
        GetGroupBuffer(): Array<ParticleGroup | null>;
        /**
         * Get the weight of each particle
         *
         * Array is length GetParticleCount()
         *
         * @return the pointer to the head of the particle positions array.
         */
        GetWeightBuffer(): number[];
        /**
         * Get the user-specified data of each particle.
         *
         * Array is length GetParticleCount()
         *
         * @return the pointer to the head of the particle user-data array.
         */
        GetUserDataBuffer<T>(): T[];
        /**
         * Get the flags for each particle. See the ParticleFlag enum.
         *
         * Array is length GetParticleCount()
         *
         * @return the pointer to the head of the particle-flags array.
         */
        GetFlagsBuffer(): ParticleFlag[];
        /**
         * Set flags for a particle. See the ParticleFlag enum.
         */
        SetParticleFlags(index: number, newFlags: ParticleFlag): void;
        /**
         * Get flags for a particle. See the ParticleFlag enum.
         */
        GetParticleFlags(index: number): ParticleFlag;
        /**
         * Set an external buffer for particle data.
         *
         * Normally, the World's block allocator is used for particle
         * data. However, sometimes you may have an OpenGL or Java
         * buffer for particle data. To avoid data duplication, you may
         * supply this external buffer.
         *
         * Note that, when World's block allocator is used, the
         * particle data buffers can grow as required. However, when
         * external buffers are used, the maximum number of particles is
         * clamped to the size of the smallest external buffer.
         *
         * @param buffer a pointer to a block of memory.
         * @param capacity the number of values in the block.
         */
        SetFlagsBuffer(buffer: ParticleFlag[]): void;
        SetPositionBuffer(buffer: Vec2[] | Float32Array): void;
        SetVelocityBuffer(buffer: TypedVec2[] | Float32Array): void;
        SetColorBuffer(buffer: Color[] | Float32Array): void;
        SetUserDataBuffer<T>(buffer: T[]): void;
        /**
         * Get contacts between particles
         * Contact data can be used for many reasons, for example to
         * trigger rendering or audio effects.
         */
        GetContacts(): ParticleContact[];
        GetContactCount(): number;
        /**
         * Get contacts between particles and bodies
         *
         * Contact data can be used for many reasons, for example to
         * trigger rendering or audio effects.
         */
        GetBodyContacts(): ParticleBodyContact[];
        GetBodyContactCount(): number;
        /**
         * Get array of particle pairs. The particles in a pair:
         *   (1) are contacting,
         *   (2) are in the same particle group,
         *   (3) are part of a rigid particle group, or are spring, elastic,
         *       or wall particles.
         *   (4) have at least one particle that is a spring or barrier
         *       particle (i.e. one of the types in k_pairFlags),
         *   (5) have at least one particle that returns true for
         *       ConnectionFilter::IsNecessary,
         *   (6) are not zombie particles.
         *
         * Essentially, this is an array of spring or barrier particles
         * that are interacting. The array is sorted by ParticlePair's
         * indexA, and then indexB. There are no duplicate entries.
         */
        GetPairs(): ParticlePair[];
        GetPairCount(): number;
        /**
         * Get array of particle triads. The particles in a triad:
         *   (1) are in the same particle group,
         *   (2) are in a Voronoi triangle together,
         *   (3) are within maxTriadDistance particle diameters of each
         *       other,
         *   (4) return true for ConnectionFilter::ShouldCreateTriad
         *   (5) have at least one particle of type elastic (i.e. one of the
         *       types in k_triadFlags),
         *   (6) are part of a rigid particle group, or are spring, elastic,
         *       or wall particles.
         *   (7) are not zombie particles.
         *
         * Essentially, this is an array of elastic particles that are
         * interacting. The array is sorted by ParticleTriad's indexA,
         * then indexB, then indexC. There are no duplicate entries.
         */
        GetTriads(): ParticleTriad[];
        GetTriadCount(): number;
        /**
         * Set an optional threshold for the maximum number of
         * consecutive particle iterations that a particle may contact
         * multiple bodies before it is considered a candidate for being
         * "stuck". Setting to zero or less disables.
         */
        SetStuckThreshold(steps: number): void;
        /**
         * Get potentially stuck particles from the last step; the user
         * must decide if they are stuck or not, and if so, delete or
         * move them
         */
        GetStuckCandidates(): number[];
        /**
         * Get the number of stuck particle candidates from the last
         * step.
         */
        GetStuckCandidateCount(): number;
        /**
         * Compute the kinetic energy that can be lost by damping force
         */
        ComputeCollisionEnergy(): number;
        static readonly ComputeCollisionEnergy_s_v: Vec2;
        /**
         * Set strict Particle/Body contact check.
         *
         * This is an option that will help ensure correct behavior if
         * there are corners in the world model where Particle/Body
         * contact is ambiguous. This option scales at n*log(n) of the
         * number of Particle/Body contacts, so it is best to only
         * enable if it is necessary for your geometry. Enable if you
         * see strange particle behavior around Body intersections.
         */
        SetStrictContactCheck(enabled: boolean): void;
        /**
         * Get the status of the strict contact check.
         */
        GetStrictContactCheck(): boolean;
        /**
         * Set the lifetime (in seconds) of a particle relative to the
         * current time.  A lifetime of less than or equal to 0.0f
         * results in the particle living forever until it's manually
         * destroyed by the application.
         */
        SetParticleLifetime(index: number, lifetime: number): void;
        /**
         * Get the lifetime (in seconds) of a particle relative to the
         * current time.  A value > 0.0f is returned if the particle is
         * scheduled to be destroyed in the future, values <= 0.0f
         * indicate the particle has an infinite lifetime.
         */
        GetParticleLifetime(index: number): number;
        /**
         * Enable / disable destruction of particles in CreateParticle()
         * when no more particles can be created due to a prior call to
         * SetMaxParticleCount().  When this is enabled, the oldest
         * particle is destroyed in CreateParticle() favoring the
         * destruction of particles with a finite lifetime over
         * particles with infinite lifetimes. This feature is enabled by
         * default when particle lifetimes are tracked.  Explicitly
         * enabling this feature using this function enables particle
         * lifetime tracking.
         */
        SetDestructionByAge(enable: boolean): void;
        /**
         * Get whether the oldest particle will be destroyed in
         * CreateParticle() when the maximum number of particles are
         * present in the system.
         */
        GetDestructionByAge(): boolean;
        /**
         * Get the array of particle expiration times indexed by
         * particle index.
         *
         * GetParticleCount() items are in the returned array.
         */
        GetExpirationTimeBuffer(): number[];
        /**
         * Convert a expiration time value in returned by
         * GetExpirationTimeBuffer() to a time in seconds relative to
         * the current simulation time.
         */
        ExpirationTimeToLifetime(expirationTime: number): number;
        /**
         * Get the array of particle indices ordered by reverse
         * lifetime. The oldest particle indexes are at the end of the
         * array with the newest at the start.  Particles with infinite
         * lifetimes (i.e expiration times less than or equal to 0) are
         * placed at the start of the array.
         * ExpirationTimeToLifetime(GetExpirationTimeBuffer()[index]) is
         * equivalent to GetParticleLifetime(index).
         *
         * GetParticleCount() items are in the returned array.
         */
        GetIndexByExpirationTimeBuffer(): number[];
        /**
         * Apply an impulse to one particle. This immediately modifies
         * the velocity. Similar to Body::ApplyLinearImpulse.
         *
         * @param index the particle that will be modified.
         * @param impulse impulse the world impulse vector, usually in N-seconds or kg-m/s.
         */
        ParticleApplyLinearImpulse(index: number, impulse: XY): void;
        /**
         * Apply an impulse to all particles between 'firstIndex' and
         * 'lastIndex'. This immediately modifies the velocity. Note
         * that the impulse is applied to the total mass of all
         * particles. So, calling ParticleApplyLinearImpulse(0, impulse)
         * and ParticleApplyLinearImpulse(1, impulse) will impart twice
         * as much velocity as calling just ApplyLinearImpulse(0, 1,
         * impulse).
         *
         * @param firstIndex the first particle to be modified.
         * @param lastIndex the last particle to be modified.
         * @param impulse the world impulse vector, usually in N-seconds or kg-m/s.
         */
        ApplyLinearImpulse(firstIndex: number, lastIndex: number, impulse: XY): void;
        static IsSignificantForce(force: XY): boolean;
        /**
         * Apply a force to the center of a particle.
         *
         * @param index the particle that will be modified.
         * @param force the world force vector, usually in Newtons (N).
         */
        ParticleApplyForce(index: number, force: XY): void;
        /**
         * Distribute a force across several particles. The particles
         * must not be wall particles. Note that the force is
         * distributed across all the particles, so calling this
         * function for indices 0..N is not the same as calling
         * ParticleApplyForce(i, force) for i in 0..N.
         *
         * @param firstIndex the first particle to be modified.
         * @param lastIndex the last particle to be modified.
         * @param force the world force vector, usually in Newtons (N).
         */
        ApplyForce(firstIndex: number, lastIndex: number, force: XY): void;
        /**
         * Get the next particle-system in the world's particle-system
         * list.
         */
        GetNext(): ParticleSystem | null;
        /**
         * Query the particle system for all particles that potentially
         * overlap the provided AABB.
         * QueryCallback::ShouldQueryParticleSystem is ignored.
         *
         * @param callback a user implemented callback class.
         * @param aabb the query box.
         */
        QueryAABB(callback: QueryCallback, aabb: AABB): void;
        /**
         * Query the particle system for all particles that potentially
         * overlap the provided shape's AABB. Calls QueryAABB
         * internally. QueryCallback::ShouldQueryParticleSystem is
         * ignored.
         *
         * @param callback a user implemented callback class.
         * @param shape the query shape
         * @param xf the transform of the AABB
         * @param childIndex
         */
        QueryShapeAABB(callback: QueryCallback, shape: Shape, xf: Transform, childIndex?: number): void;
        static readonly QueryShapeAABB_s_aabb: AABB;
        QueryPointAABB(callback: QueryCallback, point: XY, slop?: number): void;
        static readonly QueryPointAABB_s_aabb: AABB;
        /**
         * Ray-cast the particle system for all particles in the path of
         * the ray. Your callback controls whether you get the closest
         * point, any point, or n-points. The ray-cast ignores particles
         * that contain the starting point.
         * RayCastCallback::ShouldQueryParticleSystem is ignored.
         *
         * @param callback a user implemented callback class.
         * @param point1 the ray starting point
         * @param point2 the ray ending point
         */
        RayCast(callback: RayCastCallback, point1: XY, point2: XY): void;
        static readonly RayCast_s_aabb: AABB;
        static readonly RayCast_s_p: Vec2;
        static readonly RayCast_s_v: Vec2;
        static readonly RayCast_s_n: Vec2;
        static readonly RayCast_s_point: Vec2;
        /**
         * Compute the axis-aligned bounding box for all particles
         * contained within this particle system.
         * @param aabb Returns the axis-aligned bounding box of the system.
         */
        ComputeAABB(aabb: AABB): void;
        /**
         * All particle types that require creating pairs
         */
        static readonly k_pairFlags: number;
        /**
         * All particle types that require creating triads
         */
        static readonly k_triadFlags: ParticleFlag;
        /**
         * All particle types that do not produce dynamic pressure
         */
        static readonly k_noPressureFlags: number;
        /**
         * All particle types that apply extra damping force with bodies
         */
        static readonly k_extraDampingFlags: ParticleFlag;
        static readonly k_barrierWallFlags: number;
        FreeBuffer<T>(b: T[] | null, capacity: number): void;
        FreeUserOverridableBuffer<T>(b: ParticleSysteUserOverridableBuffer<T>): void;
        /**
         * Reallocate a buffer
         */
        ReallocateBuffer3<T>(oldBuffer: T[] | null, oldCapacity: number, newCapacity: number): T[];
        /**
         * Reallocate a buffer
         */
        ReallocateBuffer5<T>(buffer: T[] | null, userSuppliedCapacity: number, oldCapacity: number, newCapacity: number, deferred: boolean): T[];
        /**
         * Reallocate a buffer
         */
        ReallocateBuffer4<T>(buffer: ParticleSysteUserOverridableBuffer<any>, oldCapacity: number, newCapacity: number, deferred: boolean): T[];
        RequestBuffer<T>(buffer: T[] | null): T[];
        /**
         * Reallocate the handle / index map and schedule the allocation
         * of a new pool for handle allocation.
         */
        ReallocateHandleBuffers(newCapacity: number): void;
        ReallocateInternalAllocatedBuffers(capacity: number): void;
        CreateParticleForGroup(groupDef: IParticleGroupDef, xf: Transform, p: XY): void;
        CreateParticlesStrokeShapeForGroup(shape: Shape, groupDef: IParticleGroupDef, xf: Transform): void;
        static readonly CreateParticlesStrokeShapeForGroup_s_edge: EdgeShape;
        static readonly CreateParticlesStrokeShapeForGroup_s_d: Vec2;
        static readonly CreateParticlesStrokeShapeForGroup_s_p: Vec2;
        CreateParticlesFillShapeForGroup(shape: Shape, groupDef: IParticleGroupDef, xf: Transform): void;
        static readonly CreateParticlesFillShapeForGroup_s_aabb: AABB;
        static readonly CreateParticlesFillShapeForGroup_s_p: Vec2;
        CreateParticlesWithShapeForGroup(shape: Shape, groupDef: IParticleGroupDef, xf: Transform): void;
        CreateParticlesWithShapesForGroup(shapes: Shape[], shapeCount: number, groupDef: IParticleGroupDef, xf: Transform): void;
        CloneParticle(oldIndex: number, group: ParticleGroup): number;
        DestroyParticlesInGroup(group: ParticleGroup, callDestructionListener?: boolean): void;
        DestroyParticleGroup(group: ParticleGroup): void;
        static ParticleCanBeConnected(flags: ParticleFlag, group: ParticleGroup | null): boolean;
        UpdatePairsAndTriads(firstIndex: number, lastIndex: number, filter: ParticleSysteConnectionFilter): void;
        private static UpdatePairsAndTriads_s_dab;
        private static UpdatePairsAndTriads_s_dbc;
        private static UpdatePairsAndTriads_s_dca;
        UpdatePairsAndTriadsWithReactiveParticles(): void;
        static ComparePairIndices(a: ParticlePair, b: ParticlePair): boolean;
        static MatchPairIndices(a: ParticlePair, b: ParticlePair): boolean;
        static CompareTriadIndices(a: ParticleTriad, b: ParticleTriad): boolean;
        static MatchTriadIndices(a: ParticleTriad, b: ParticleTriad): boolean;
        static InitializeParticleLists(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[]): void;
        MergeParticleListsInContact(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[]): void;
        static MergeParticleLists(listA: ParticleSysteParticleListNode, listB: ParticleSysteParticleListNode): void;
        static FindLongestParticleList(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[]): ParticleSysteParticleListNode;
        MergeZombieParticleListNodes(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[], survivingList: ParticleSysteParticleListNode): void;
        static MergeParticleListAndNode(list: ParticleSysteParticleListNode, node: ParticleSysteParticleListNode): void;
        CreateParticleGroupsFromParticleList(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[], survivingList: ParticleSysteParticleListNode): void;
        UpdatePairsAndTriadsWithParticleList(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[]): void;
        ComputeDepth(): void;
        GetInsideBoundsEnumerator(aabb: AABB): ParticleSysteInsideBoundsEnumerator;
        UpdateAllParticleFlags(): void;
        UpdateAllGroupFlags(): void;
        AddContact(a: number, b: number, contacts: GrowableBuffer<ParticleContact>): void;
        static readonly AddContact_s_d: Vec2;
        FindContacts_Reference(contacts: GrowableBuffer<ParticleContact>): void;
        FindContacts(contacts: GrowableBuffer<ParticleContact>): void;
        UpdateProxies_Reference(proxies: GrowableBuffer<ParticleSysteProxy>): void;
        UpdateProxies(proxies: GrowableBuffer<ParticleSysteProxy>): void;
        SortProxies(proxies: GrowableBuffer<ParticleSysteProxy>): void;
        FilterContacts(contacts: GrowableBuffer<ParticleContact>): void;
        NotifyContactListenerPreContact(particlePairs: ParticlePairSet): void;
        NotifyContactListenerPostContact(particlePairs: ParticlePairSet): void;
        static ParticleContactIsZombie(contact: ParticleContact): boolean;
        UpdateContacts(exceptZombie: boolean): void;
        NotifyBodyContactListenerPreContact(fixtureSet: ParticleSysteFixtureParticleSet): void;
        NotifyBodyContactListenerPostContact(fixtureSet: ParticleSysteFixtureParticleSet): void;
        UpdateBodyContacts(): void;
        static readonly UpdateBodyContacts_s_aabb: AABB;
        UpdateBodyContacts_callback: ParticleSysteUpdateBodyContactsCallback | null;
        Solve(step: TimeStep): void;
        static readonly Solve_s_subStep: TimeStep;
        SolveCollision(step: TimeStep): void;
        static readonly SolveCollision_s_aabb: AABB;
        SolveCollision_callback: ParticleSysteSolveCollisionCallback | null;
        LimitVelocity(step: TimeStep): void;
        SolveGravity(step: TimeStep): void;
        static readonly SolveGravity_s_gravity: Vec2;
        SolveBarrier(step: TimeStep): void;
        static readonly SolveBarrier_s_aabb: AABB;
        static readonly SolveBarrier_s_va: Vec2;
        static readonly SolveBarrier_s_vb: Vec2;
        static readonly SolveBarrier_s_pba: Vec2;
        static readonly SolveBarrier_s_vba: Vec2;
        static readonly SolveBarrier_s_vc: Vec2;
        static readonly SolveBarrier_s_pca: Vec2;
        static readonly SolveBarrier_s_vca: Vec2;
        static readonly SolveBarrier_s_qba: Vec2;
        static readonly SolveBarrier_s_qca: Vec2;
        static readonly SolveBarrier_s_dv: Vec2;
        static readonly SolveBarrier_s_f: Vec2;
        SolveStaticPressure(step: TimeStep): void;
        ComputeWeight(): void;
        SolvePressure(step: TimeStep): void;
        static readonly SolvePressure_s_f: Vec2;
        SolveDamping(step: TimeStep): void;
        static readonly SolveDamping_s_v: Vec2;
        static readonly SolveDamping_s_f: Vec2;
        SolveRigidDamping(): void;
        static readonly SolveRigidDamping_s_t0: Vec2;
        static readonly SolveRigidDamping_s_t1: Vec2;
        static readonly SolveRigidDamping_s_p: Vec2;
        static readonly SolveRigidDamping_s_v: Vec2;
        SolveExtraDamping(): void;
        static readonly SolveExtraDamping_s_v: Vec2;
        static readonly SolveExtraDamping_s_f: Vec2;
        SolveWall(): void;
        SolveRigid(step: TimeStep): void;
        static readonly SolveRigid_s_position: Vec2;
        static readonly SolveRigid_s_rotation: Rot;
        static readonly SolveRigid_s_transform: Transform;
        static readonly SolveRigid_s_velocityTransform: Transform;
        SolveElastic(step: TimeStep): void;
        static readonly SolveElastic_s_pa: Vec2;
        static readonly SolveElastic_s_pb: Vec2;
        static readonly SolveElastic_s_pc: Vec2;
        static readonly SolveElastic_s_r: Rot;
        static readonly SolveElastic_s_t0: Vec2;
        SolveSpring(step: TimeStep): void;
        static readonly SolveSpring_s_pa: Vec2;
        static readonly SolveSpring_s_pb: Vec2;
        static readonly SolveSpring_s_d: Vec2;
        static readonly SolveSpring_s_f: Vec2;
        SolveTensile(step: TimeStep): void;
        static readonly SolveTensile_s_weightedNormal: Vec2;
        static readonly SolveTensile_s_s: Vec2;
        static readonly SolveTensile_s_f: Vec2;
        SolveViscous(): void;
        static readonly SolveViscous_s_v: Vec2;
        static readonly SolveViscous_s_f: Vec2;
        SolveRepulsive(step: TimeStep): void;
        static readonly SolveRepulsive_s_f: Vec2;
        SolvePowder(step: TimeStep): void;
        static readonly SolvePowder_s_f: Vec2;
        SolveSolid(step: TimeStep): void;
        static readonly SolveSolid_s_f: Vec2;
        SolveForce(step: TimeStep): void;
        SolveColorMixing(): void;
        SolveZombie(): void;
        /**
         * Destroy all particles which have outlived their lifetimes set
         * by SetParticleLifetime().
         */
        SolveLifetimes(step: TimeStep): void;
        RotateBuffer(start: number, mid: number, end: number): void;
        GetCriticalVelocity(step: TimeStep): number;
        GetCriticalVelocitySquared(step: TimeStep): number;
        GetCriticalPressure(step: TimeStep): number;
        GetParticleStride(): number;
        GetParticleMass(): number;
        GetParticleInvMass(): number;
        /**
         * Get the world's contact filter if any particles with the
         * contactFilterParticle flag are present in the system.
         */
        GetFixtureContactFilter(): ContactFilter | null;
        /**
         * Get the world's contact filter if any particles with the
         * particleContactFilterParticle flag are present in the
         * system.
         */
        GetParticleContactFilter(): ContactFilter | null;
        /**
         * Get the world's contact listener if any particles with the
         * fixtureContactListenerParticle flag are present in the
         * system.
         */
        GetFixtureContactListener(): ContactListener | null;
        /**
         * Get the world's contact listener if any particles with the
         * particleContactListenerParticle flag are present in the
         * system.
         */
        GetParticleContactListener(): ContactListener | null;
        SetUserOverridableBuffer<T>(buffer: ParticleSysteUserOverridableBuffer<T>, data: T[]): void;
        SetGroupFlags(group: ParticleGroup, newFlags: ParticleGroupFlag): void;
        static BodyContactCompare(lhs: ParticleBodyContact, rhs: ParticleBodyContact): boolean;
        RemoveSpuriousBodyContacts(): void;
        private static RemoveSpuriousBodyContacts_s_n;
        private static RemoveSpuriousBodyContacts_s_pos;
        private static RemoveSpuriousBodyContacts_s_normal;
        DetectStuckParticle(particle: number): void;
        /**
         * Determine whether a particle index is valid.
         */
        ValidateParticleIndex(index: number): boolean;
        /**
         * Get the time elapsed in
         * ParticleSystemDef::lifetimeGranularity.
         */
        GetQuantizedTimeElapsed(): number;
        /**
         * Convert a lifetime in seconds to an expiration time.
         */
        LifetimeToExpirationTime(lifetime: number): number;
        ForceCanBeApplied(flags: ParticleFlag): boolean;
        PrepareForceBuffer(): void;
        IsRigidGroup(group: ParticleGroup | null): boolean;
        GetLinearVelocity(group: ParticleGroup | null, particleIndex: number, point: Vec2, out: Vec2): Vec2;
        InitDampingParameter(invMass: number[], invInertia: number[], tangentDistance: number[], mass: number, inertia: number, center: Vec2, point: Vec2, normal: Vec2): void;
        InitDampingParameterWithRigidGroupOrParticle(invMass: number[], invInertia: number[], tangentDistance: number[], isRigidGroup: boolean, group: ParticleGroup | null, particleIndex: number, point: Vec2, normal: Vec2): void;
        ComputeDampingImpulse(invMassA: number, invInertiaA: number, tangentDistanceA: number, invMassB: number, invInertiaB: number, tangentDistanceB: number, normalVelocity: number): number;
        ApplyDamping(invMass: number, invInertia: number, tangentDistance: number, isRigidGroup: boolean, group: ParticleGroup | null, particleIndex: number, impulse: number, normal: Vec2): void;
    }
    class ParticleSysteUserOverridableBuffer<T> {
        _data: T[] | null;
        data: T[];
        userSuppliedCapacity: number;
    }
    class ParticleSysteProxy {
        index: number;
        tag: number;
        static CompareProxyProxy(a: ParticleSysteProxy, b: ParticleSysteProxy): boolean;
        static CompareTagProxy(a: number, b: ParticleSysteProxy): boolean;
        static CompareProxyTag(a: ParticleSysteProxy, b: number): boolean;
    }
    class ParticleSysteInsideBoundsEnumerator {
        system: ParticleSystem;
        xLower: number;
        xUpper: number;
        yLower: number;
        yUpper: number;
        first: number;
        last: number;
        /**
         * InsideBoundsEnumerator enumerates all particles inside the
         * given bounds.
         *
         * Construct an enumerator with bounds of tags and a range of
         * proxies.
         */
        constructor(system: ParticleSystem, lower: number, upper: number, first: number, last: number);
        /**
         * Get index of the next particle. Returns
         * invalidParticleIndex if there are no more particles.
         */
        GetNext(): number;
    }
    class ParticleSysteParticleListNode {
        /**
         * The head of the list.
         */
        list: ParticleSysteParticleListNode;
        /**
         * The next node in the list.
         */
        next: ParticleSysteParticleListNode | null;
        /**
         * Number of entries in the list. Valid only for the node at the
         * head of the list.
         */
        count: number;
        /**
         * Particle index.
         */
        index: number;
    }
    /**
     * @constructor
     */
    class ParticleSysteFixedSetAllocator<T> {
        Allocate(itemSize: number, count: number): number;
        Clear(): void;
        GetCount(): number;
        Invalidate(itemIndex: number): void;
        GetValidBuffer(): boolean[];
        GetBuffer(): T[];
        SetCount(count: number): void;
    }
    class ParticleSysteFixtureParticle {
        first: Fixture;
        second: number;
        constructor(fixture: Fixture, particle: number);
    }
    class ParticleSysteFixtureParticleSet extends ParticleSysteFixedSetAllocator<ParticleSysteFixtureParticle> {
        Initialize(bodyContactBuffer: GrowableBuffer<ParticleBodyContact>, flagsBuffer: ParticleSysteUserOverridableBuffer<ParticleFlag>): void;
        Find(pair: ParticleSysteFixtureParticle): number;
    }
    class ParticleSysteParticlePair {
        first: number;
        second: number;
        constructor(particleA: number, particleB: number);
    }
    class ParticlePairSet extends ParticleSysteFixedSetAllocator<ParticleSysteParticlePair> {
        Initialize(contactBuffer: GrowableBuffer<ParticleContact>, flagsBuffer: ParticleSysteUserOverridableBuffer<ParticleFlag>): void;
        Find(pair: ParticleSysteParticlePair): number;
    }
    class ParticleSysteConnectionFilter {
        /**
         * Is the particle necessary for connection?
         * A pair or a triad should contain at least one 'necessary'
         * particle.
         */
        IsNecessary(index: number): boolean;
        /**
         * An additional condition for creating a pair.
         */
        ShouldCreatePair(a: number, b: number): boolean;
        /**
         * An additional condition for creating a triad.
         */
        ShouldCreateTriad(a: number, b: number, c: number): boolean;
    }
    class ParticleSysteDestroyParticlesInShapeCallback extends QueryCallback {
        system: ParticleSystem;
        shape: Shape;
        xf: Transform;
        callDestructionListener: boolean;
        destroyed: number;
        constructor(system: ParticleSystem, shape: Shape, xf: Transform, callDestructionListener: boolean);
        ReportFixture(fixture: Fixture): boolean;
        ReportParticle(particleSystem: ParticleSystem, index: number): boolean;
        Destroyed(): number;
    }
    class ParticleSysteJoinParticleGroupsFilter extends ParticleSysteConnectionFilter {
        threshold: number;
        constructor(threshold: number);
        /**
         * An additional condition for creating a pair.
         */
        ShouldCreatePair(a: number, b: number): boolean;
        /**
         * An additional condition for creating a triad.
         */
        ShouldCreateTriad(a: number, b: number, c: number): boolean;
    }
    class ParticleSysteCompositeShape extends Shape {
        constructor(shapes: Shape[], shapeCount?: number);
        shapes: Shape[];
        shapeCount: number;
        Clone(): Shape;
        GetChildCount(): number;
        /**
         * @see Shape::TestPoint
         */
        TestPoint(xf: Transform, p: XY): boolean;
        /**
         * @see Shape::ComputeDistance
         */
        ComputeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number;
        /**
         * Implement Shape.
         */
        RayCast(output: RayCastOutput, input: RayCastInput, xf: Transform, childIndex: number): boolean;
        /**
         * @see Shape::ComputeAABB
         */
        ComputeAABB(aabb: AABB, xf: Transform, childIndex: number): void;
        /**
         * @see Shape::ComputeMass
         */
        ComputeMass(massData: MassData, density: number): void;
        SetupDistanceProxy(proxy: DistanceProxy, index: number): void;
        ComputeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number;
        Dump(log: (format: string, ...args: any[]) => void): void;
    }
    class ParticleSysteReactiveFilter extends ParticleSysteConnectionFilter {
        flagsBuffer: ParticleSysteUserOverridableBuffer<ParticleFlag>;
        constructor(flagsBuffer: ParticleSysteUserOverridableBuffer<ParticleFlag>);
        IsNecessary(index: number): boolean;
    }
    class ParticleSysteUpdateBodyContactsCallback extends FixtureParticleQueryCallback {
        contactFilter: ContactFilter | null;
        constructor(system: ParticleSystem, contactFilter?: ContactFilter | null);
        ShouldCollideFixtureParticle(fixture: Fixture, particleSystem: ParticleSystem, particleIndex: number): boolean;
        ReportFixtureAndParticle(fixture: Fixture, childIndex: number, a: number): void;
        static readonly ReportFixtureAndParticle_s_n: Vec2;
        static readonly ReportFixtureAndParticle_s_rp: Vec2;
    }
    class ParticleSysteSolveCollisionCallback extends FixtureParticleQueryCallback {
        step: TimeStep;
        constructor(system: ParticleSystem, step: TimeStep);
        ReportFixtureAndParticle(fixture: Fixture, childIndex: number, a: number): void;
        static readonly ReportFixtureAndParticle_s_p1: Vec2;
        static readonly ReportFixtureAndParticle_s_output: RayCastOutput;
        static readonly ReportFixtureAndParticle_s_input: RayCastInput;
        static readonly ReportFixtureAndParticle_s_p: Vec2;
        static readonly ReportFixtureAndParticle_s_v: Vec2;
        static readonly ReportFixtureAndParticle_s_f: Vec2;
        ReportParticle(system: ParticleSystem, index: number): boolean;
    }
}
declare namespace b2 {
    class StackQueue<T> {
        readonly buffer: Array<T | null>;
        front: number;
        back: number;
        readonly capacity: number;
        constructor(capacity: number);
        Push(item: T): void;
        Pop(): void;
        Empty(): boolean;
        Front(): T;
    }
}
/**
 * A field representing the nearest generator from each point.
 */
declare namespace b2 {
    class VoronoiDiagram {
        generatorBuffer: VoronoiDiagraGenerator[];
        generatorCapacity: number;
        generatorCount: number;
        countX: number;
        countY: number;
        diagram: VoronoiDiagraGenerator[];
        constructor(generatorCapacity: number);
        /**
         * Add a generator.
         *
         * @param center the position of the generator.
         * @param tag a tag used to identify the generator in callback functions.
         * @param necessary whether to callback for nodes associated with the generator.
         */
        AddGenerator(center: Vec2, tag: number, necessary: boolean): void;
        /**
         * Generate the Voronoi diagram. It is rasterized with a given
         * interval in the same range as the necessary generators exist.
         *
         * @param radius the interval of the diagram.
         * @param margin margin for which the range of the diagram is extended.
         */
        Generate(radius: number, margin: number): void;
        /**
         * Enumerate all nodes that contain at least one necessary
         * generator.
         */
        GetNodes(callback: VoronoiDiagraNodeCallback): void;
    }
    /**
     * Callback used by GetNodes().
     *
     * Receive tags for generators associated with a node.
     */
    type VoronoiDiagraNodeCallback = (a: number, b: number, c: number) => void;
    class VoronoiDiagraGenerator {
        center: Vec2;
        tag: number;
        necessary: boolean;
    }
    class VoronoiDiagraTask {
        x: number;
        y: number;
        i: number;
        generator: VoronoiDiagraGenerator;
        constructor(x: number, y: number, i: number, g: VoronoiDiagraGenerator);
    }
}
declare namespace b2 {
    enum StretchingModel {
        pbdStretchingModel = 0,
        xpbdStretchingModel = 1
    }
    enum BendingModel {
        springAngleBendingModel = 0,
        pbdAngleBendingModel = 1,
        xpbdAngleBendingModel = 2,
        pbdDistanceBendingModel = 3,
        pbdHeightBendingModel = 4,
        pbdTriangleBendingModel = 5
    }
    class RopeTuning {
        stretchingModel: StretchingModel;
        bendingModel: BendingModel;
        damping: number;
        stretchStiffness: number;
        stretchHertz: number;
        stretchDamping: number;
        bendStiffness: number;
        bendHertz: number;
        bendDamping: number;
        isometric: boolean;
        fixedEffectiveMass: boolean;
        warmStart: boolean;
        Copy(other: RopeTuning): this;
    }
    class RopeDef {
        readonly position: Vec2;
        readonly vertices: Vec2[];
        count: number;
        readonly masses: number[];
        readonly gravity: Vec2;
        readonly tuning: RopeTuning;
    }
    class Rope {
        private readonly position;
        private count;
        private stretchCount;
        private bendCount;
        private readonly stretchConstraints;
        private readonly bendConstraints;
        private readonly bindPositions;
        private readonly ps;
        private readonly p0s;
        private readonly vs;
        private readonly invMasses;
        private readonly gravity;
        private readonly tuning;
        Create(def: RopeDef): void;
        SetTuning(tuning: RopeTuning): void;
        Step(dt: number, iterations: number, position: Vec2): void;
        Reset(position: Vec2): void;
        Draw(draw: Draw): void;
        private SolveStretch_PBD;
        private SolveStretch_XPBD;
        private SolveBend_PBD_Angle;
        private SolveBend_XPBD_Angle;
        private SolveBend_PBD_Distance;
        private SolveBend_PBD_Height;
        private SolveBend_PBD_Triangle;
        private ApplyBendForces;
    }
}
