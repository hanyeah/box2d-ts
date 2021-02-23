declare namespace b2 {
    function assert(condition: boolean, ...args: any[]): void;
    function maybe<T>(value: T, def: T): T;
    const maxFloat: number;
    const epsilon: number;
    const epsilonSq: number;
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
    function MakeNullArray<T>(length: number): Array<T>;
    function MakeNumberArray(length: number, init?: number): number[];
}
declare namespace b2 {
    const pi_over_180: number;
    const _180_over_pi: number;
    const two_pi: number;
    function clamp(a: number, lo: number, hi: number): number;
    function swap<T>(a: T[], b: T[]): void;
    const isValid: typeof isFinite;
    function sqrt(n: number): number;
    function invSqrt(n: number): number;
    function degToRad(degrees: number): number;
    function radToDeg(radians: number): number;
    const Pow: (x: number, y: number) => number;
    const Sqrt: (x: number) => number;
    const Abs: (x: number) => number;
    const Cos: (x: number) => number;
    const Sin: (x: number) => number;
    const Acos: (x: number) => number;
    const Asin: (x: number) => number;
    const Atan2: (y: number, x: number) => number;
    const Random: () => number;
    const Min: (...values: number[]) => number;
    const Max: (...values: number[]) => number;
    function nextPowerOfTwo(x: number): number;
    function isPowerOfTwo(x: number): boolean;
    function randomRange(lo: number, hi: number): number;
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
        clone(): Vec2;
        setZero(): this;
        set(x: number, y: number): this;
        copy(other: XY): this;
        selfAdd(v: XY): this;
        selfAddXY(x: number, y: number): this;
        selfSub(v: XY): this;
        selfSubXY(x: number, y: number): this;
        selfMul(s: number): this;
        selfMulAdd(s: number, v: XY): this;
        selfMulSub(s: number, v: XY): this;
        dot(v: XY): number;
        cross(v: XY): number;
        length(): number;
        lengthSquared(): number;
        normalize(): number;
        selfNormalize(): this;
        selfRotate(radians: number): this;
        selfRotateCosSin(c: number, s: number): this;
        isValid(): boolean;
        selfCrossVS(s: number): this;
        selfCrossSV(s: number): this;
        selfMinV(v: XY): this;
        selfMaxV(v: XY): this;
        selfAbs(): this;
        selfNeg(): this;
        selfSkew(): this;
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
        clone(): TypedVec2;
        setZero(): this;
        set(x: number, y: number): this;
        copy(other: XY): this;
        selfAdd(v: XY): this;
        selfAddXY(x: number, y: number): this;
        selfSub(v: XY): this;
        selfSubXY(x: number, y: number): this;
        selfMul(s: number): this;
        selfMulAdd(s: number, v: XY): this;
        selfMulSub(s: number, v: XY): this;
        dot(v: XY): number;
        cross(v: XY): number;
        length(): number;
        lengthSquared(): number;
        normalize(): number;
        selfNormalize(): this;
        selfRotate(radians: number): this;
        selfRotateCosSin(c: number, s: number): this;
        isValid(): boolean;
        selfCrossVS(s: number): this;
        selfCrossSV(s: number): this;
        selfMinV(v: XY): this;
        selfMaxV(v: XY): this;
        selfAbs(): this;
        selfNeg(): this;
        selfSkew(): this;
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
        clone(): Vec3;
        setZero(): this;
        setXYZ(x: number, y: number, z: number): this;
        copy(other: XYZ): this;
        selfNeg(): this;
        selfAdd(v: XYZ): this;
        selfAddXYZ(x: number, y: number, z: number): this;
        selfSub(v: XYZ): this;
        selfSubXYZ(x: number, y: number, z: number): this;
        selfMul(s: number): this;
        static dotV3V3(a: XYZ, b: XYZ): number;
        static crossV3V3<T extends XYZ>(a: XYZ, b: XYZ, out: T): T;
    }
    class Mat22 {
        static readonly IDENTITY: Mat22;
        readonly ex: Vec2;
        readonly ey: Vec2;
        clone(): Mat22;
        static fromVV(c1: XY, c2: XY): Mat22;
        static fromSSSS(r1c1: number, r1c2: number, r2c1: number, r2c2: number): Mat22;
        static fromAngle(radians: number): Mat22;
        setSSSS(r1c1: number, r1c2: number, r2c1: number, r2c2: number): this;
        setVV(c1: XY, c2: XY): this;
        setAngle(radians: number): this;
        copy(other: Mat22): this;
        setIdentity(): this;
        setZero(): this;
        getAngle(): number;
        getInverse(out: Mat22): Mat22;
        solve<T extends XY>(b_x: number, b_y: number, out: T): T;
        selfAbs(): this;
        selfInv(): this;
        selfAddM(M: Mat22): this;
        selfSubM(M: Mat22): this;
        static absM(M: Mat22, out: Mat22): Mat22;
        static mulMV<T extends XY>(M: Mat22, v: XY, out: T): T;
        static mulTMV<T extends XY>(M: Mat22, v: XY, out: T): T;
        static addMM(A: Mat22, B: Mat22, out: Mat22): Mat22;
        static mulMM(A: Mat22, B: Mat22, out: Mat22): Mat22;
        static mulTMM(A: Mat22, B: Mat22, out: Mat22): Mat22;
    }
    class Mat33 {
        static readonly IDENTITY: Mat33;
        readonly data: Float32Array;
        readonly ex: Vec3;
        readonly ey: Vec3;
        readonly ez: Vec3;
        clone(): Mat33;
        setVVV(c1: XYZ, c2: XYZ, c3: XYZ): this;
        copy(other: Mat33): this;
        setIdentity(): this;
        setZero(): this;
        selfAddM(M: Mat33): this;
        solve33<T extends XYZ>(b_x: number, b_y: number, b_z: number, out: T): T;
        solve22<T extends XY>(b_x: number, b_y: number, out: T): T;
        getInverse22(M: Mat33): void;
        getSymInverse33(M: Mat33): void;
        static mulM33V3<T extends XYZ>(A: Mat33, v: XYZ, out: T): T;
        static mulM33XYZ<T extends XYZ>(A: Mat33, x: number, y: number, z: number, out: T): T;
        static mulM33V2<T extends XY>(A: Mat33, v: XY, out: T): T;
        static mulM33XY<T extends XY>(A: Mat33, x: number, y: number, out: T): T;
    }
    class Rot {
        static readonly IDENTITY: Rot;
        s: number;
        c: number;
        constructor(angle?: number);
        clone(): Rot;
        copy(other: Rot): this;
        setAngle(angle: number): this;
        setIdentity(): this;
        getAngle(): number;
        getXAxis<T extends XY>(out: T): T;
        getYAxis<T extends XY>(out: T): T;
        static mulRR(q: Rot, r: Rot, out: Rot): Rot;
        static mulTRR(q: Rot, r: Rot, out: Rot): Rot;
        static mulRV<T extends XY>(q: Rot, v: XY, out: T): T;
        static mulTRV<T extends XY>(q: Rot, v: XY, out: T): T;
    }
    class Transform {
        static readonly IDENTITY: Transform;
        readonly p: Vec2;
        readonly q: Rot;
        clone(): Transform;
        copy(other: Transform): this;
        setIdentity(): this;
        setPositionRotation(position: XY, q: Rot): this;
        setPositionAngle(pos: XY, a: number): this;
        setPosition(position: XY): this;
        setPositionXY(x: number, y: number): this;
        setRotation(rotation: Rot): this;
        setRotationAngle(radians: number): this;
        getPosition(): Vec2;
        getRotation(): Rot;
        getRotationAngle(): number;
        getAngle(): number;
        static mulXV<T extends XY>(T: Transform, v: XY, out: T): T;
        static mulTXV<T extends XY>(T: Transform, v: XY, out: T): T;
        static mulXX(A: Transform, B: Transform, out: Transform): Transform;
        static mulTXX(A: Transform, B: Transform, out: Transform): Transform;
    }
    class Sweep {
        readonly localCenter: Vec2;
        readonly c0: Vec2;
        readonly c: Vec2;
        a0: number;
        a: number;
        alpha0: number;
        clone(): Sweep;
        copy(other: Sweep): this;
        getTransform(xf: Transform, beta: number): Transform;
        advance(alpha: number): void;
        normalize(): void;
    }
}
declare namespace b2 {
    class DistanceProxy {
        readonly buffer: Vec2[];
        vertices: Vec2[];
        count: number;
        radius: number;
        copy(other: DistanceProxy): this;
        reset(): DistanceProxy;
        setShape(shape: Shape, index: number): void;
        setVerticesRadius(vertices: Vec2[], count: number, radius: number): void;
        getSupport(d: Vec2): number;
        getSupportVertex(d: Vec2): Vec2;
        getVertexCount(): number;
        getVertex(index: number): Vec2;
    }
    class SimplexCache {
        metric: number;
        count: number;
        readonly indexA: [number, number, number];
        readonly indexB: [number, number, number];
        reset(): SimplexCache;
    }
    class DistanceInput {
        readonly proxyA: DistanceProxy;
        readonly proxyB: DistanceProxy;
        readonly transformA: Transform;
        readonly transformB: Transform;
        useRadii: boolean;
        reset(): DistanceInput;
    }
    class DistanceOutput {
        readonly pointA: Vec2;
        readonly pointB: Vec2;
        distance: number;
        iterations: number;
        reset(): DistanceOutput;
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
    function gjkReset(): void;
    class SimplexVertex {
        readonly wA: Vec2;
        readonly wB: Vec2;
        readonly w: Vec2;
        a: number;
        indexA: number;
        indexB: number;
        copy(other: SimplexVertex): SimplexVertex;
    }
    class Simplex {
        readonly v1: SimplexVertex;
        readonly v2: SimplexVertex;
        readonly v3: SimplexVertex;
        readonly vertices: SimplexVertex[];
        count: number;
        constructor();
        readCache(cache: SimplexCache, proxyA: DistanceProxy, transformA: Transform, proxyB: DistanceProxy, transformB: Transform): void;
        writeCache(cache: SimplexCache): void;
        getSearchDirection(out: Vec2): Vec2;
        getClosestPoint(out: Vec2): Vec2;
        getWitnessPoints(pA: Vec2, pB: Vec2): void;
        getMetric(): number;
        solve2(): void;
        solve3(): void;
        private static s_e12;
        private static s_e13;
        private static s_e23;
    }
    function distance(output: DistanceOutput, cache: SimplexCache, input: DistanceInput): void;
    function shapeCast(output: ShapeCastOutput, input: ShapeCastInput): boolean;
}
declare namespace b2 {
    class Timer {
        start: number;
        reset(): Timer;
        getMilliseconds(): number;
    }
    class Counter {
        count: number;
        minCount: number;
        maxCount: number;
        increment(): void;
        decrement(): void;
    }
}
declare namespace b2 {
    enum ContactFeatureType {
        Vertex = 0,
        Face = 1
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
        copy(o: ContactID): ContactID;
        clone(): ContactID;
        key: number;
    }
    class ManifoldPoint {
        readonly localPoint: Vec2;
        normalImpulse: number;
        tangentImpulse: number;
        readonly id: ContactID;
        static makeArray(length: number): ManifoldPoint[];
        reset(): void;
        copy(o: ManifoldPoint): ManifoldPoint;
    }
    enum ManifoldType {
        Unknown = -1,
        Circles = 0,
        FaceA = 1,
        FaceB = 2
    }
    class Manifold {
        readonly points: ManifoldPoint[];
        readonly localNormal: Vec2;
        readonly localPoint: Vec2;
        type: ManifoldType;
        pointCount: number;
        reset(): void;
        copy(o: Manifold): Manifold;
        clone(): Manifold;
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
        initialize(manifold: Manifold, xfA: Transform, radiusA: number, xfB: Transform, radiusB: number): void;
    }
    enum PointState {
        NullState = 0,
        AddState = 1,
        PersistState = 2,
        RemoveState = 3
    }
    function getPointStates(state1: PointState[], state2: PointState[], manifold1: Manifold, manifold2: Manifold): void;
    class ClipVertex {
        readonly v: Vec2;
        readonly id: ContactID;
        static makeArray(length: number): ClipVertex[];
        copy(other: ClipVertex): ClipVertex;
    }
    class RayCastInput {
        readonly p1: Vec2;
        readonly p2: Vec2;
        maxFraction: number;
        copy(o: RayCastInput): RayCastInput;
    }
    class RayCastOutput {
        readonly normal: Vec2;
        fraction: number;
        copy(o: RayCastOutput): RayCastOutput;
    }
    class AABB {
        readonly lowerBound: Vec2;
        readonly upperBound: Vec2;
        private readonly cache_center;
        private readonly cache_extent;
        copy(o: AABB): AABB;
        isValid(): boolean;
        getCenter(): Vec2;
        getExtents(): Vec2;
        getPerimeter(): number;
        combine1(aabb: AABB): AABB;
        combine2(aabb1: AABB, aab: AABB): AABB;
        static combine(aabb1: AABB, aab: AABB, out: AABB): AABB;
        contains(aabb: AABB): boolean;
        rayCast(output: RayCastOutput, input: RayCastInput): boolean;
        testContain(point: XY): boolean;
        testOverlap(other: AABB): boolean;
    }
    function testOverlapAABB(a: AABB, b: AABB): boolean;
    function clipSegmentToLine(vOut: [ClipVertex, ClipVertex], vIn: [ClipVertex, ClipVertex], normal: Vec2, offset: number, vertexIndexA: number): number;
    function testOverlapShape(shapeA: Shape, indexA: number, shapeB: Shape, indexB: number, xfA: Transform, xfB: Transform): boolean;
}
declare namespace b2 {
    class MassData {
        mass: number;
        readonly center: Vec2;
        I: number;
    }
    enum ShapeType {
        Unknown = -1,
        CircleShape = 0,
        EdgeShape = 1,
        PolygonShape = 2,
        ChainShape = 3,
        ShapeTypeCount = 4
    }
    abstract class Shape {
        readonly type: ShapeType;
        radius: number;
        constructor(type: ShapeType, radius: number);
        abstract clone(): Shape;
        copy(other: Shape): Shape;
        getType(): ShapeType;
        abstract getChildCount(): number;
        abstract testPoint(xf: Transform, p: XY): boolean;
        abstract computeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number;
        abstract rayCast(output: RayCastOutput, input: RayCastInput, transform: Transform, childIndex: number): boolean;
        abstract computeAABB(aabb: AABB, xf: Transform, childIndex: number): void;
        abstract computeMass(massData: MassData, density: number): void;
        abstract setupDistanceProxy(proxy: DistanceProxy, index: number): void;
        abstract computeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number;
        abstract dump(log: (format: string, ...args: any[]) => void): void;
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
    function toiReset(): void;
    class TOIInput {
        readonly proxyA: DistanceProxy;
        readonly proxyB: DistanceProxy;
        readonly sweepA: Sweep;
        readonly sweepB: Sweep;
        tMax: number;
    }
    enum TOIOutputState {
        Unknown = 0,
        Failed = 1,
        Overlapped = 2,
        Touching = 3,
        Separated = 4
    }
    class TOIOutput {
        state: TOIOutputState;
        t: number;
    }
    enum SeparationFunctionType {
        Unknown = -1,
        Points = 0,
        FaceA = 1,
        FaceB = 2
    }
    class SeparationFunction {
        proxyA: DistanceProxy;
        proxyB: DistanceProxy;
        readonly sweepA: Sweep;
        readonly sweepB: Sweep;
        type: SeparationFunctionType;
        readonly localPoint: Vec2;
        readonly axis: Vec2;
        initialize(cache: SimplexCache, proxyA: DistanceProxy, sweepA: Sweep, proxyB: DistanceProxy, sweepB: Sweep, t1: number): number;
        findMinSeparation(indexA: [number], indexB: [number], t: number): number;
        evaluate(indexA: number, indexB: number, t: number): number;
    }
    function timeOfImpact(output: TOIOutput, input: TOIInput): void;
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
        clone(): Color;
        copy(other: RGBA): this;
        isEqual(color: RGBA): boolean;
        isZero(): boolean;
        set(r: number, g: number, b: number, a?: number): void;
        setByteRGB(r: number, g: number, b: number): this;
        setByteRGBA(r: number, g: number, b: number, a: number): this;
        setRGB(rr: number, gg: number, bb: number): this;
        setRGBA(rr: number, gg: number, bb: number, aa: number): this;
        selfAdd(color: RGBA): this;
        add<T extends RGBA>(color: RGBA, out: T): T;
        selfSub(color: RGBA): this;
        sub<T extends RGBA>(color: RGBA, out: T): T;
        selfMul(s: number): this;
        mul<T extends RGBA>(s: number, out: T): T;
        mix(mixColor: RGBA, strength: number): void;
        static mixColors(colorA: RGBA, colorB: RGBA, strength: number): void;
        makeStyleString(alpha?: number): string;
        static makeStyleString(r: number, g: number, b: number, a?: number): string;
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
        clone(): TypedColor;
        copy(other: RGBA): this;
        isEqual(color: RGBA): boolean;
        isZero(): boolean;
        set(r: number, g: number, b: number, a?: number): void;
        setByteRGB(r: number, g: number, b: number): this;
        setByteRGBA(r: number, g: number, b: number, a: number): this;
        setRGB(rr: number, gg: number, bb: number): this;
        setRGBA(rr: number, gg: number, bb: number, aa: number): this;
        selfAdd(color: RGBA): this;
        add<T extends RGBA>(color: RGBA, out: T): T;
        selfSub(color: RGBA): this;
        sub<T extends RGBA>(color: RGBA, out: T): T;
        selfMul(s: number): this;
        mul<T extends RGBA>(s: number, out: T): T;
        mix(mixColor: RGBA, strength: number): void;
        makeStyleString(alpha?: number): string;
    }
    enum DrawFlags {
        None = 0,
        ShapeBit = 1,
        JointBit = 2,
        AABBBit = 4,
        PairBit = 8,
        CenterOfMassBit = 16,
        ParticleBit = 32,
        ControllerBit = 64,
        All = 63
    }
    class Draw {
        drawFlags: DrawFlags;
        setFlags(flags: DrawFlags): void;
        getFlags(): DrawFlags;
        appendFlags(flags: DrawFlags): void;
        clearFlags(flags: DrawFlags): void;
        pushTransform(xf: Transform): void;
        popTransform(xf: Transform): void;
        drawPolygon(vertices: XY[], vertexCount: number, color: RGBA): void;
        drawSolidPolygon(vertices: XY[], vertexCount: number, color: RGBA): void;
        drawCircle(center: XY, radius: number, color: RGBA): void;
        drawSolidCircle(center: XY, radius: number, axis: XY, color: RGBA): void;
        drawParticles(centers: XY[], radius: number, colors: RGBA[], count: number): void;
        drawSegment(p1: XY, p2: XY, color: RGBA): void;
        drawTransform(xf: Transform): void;
        drawPoint(p: XY, size: number, color: RGBA): void;
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
        reset(): this;
    }
    class TimeStep {
        dt: number;
        inv_dt: number;
        dtRatio: number;
        velocityIterations: number;
        positionIterations: number;
        particleIterations: number;
        warmStarting: boolean;
        copy(step: TimeStep): TimeStep;
    }
    class Position {
        readonly c: Vec2;
        a: number;
        static makeArray(length: number): Position[];
    }
    class Velocity {
        readonly v: Vec2;
        w: number;
        static makeArray(length: number): Velocity[];
    }
    class SolverData {
        readonly step: TimeStep;
        positions: Position[];
        velocities: Velocity[];
    }
}
declare namespace b2 {
    class EdgeShape extends Shape {
        readonly vertex0: Vec2;
        readonly vertex1: Vec2;
        readonly vertex2: Vec2;
        readonly vertex3: Vec2;
        oneSided: boolean;
        constructor();
        setOneSided(v0: XY, v1: XY, v2: XY, v3: XY): EdgeShape;
        setTwoSided(v1: XY, v2: XY): EdgeShape;
        clone(): EdgeShape;
        copy(other: EdgeShape): EdgeShape;
        getChildCount(): number;
        testPoint(xf: Transform, p: XY): boolean;
        private static ComputeDistance_s_v1;
        private static ComputeDistance_s_v2;
        private static ComputeDistance_s_d;
        private static ComputeDistance_s_s;
        computeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number;
        private static rayCast_s_p1;
        private static rayCast_s_p2;
        private static rayCast_s_d;
        private static rayCast_s_e;
        private static rayCast_s_q;
        private static rayCast_s_r;
        rayCast(output: RayCastOutput, input: RayCastInput, xf: Transform, childIndex: number): boolean;
        private static computeAABB_s_v1;
        private static computeAABB_s_v2;
        computeAABB(aabb: AABB, xf: Transform, childIndex: number): void;
        computeMass(massData: MassData, density: number): void;
        setupDistanceProxy(proxy: DistanceProxy, index: number): void;
        computeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number;
        dump(log: (format: string, ...args: any[]) => void): void;
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
        prevBody: ControllerEdge;
        nextBody: ControllerEdge;
        prevController: ControllerEdge;
        nextController: ControllerEdge;
        constructor(controller: Controller, body: Body);
    }
    /**
     * Base class for controllers. Controllers are a convience for
     * encapsulating common per-step functionality.
     */
    abstract class Controller {
        bodyList: ControllerEdge;
        bodyCount: number;
        prev: Controller;
        next: Controller;
        /**
         * Controllers override this to implement per-step functionality.
         */
        abstract step(step: TimeStep): void;
        /**
         * Controllers override this to provide debug drawing.
         */
        draw(debugDraw: Draw): void;
        /**
         * Adds a body to the controller list.
         */
        addBody(body: Body): void;
        /**
         * Removes a body from the controller list.
         */
        removeBody(body: Body): void;
        /**
         * Removes all bodies from the controller list.
         */
        clear(): void;
    }
}
declare namespace b2 {
    function mixFriction(friction1: number, friction2: number): number;
    function mixRestitution(restitution1: number, restitution2: number): number;
    function mixRestitutionThreshold(threshold1: number, threshold2: number): number;
    class ContactEdge {
        private _other;
        other: Body;
        readonly contact: Contact;
        prev: ContactEdge;
        next: ContactEdge;
        constructor(contact: Contact);
        reset(): void;
    }
    abstract class Contact<A extends Shape = Shape, B extends Shape = Shape> {
        islandFlag: boolean;
        touchingFlag: boolean;
        enabledFlag: boolean;
        filterFlag: boolean;
        bulletHitFlag: boolean;
        toiFlag: boolean;
        prev: Contact;
        next: Contact;
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
        getManifold(): Manifold;
        getWorldManifold(worldManifold: WorldManifold): void;
        isTouching(): boolean;
        setEnabled(flag: boolean): void;
        isEnabled(): boolean;
        getNext(): Contact;
        getFixtureA(): Fixture;
        getChildIndexA(): number;
        getShapeA(): A;
        getFixtureB(): Fixture;
        getChildIndexB(): number;
        getShapeB(): B;
        abstract evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
        glagForFiltering(): void;
        setFriction(friction: number): void;
        getFriction(): number;
        resetFriction(): void;
        setRestitution(restitution: number): void;
        getRestitution(): number;
        resetRestitution(): void;
        setRestitutionThreshold(threshold: number): void;
        getRestitutionThreshold(): number;
        resetRestitutionThreshold(): void;
        setTangentSpeed(speed: number): void;
        getTangentSpeed(): number;
        reset(fixtureA: Fixture, indexA: number, fixtureB: Fixture, indexB: number): void;
        update(listener: ContactListener): void;
        private static computeTOI_s_input;
        private static computeTOI_s_output;
        computeTOI(sweepA: Sweep, sweepB: Sweep): number;
    }
}
declare namespace b2 {
    let gBlockSolve: boolean;
    class VelocityConstraintPoint {
        readonly rA: Vec2;
        readonly rB: Vec2;
        normalImpulse: number;
        tangentImpulse: number;
        normalMass: number;
        tangentMass: number;
        velocityBias: number;
        static makeArray(length: number): VelocityConstraintPoint[];
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
        static makeArray(length: number): ContactVelocityConstraint[];
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
        static makeArray(length: number): ContactPositionConstraint[];
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
        initialize(pc: ContactPositionConstraint, xfA: Transform, xfB: Transform, index: number): void;
    }
    class ContactSolver {
        readonly step: TimeStep;
        positions: Position[];
        velocities: Velocity[];
        readonly positionConstraints: ContactPositionConstraint[];
        readonly velocityConstraints: ContactVelocityConstraint[];
        contacts: Contact[];
        count: number;
        initialize(def: ContactSolverDef): ContactSolver;
        private static initializeVelocityConstraints_s_xfA;
        private static initializeVelocityConstraints_s_xfB;
        private static initializeVelocityConstraints_s_worldManifold;
        initializeVelocityConstraints(): void;
        private static warmStart_s_P;
        warmStart(): void;
        private static solveVelocityConstraints_s_dv;
        private static solveVelocityConstraints_s_dv1;
        private static solveVelocityConstraints_s_dv2;
        private static solveVelocityConstraints_s_P;
        private static solveVelocityConstraints_s_a;
        private static solveVelocityConstraints_s_b;
        private static solveVelocityConstraints_s_x;
        private static solveVelocityConstraints_s_d;
        private static solveVelocityConstraints_s_P1;
        private static solveVelocityConstraints_s_P2;
        private static solveVelocityConstraints_s_P1P2;
        solveVelocityConstraints(): void;
        storeImpulses(): void;
        private static solvePositionConstraints_s_xfA;
        private static solvePositionConstraints_s_xfB;
        private static solvePositionConstraints_s_psm;
        private static solvePositionConstraints_s_rA;
        private static solvePositionConstraints_s_rB;
        private static solvePositionConstraints_s_P;
        solvePositionConstraints(): boolean;
        private static solveTOIPositionConstraints_s_xfA;
        private static solveTOIPositionConstraints_s_xfB;
        private static solveTOIPositionConstraints_s_psm;
        private static solveTOIPositionConstraints_s_rA;
        private static solveTOIPositionConstraints_s_rB;
        private static solveTOIPositionConstraints_s_P;
        solveTOIPositionConstraints(toiIndexA: number, toiIndexB: number): boolean;
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
        clone(): Filter;
        copy(other: IFilter): this;
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
        reset(): void;
        touch(): void;
        private static synchronize_s_aabb1;
        private static synchronize_s_aab;
        private static synchronize_s_displacement;
        synchronize(transform1: Transform, transform2: Transform): void;
    }
    class Fixture {
        density: number;
        next: Fixture;
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
        reset(): void;
        getType(): ShapeType;
        getShape(): Shape;
        setSensor(sensor: boolean): void;
        setFilterData(filter: Filter): void;
        getFilterData(): Filter;
        refilter(): void;
        getBody(): Body;
        testPoint(p: XY): boolean;
        computeDistance(p: Vec2, normal: Vec2, childIndex: number): number;
        rayCast(output: RayCastOutput, input: RayCastInput, childIndex: number): boolean;
        getMassData(massData?: MassData): MassData;
        setDensity(density: number): void;
        getDensity(): number;
        getFriction(): number;
        setFriction(friction: number): void;
        getRestitution(): number;
        setRestitution(restitution: number): void;
        getRestitutionThreshold(): number;
        setRestitutionThreshold(threshold: number): void;
        getAABB(childIndex: number): AABB;
        dump(log: (format: string, ...args: any[]) => void, bodyIndex: number): void;
        createProxies(): void;
        destroyProxies(): void;
        touchProxies(): void;
        synchronizeProxies(transform1: Transform, transform2: Transform): void;
    }
}
declare namespace b2 {
    enum JointType {
        UnknownJoint = 0,
        RevoluteJoint = 1,
        PrismaticJoint = 2,
        DistanceJoint = 3,
        PulleyJoint = 4,
        MouseJoint = 5,
        GearJoint = 6,
        WheelJoint = 7,
        WeldJoint = 8,
        FrictionJoint = 9,
        RopeJoint = 10,
        MotorJoint = 11,
        AreaJoint = 12
    }
    class Jacobian {
        readonly linear: Vec2;
        angularA: number;
        angularB: number;
        setZero(): Jacobian;
        set(x: XY, a1: number, a2: number): Jacobian;
    }
    class JointEdge {
        private _other;
        other: Body;
        readonly joint: Joint;
        prev: JointEdge;
        next: JointEdge;
        constructor(joint: Joint);
        reset(): void;
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
    function linearStiffness(def: {
        stiffness: number;
        damping: number;
    }, frequencyHertz: number, dampingRatio: number, bodyA: Body, bodyB: Body): void;
    function angularStiffness(def: {
        stiffness: number;
        damping: number;
    }, frequencyHertz: number, dampingRatio: number, bodyA: Body, bodyB: Body): void;
    abstract class Joint {
        readonly type: JointType;
        prev: Joint;
        next: Joint;
        readonly edgeA: JointEdge;
        readonly edgeB: JointEdge;
        bodyA: Body;
        bodyB: Body;
        index: number;
        islandFlag: boolean;
        collideConnected: boolean;
        userData: any;
        constructor(def: IJointDef);
        getType(): JointType;
        getBodyA(): Body;
        getBodyB(): Body;
        abstract getAnchorA<T extends XY>(out: T): T;
        abstract getAnchorB<T extends XY>(out: T): T;
        abstract getReactionForce<T extends XY>(inv_dt: number, out: T): T;
        abstract getReactionTorque(inv_dt: number): number;
        isEnabled(): boolean;
        getCollideConnected(): boolean;
        dump(log: (format: string, ...args: any[]) => void): void;
        shiftOrigin(newOrigin: XY): void;
        private static draw_s_p1;
        private static draw_s_p2;
        private static draw_s_color;
        private static draw_s_c;
        draw(draw: Draw): void;
        abstract initVelocityConstraints(data: SolverData): void;
        abstract solveVelocityConstraints(data: SolverData): void;
        abstract solvePositionConstraints(data: SolverData): boolean;
    }
}
declare namespace b2 {
    class DestructionListener {
        sayGoodbyeJoint(joint: Joint): void;
        sayGoodbyeFixture(fixture: Fixture): void;
        sayGoodbyeParticleGroup(group: ParticleGroup): void;
        sayGoodbyeParticle(system: ParticleSystem, index: number): void;
    }
    class ContactFilter {
        ShouldCollide(fixtureA: Fixture, fixtureB: Fixture): boolean;
        shouldCollideFixtureParticle(fixture: Fixture, system: ParticleSystem, index: number): boolean;
        shouldCollideParticleParticle(system: ParticleSystem, indexA: number, indexB: number): boolean;
        static readonly defaultFilter: ContactFilter;
    }
    class ContactImpulse {
        normalImpulses: number[];
        tangentImpulses: number[];
        count: number;
    }
    class ContactListener {
        beginContact(contact: Contact): void;
        endContact(contact: Contact): void;
        beginContactFixtureParticle(system: ParticleSystem, contact: ParticleBodyContact): void;
        endContactFixtureParticle(system: ParticleSystem, contact: ParticleBodyContact): void;
        beginContactParticleParticle(system: ParticleSystem, contact: ParticleContact): void;
        endContactParticleParticle(system: ParticleSystem, contact: ParticleContact): void;
        preSolve(contact: Contact, oldManifold: Manifold): void;
        postSolve(contact: Contact, impulse: ContactImpulse): void;
        static readonly defaultListener: ContactListener;
    }
    class QueryCallback {
        reportFixture(fixture: Fixture): boolean;
        reportParticle(system: ParticleSystem, index: number): boolean;
        shouldQueryParticleSystem(system: ParticleSystem): boolean;
    }
    type QueryCallbackFunction = (fixture: Fixture) => boolean;
    class RayCastCallback {
        reportFixture(fixture: Fixture, point: Vec2, normal: Vec2, fraction: number): number;
        reportParticle(system: ParticleSystem, index: number, point: Vec2, normal: Vec2, fraction: number): number;
        shouldQueryParticleSystem(system: ParticleSystem): boolean;
    }
    type RayCastCallbackFunction = (fixture: Fixture, point: Vec2, normal: Vec2, fraction: number) => number;
}
declare namespace b2 {
    /**
     * The particle type. Can be combined with the | operator.
     */
    enum ParticleFlag {
        WaterParticle = 0,
        ZombieParticle = 2,
        WallParticle = 4,
        SpringParticle = 8,
        ElasticParticle = 16,
        ViscousParticle = 32,
        PowderParticle = 64,
        TensileParticle = 128,
        ColorMixingParticle = 256,
        DestructionListenerParticle = 512,
        BarrierParticle = 1024,
        StaticPressureParticle = 2048,
        ReactiveParticle = 4096,
        RepulsiveParticle = 8192,
        FixtureContactListenerParticle = 16384,
        ParticleContactListenerParticle = 32768,
        FixtureContactFilterParticle = 65536,
        ParticleContactFilterParticle = 131072
    }
    interface IParticleDef {
        flags?: ParticleFlag;
        position?: XY;
        velocity?: XY;
        color?: RGBA;
        lifetime?: number;
        userData?: any;
        group?: ParticleGroup;
    }
    class ParticleDef implements IParticleDef {
        flags: ParticleFlag;
        readonly position: Vec2;
        readonly velocity: Vec2;
        readonly color: Color;
        lifetime: number;
        userData: any;
        group: ParticleGroup;
    }
    function CalculateParticleIterations(gravity: number, radius: number, timeStep: number): number;
    class ParticleHandle {
        index: number;
        getIndex(): number;
        setIndex(index: number): void;
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
        readonly moveBuffer: Array<TreeNode<T>>;
        pairCount: number;
        readonly pairBuffer: Array<Pair<T>>;
        createProxy(aabb: AABB, userData: T): TreeNode<T>;
        destroyProxy(proxy: TreeNode<T>): void;
        moveProxy(proxy: TreeNode<T>, aabb: AABB, displacement: Vec2): void;
        touchProxy(proxy: TreeNode<T>): void;
        getProxyCount(): number;
        updatePairs(callback: (a: T, b: T) => void): void;
        query(aabb: AABB, callback: (node: TreeNode<T>) => boolean): void;
        queryPoint(point: XY, callback: (node: TreeNode<T>) => boolean): void;
        rayCast(input: RayCastInput, callback: (input: RayCastInput, node: TreeNode<T>) => number): void;
        getTreeHeight(): number;
        getTreeBalance(): number;
        getTreeQuality(): number;
        shiftOrigin(newOrigin: XY): void;
        bufferMove(proxy: TreeNode<T>): void;
        unBufferMove(proxy: TreeNode<T>): void;
    }
}
declare namespace b2 {
    class ChainShape extends Shape {
        vertices: Vec2[];
        count: number;
        readonly prevVertex: Vec2;
        readonly nextVertex: Vec2;
        constructor();
        createLoop(vertices: XY[]): ChainShape;
        createLoop(vertices: XY[], count: number): ChainShape;
        createLoop(vertices: number[]): ChainShape;
        private _createLoop;
        createChain(vertices: XY[], prevVertex: XY, nextVertex: XY): ChainShape;
        createChain(vertices: XY[], count: number, prevVertex: XY, nextVertex: XY): ChainShape;
        createChain(vertices: number[], prevVertex: XY, nextVertex: XY): ChainShape;
        private _createChain;
        clone(): ChainShape;
        copy(other: ChainShape): ChainShape;
        getChildCount(): number;
        getChildEdge(edge: EdgeShape, index: number): void;
        testPoint(xf: Transform, p: XY): boolean;
        private static ComputeDistance_s_edgeShape;
        computeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number;
        private static rayCast_s_edgeShape;
        rayCast(output: RayCastOutput, input: RayCastInput, xf: Transform, childIndex: number): boolean;
        private static computeAABB_s_v1;
        private static computeAABB_s_v2;
        private static computeAABB_s_lower;
        private static computeAABB_s_upper;
        computeAABB(aabb: AABB, xf: Transform, childIndex: number): void;
        computeMass(massData: MassData, density: number): void;
        setupDistanceProxy(proxy: DistanceProxy, index: number): void;
        computeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number;
        dump(log: (format: string, ...args: any[]) => void): void;
    }
}
declare namespace b2 {
    class CircleShape extends Shape {
        readonly p: Vec2;
        constructor(radius?: number);
        set(position: XY, radius?: number): this;
        clone(): CircleShape;
        copy(other: CircleShape): CircleShape;
        getChildCount(): number;
        private static TestPoint_s_center;
        private static TestPoint_s_d;
        testPoint(transform: Transform, p: XY): boolean;
        private static ComputeDistance_s_center;
        computeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number;
        private static rayCast_s_position;
        private static rayCast_s_s;
        private static rayCast_s_r;
        rayCast(output: RayCastOutput, input: RayCastInput, transform: Transform, childIndex: number): boolean;
        private static computeAABB_s_p;
        computeAABB(aabb: AABB, transform: Transform, childIndex: number): void;
        computeMass(massData: MassData, density: number): void;
        setupDistanceProxy(proxy: DistanceProxy, index: number): void;
        computeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number;
        dump(log: (format: string, ...args: any[]) => void): void;
    }
}
declare namespace b2 {
    function collideCircles(manifold: Manifold, circleA: CircleShape, xfA: Transform, circleB: CircleShape, xfB: Transform): void;
    function collidePolygonAndCircle(manifold: Manifold, polygonA: PolygonShape, xfA: Transform, circleB: CircleShape, xfB: Transform): void;
}
declare namespace b2 {
    function collideEdgeAndCircle(manifold: Manifold, edgeA: EdgeShape, xfA: Transform, circleB: CircleShape, xfB: Transform): void;
    function collideEdgeAndPolygon(manifold: Manifold, edgeA: EdgeShape, xfA: Transform, polygonB: PolygonShape, xfB: Transform): void;
}
declare namespace b2 {
    function collidePolygons(manifold: Manifold, polyA: PolygonShape, xfA: Transform, polyB: PolygonShape, xfB: Transform): void;
}
declare namespace b2 {
    class TreeNode<T> {
        readonly id: number;
        readonly aabb: AABB;
        private _userData;
        userData: T;
        parent: TreeNode<T>;
        child1: TreeNode<T>;
        child2: TreeNode<T>;
        height: number;
        moved: boolean;
        constructor(id?: number);
        reset(): void;
        isLeaf(): boolean;
    }
    class DynamicTree<T> {
        root: TreeNode<T>;
        freeList: TreeNode<T>;
        insertionCount: number;
        readonly stack: GrowableStack<TreeNode<T>>;
        static readonly sR: Vec2;
        static readonly sV: Vec2;
        static readonly absV: Vec2;
        static readonly segmentAABB: AABB;
        static readonly subInput: RayCastInput;
        static readonly combinedAABB: AABB;
        static readonly sAabb: AABB;
        query(aabb: AABB, callback: (node: TreeNode<T>) => boolean): void;
        queryPoint(point: XY, callback: (node: TreeNode<T>) => boolean): void;
        rayCast(input: RayCastInput, callback: (input: RayCastInput, node: TreeNode<T>) => number): void;
        static sNodeId: number;
        allocateNode(): TreeNode<T>;
        freeNode(node: TreeNode<T>): void;
        createProxy(aabb: AABB, userData: T): TreeNode<T>;
        destroyProxy(node: TreeNode<T>): void;
        private static moveProxy_s_fatAABB;
        private static moveProxy_s_hugeAABB;
        moveProxy(node: TreeNode<T>, aabb: AABB, displacement: Vec2): boolean;
        insertLeaf(leaf: TreeNode<T>): void;
        removeLeaf(leaf: TreeNode<T>): void;
        balance(A: TreeNode<T>): TreeNode<T>;
        getHeight(): number;
        private static getAreaNode;
        getAreaRatio(): number;
        static computeHeightNode<T>(node: TreeNode<T>): number;
        computeHeight(): number;
        validateStructure(node: TreeNode<T>): void;
        validateMetrics(node: TreeNode<T>): void;
        validate(): void;
        private static getMaxBalanceNode;
        getMaxBalance(): number;
        rebuildBottomUp(): void;
        private static shiftOriginNode;
        shiftOrigin(newOrigin: XY): void;
    }
}
declare namespace b2 {
    class PolygonShape extends Shape {
        readonly centroid: Vec2;
        vertices: Vec2[];
        normals: Vec2[];
        count: number;
        constructor();
        clone(): PolygonShape;
        copy(other: PolygonShape): PolygonShape;
        getChildCount(): number;
        private static Set_s_r;
        private static Set_s_v;
        set(vertices: XY[]): PolygonShape;
        set(vertices: XY[], count: number): PolygonShape;
        set(vertices: number[]): PolygonShape;
        _set(vertices: (index: number) => XY, count: number): PolygonShape;
        setAsBox(hx: number, hy: number, center?: XY, angle?: number): PolygonShape;
        private static testPointSPLocal;
        testPoint(xf: Transform, p: XY): boolean;
        private static computeDistance_s_pLocal;
        private static computeDistance_s_normalForMaxDistance;
        private static computeDistance_s_minDistance;
        private static computeDistance_s_distance;
        computeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number;
        private static rayCast_s_p1;
        private static rayCast_s_p2;
        private static rayCast_s_d;
        rayCast(output: RayCastOutput, input: RayCastInput, xf: Transform, childIndex: number): boolean;
        private static computeAABB_s_v;
        computeAABB(aabb: AABB, xf: Transform, childIndex: number): void;
        private static computeMass_s_center;
        private static computeMass_s_s;
        private static computeMass_s_e1;
        private static computeMass_s_e2;
        computeMass(massData: MassData, density: number): void;
        private static validate_s_e;
        private static validate_s_v;
        validate(): boolean;
        setupDistanceProxy(proxy: DistanceProxy, index: number): void;
        private static ComputeSubmergedArea_s_normalL;
        private static ComputeSubmergedArea_s_md;
        private static ComputeSubmergedArea_s_intoVec;
        private static ComputeSubmergedArea_s_outoVec;
        private static ComputeSubmergedArea_s_center;
        computeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number;
        dump(log: (format: string, ...args: any[]) => void): void;
        private static computeCentroid_s_s;
        private static computeCentroid_s_p1;
        private static computeCentroid_s_p2;
        private static computeCentroid_s_p3;
        private static computeCentroid_s_e1;
        private static computeCentroid_s_e2;
        static computeCentroid(vs: Vec2[], count: number, out: Vec2): Vec2;
    }
}
declare namespace b2 {
    class BlockAllocator {
    }
}
declare namespace b2 {
    class GrowableStack<T> {
        stack: Array<T>;
        count: number;
        constructor(N: number);
        reset(): this;
        push(element: T): void;
        pop(): T;
        getCount(): number;
    }
}
declare namespace b2 {
    function log(message: string, ...args: any[]): void;
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
        step(step: TimeStep): void;
        draw(debugDraw: Draw): void;
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
        step(step: TimeStep): void;
        private static step_s_dtA;
        draw(draw: Draw): void;
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
        step(step: TimeStep): void;
        draw(draw: Draw): void;
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
        step(step: TimeStep): void;
        private static step_s_f;
        draw(draw: Draw): void;
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
        step(step: TimeStep): void;
        private static step_s_damping;
        draw(draw: Draw): void;
        /**
         * Sets damping independantly along the x and y axes
         */
        setAxisAligned(xDamping: number, yDamping: number): void;
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
        addBody(body: Body): void;
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
        getAnchorA<T extends XY>(out: T): T;
        getAnchorB<T extends XY>(out: T): T;
        getReactionForce<T extends XY>(inv_dt: number, out: T): T;
        getReactionTorque(inv_dt: number): number;
        setStiffness(stiffness: number): void;
        getStiffness(): number;
        setDamping(damping: number): void;
        getDamping(): number;
        dump(log: (format: string, ...args: any[]) => void): void;
        initVelocityConstraints(data: SolverData): void;
        solveVelocityConstraints(data: SolverData): void;
        solvePositionConstraints(data: SolverData): boolean;
    }
}
declare namespace b2 {
    enum BodyType {
        Unknown = -1,
        StaticBody = 0,
        KinematicBody = 1,
        DynamicBody = 2
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
        prev: Body;
        next: Body;
        fixtureList: Fixture;
        fixtureCount: number;
        jointList: JointEdge;
        contactList: ContactEdge;
        mass: number;
        invMass: number;
        I: number;
        invI: number;
        linearDamping: number;
        angularDamping: number;
        gravityScale: number;
        sleepTime: number;
        userData: any;
        controllerList: ControllerEdge;
        controllerCount: number;
        constructor(bd: IBodyDef, world: World);
        createFixture(def: IFixtureDef): Fixture;
        createFixture(shape: Shape): Fixture;
        createFixture(shape: Shape, density: number): Fixture;
        addFixture(fixture: Fixture): void;
        createFixtureDef(def: IFixtureDef): Fixture;
        private static createFixtureShapeDensity_s_def;
        createFixtureShapeDensity(shape: Shape, density?: number): Fixture;
        destroyFixture(fixture: Fixture): void;
        setTransformVec(position: XY, angle: number): void;
        setTransformXY(x: number, y: number, angle: number): void;
        setTransform(xf: Transform): void;
        getTransform(): Transform;
        getPosition(): Vec2;
        setPosition(position: XY): void;
        setPositionXY(x: number, y: number): void;
        getAngle(): number;
        setAngle(angle: number): void;
        getWorldCenter(): Vec2;
        getLocalCenter(): Vec2;
        setLinearVelocity(v: XY): void;
        getLinearVelocity(): Vec2;
        setAngularVelocity(w: number): void;
        getAngularVelocity(): number;
        getDefinition(bd: BodyDef): BodyDef;
        applyForce(force: XY, point: XY, wake?: boolean): void;
        applyForceToCenter(force: XY, wake?: boolean): void;
        applyTorque(torque: number, wake?: boolean): void;
        applyLinearImpulse(impulse: XY, point: XY, wake?: boolean): void;
        applyLinearImpulseToCenter(impulse: XY, wake?: boolean): void;
        applyAngularImpulse(impulse: number, wake?: boolean): void;
        getMass(): number;
        getInertia(): number;
        getMassData(data: MassData): MassData;
        private static setMassData_s_oldCenter;
        setMassData(massData: MassData): void;
        private static resetMassData_s_localCenter;
        private static resetMassData_s_oldCenter;
        private static resetMassData_s_massData;
        resetMassData(): void;
        getWorldPoint<T extends XY>(localPoint: XY, out: T): T;
        getWorldVector<T extends XY>(localVector: XY, out: T): T;
        getLocalPoint<T extends XY>(worldPoint: XY, out: T): T;
        getLocalVector<T extends XY>(worldVector: XY, out: T): T;
        getLinearVelocityFromWorldPoint<T extends XY>(worldPoint: XY, out: T): T;
        getLinearVelocityFromLocalPoint<T extends XY>(localPoint: XY, out: T): T;
        getLinearDamping(): number;
        setLinearDamping(linearDamping: number): void;
        getAngularDamping(): number;
        setAngularDamping(angularDamping: number): void;
        getGravityScale(): number;
        setGravityScale(scale: number): void;
        setType(type: BodyType): void;
        getType(): BodyType;
        setBullet(flag: boolean): void;
        isBullet(): boolean;
        setSleepingAllowed(flag: boolean): void;
        isSleepingAllowed(): boolean;
        setAwake(flag: boolean): void;
        isAwake(): boolean;
        setEnabled(flag: boolean): void;
        isEnabled(): boolean;
        setFixedRotation(flag: boolean): void;
        isFixedRotation(): boolean;
        getFixtureList(): Fixture;
        getJointList(): JointEdge;
        getContactList(): ContactEdge;
        getNext(): Body;
        getUserData(): any;
        setUserData(data: any): void;
        getWorld(): World;
        dump(log: (format: string, ...args: any[]) => void): void;
        private static synchronizeFixtures_s_xf1;
        synchronizeFixtures(): void;
        synchronizeTransform(): void;
        shouldCollide(other: Body): boolean;
        shouldCollideConnected(other: Body): boolean;
        advance(alpha: number): void;
        getControllerList(): ControllerEdge;
        getControllerCount(): number;
    }
}
declare namespace b2 {
    class ChainAndCircleContact extends Contact<ChainShape, CircleShape> {
        static create(): Contact;
        static destroy(contact: Contact): void;
        private static evaluate_s_edge;
        evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
    }
}
declare namespace b2 {
    class ChainAndPolygonContact extends Contact<ChainShape, PolygonShape> {
        static create(): Contact;
        static destroy(contact: Contact): void;
        private static evaluate_s_edge;
        evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
    }
}
declare namespace b2 {
    class CircleContact extends Contact<CircleShape, CircleShape> {
        static create(): Contact;
        static destroy(contact: Contact): void;
        evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
    }
}
declare namespace b2 {
    class ContactRegister {
        pool: Contact[];
        createFcn: (() => Contact);
        destroyFcn: ((contact: Contact) => void);
        primary: boolean;
    }
    class ContactFactory {
        readonly registers: ContactRegister[][];
        constructor();
        private addType;
        private initializeRegisters;
        create(fixtureA: Fixture, indexA: number, fixtureB: Fixture, indexB: number): Contact;
        destroy(contact: Contact): void;
    }
}
declare namespace b2 {
    class ContactManager {
        readonly broadPhase: BroadPhase<FixtureProxy>;
        contactList: Contact;
        contactCount: number;
        contactFilter: ContactFilter;
        contactListener: ContactListener;
        readonly contactFactory: ContactFactory;
        addPair(proxyA: FixtureProxy, proxyB: FixtureProxy): void;
        findNewContacts(): void;
        destroy(c: Contact): void;
        collide(): void;
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
        getAnchorA<T extends XY>(out: T): T;
        getAnchorB<T extends XY>(out: T): T;
        getReactionForce<T extends XY>(inv_dt: number, out: T): T;
        getReactionTorque(inv_dt: number): number;
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
        dump(log: (format: string, ...args: any[]) => void): void;
        private static initVelocityConstraints_s_P;
        initVelocityConstraints(data: SolverData): void;
        private static solveVelocityConstraints_s_vpA;
        private static solveVelocityConstraints_s_vpB;
        private static solveVelocityConstraints_s_P;
        solveVelocityConstraints(data: SolverData): void;
        private static solvePositionConstraints_s_P;
        solvePositionConstraints(data: SolverData): boolean;
        private static draw_s_pA;
        private static draw_s_pB;
        private static draw_s_axis;
        private static draw_s_c1;
        private static draw_s_c2;
        private static draw_s_c3;
        private static draw_s_c4;
        private static draw_s_pRest;
        private static draw_s_pMin;
        private static draw_s_pMax;
        draw(draw: Draw): void;
    }
}
declare namespace b2 {
    class EdgeAndCircleContact extends Contact<EdgeShape, CircleShape> {
        static create(): Contact;
        static destroy(contact: Contact): void;
        evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
    }
}
declare namespace b2 {
    class EdgeAndPolygonContact extends Contact<EdgeShape, PolygonShape> {
        static create(): Contact;
        static destroy(contact: Contact): void;
        evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
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
        initialize(bA: Body, bB: Body, anchor: Vec2): void;
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
        initVelocityConstraints(data: SolverData): void;
        private static solveVelocityConstraints_s_Cdot_v2;
        private static solveVelocityConstraints_s_impulseV;
        private static solveVelocityConstraints_s_oldImpulseV;
        solveVelocityConstraints(data: SolverData): void;
        solvePositionConstraints(data: SolverData): boolean;
        getAnchorA<T extends XY>(out: T): T;
        getAnchorB<T extends XY>(out: T): T;
        getReactionForce<T extends XY>(inv_dt: number, out: T): T;
        getReactionTorque(inv_dt: number): number;
        dump(log: (format: string, ...args: any[]) => void): void;
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
        private static initVelocityConstraints_s_u;
        private static initVelocityConstraints_s_rA;
        private static initVelocityConstraints_s_rB;
        private static initVelocityConstraints_s_rC;
        private static initVelocityConstraints_s_rD;
        initVelocityConstraints(data: SolverData): void;
        solveVelocityConstraints(data: SolverData): void;
        private static solvePositionConstraints_s_u;
        private static solvePositionConstraints_s_rA;
        private static solvePositionConstraints_s_rB;
        private static solvePositionConstraints_s_rC;
        private static solvePositionConstraints_s_rD;
        solvePositionConstraints(data: SolverData): boolean;
        getAnchorA<T extends XY>(out: T): T;
        getAnchorB<T extends XY>(out: T): T;
        getReactionForce<T extends XY>(inv_dt: number, out: T): T;
        getReactionTorque(inv_dt: number): number;
        dump(log: (format: string, ...args: any[]) => void): void;
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
        initialize(bodyCapacity: number, contactCapacity: number, jointCapacity: number, listener: ContactListener): void;
        clear(): void;
        addBody(body: Body): void;
        addContact(contact: Contact): void;
        addJoint(joint: Joint): void;
        private static s_timer;
        private static s_solverData;
        private static s_contactSolverDef;
        private static s_contactSolver;
        private static s_translation;
        solve(profile: Profile, step: TimeStep, gravity: Vec2, allowSleep: boolean): void;
        solveTOI(subStep: TimeStep, toiIndexA: number, toiIndexB: number): void;
        private static s_impulse;
        report(constraints: ContactVelocityConstraint[]): void;
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
        initialize(bA: Body, bB: Body): void;
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
        getAnchorA<T extends XY>(out: T): T;
        getAnchorB<T extends XY>(out: T): T;
        getReactionForce<T extends XY>(inv_dt: number, out: T): T;
        getReactionTorque(inv_dt: number): number;
        setLinearOffset(linearOffset: Vec2): void;
        getLinearOffset(): Vec2;
        setAngularOffset(angularOffset: number): void;
        getAngularOffset(): number;
        setMaxForce(force: number): void;
        getMaxForce(): number;
        setMaxTorque(torque: number): void;
        getMaxTorque(): number;
        initVelocityConstraints(data: SolverData): void;
        private static solveVelocityConstraints_s_Cdot_v2;
        private static solveVelocityConstraints_s_impulse_v2;
        private static solveVelocityConstraints_s_oldImpulse_v2;
        solveVelocityConstraints(data: SolverData): void;
        solvePositionConstraints(data: SolverData): boolean;
        dump(log: (format: string, ...args: any[]) => void): void;
    }
}
declare namespace b2 {
    interface IMouseJointDef extends IJointDef {
        target?: XY;
        maxForce?: number;
        stiffness?: number;
        damping?: number;
        localAnchorB?: XY;
    }
    class MouseJointDef extends JointDef implements IMouseJointDef {
        target: Vec2;
        maxForce: number;
        stiffness: number;
        damping: number;
        localAnchorB: Vec2;
        constructor();
        initialize(bB: Body, anchorB: XY, target: XY): void;
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
        setTarget(target: XY): void;
        getTarget(): Vec2;
        setMaxForce(maxForce: number): void;
        getMaxForce(): number;
        setStiffness(stiffness: number): void;
        getStiffness(): number;
        setDamping(damping: number): void;
        getDamping(): number;
        initVelocityConstraints(data: SolverData): void;
        private static solveVelocityConstraints_s_Cdot;
        private static solveVelocityConstraints_s_impulse;
        private static solveVelocityConstraints_s_oldImpulse;
        solveVelocityConstraints(data: SolverData): void;
        solvePositionConstraints(data: SolverData): boolean;
        getAnchorA<T extends XY>(out: T): T;
        getAnchorB<T extends XY>(out: T): T;
        getReactionForce<T extends XY>(inv_dt: number, out: T): T;
        getReactionTorque(inv_dt: number): number;
        dump(log: (format: string, ...args: any[]) => void): void;
        shiftOrigin(newOrigin: Vec2): void;
    }
}
declare namespace b2 {
    class PolygonAndCircleContact extends Contact<PolygonShape, CircleShape> {
        static create(): Contact;
        static destroy(contact: Contact): void;
        evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
    }
}
declare namespace b2 {
    class PolygonContact extends Contact<PolygonShape, PolygonShape> {
        static create(): Contact;
        static destroy(contact: Contact): void;
        evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void;
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
        initialize(bA: Body, bB: Body, anchor1: XY, anchor2: XY, axis: XY): void;
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
        private static initVelocityConstraints_s_d;
        private static initVelocityConstraints_s_P;
        initVelocityConstraints(data: SolverData): void;
        private static solveVelocityConstraints_s_P;
        private static solveVelocityConstraints_s_df;
        solveVelocityConstraints(data: SolverData): void;
        private static solvePositionConstraints_s_d;
        private static solvePositionConstraints_s_impulse;
        private static solvePositionConstraints_s_impulse1;
        private static solvePositionConstraints_s_P;
        solvePositionConstraints(data: SolverData): boolean;
        getAnchorA<T extends XY>(out: T): T;
        getAnchorB<T extends XY>(out: T): T;
        getReactionForce<T extends XY>(inv_dt: number, out: T): T;
        getReactionTorque(inv_dt: number): number;
        private static getJointTranslation_s_pA;
        private static getJointTranslation_s_pB;
        private static getJointTranslation_s_d;
        private static getJointTranslation_s_axis;
        getJointTranslation(): number;
        getJointSpeed(): number;
        isLimitEnabled(): boolean;
        setEnableLimit(flag: boolean): void;
        getLowerLimit(): number;
        getUpperLimit(): number;
        setLimits(lower: number, upper: number): void;
        isMotorEnabled(): boolean;
        setEnableMotor(flag: boolean): void;
        setMotorSpeed(speed: number): void;
        getMotorSpeed(): number;
        setMaxMotorForce(force: number): void;
        getMaxMotorForce(): number;
        getMotorForce(inv_dt: number): number;
        dump(log: (format: string, ...args: any[]) => void): void;
        private static draw_s_pA;
        private static draw_s_pB;
        private static draw_s_axis;
        private static draw_s_c1;
        private static draw_s_c2;
        private static draw_s_c3;
        private static draw_s_c4;
        private static draw_s_c5;
        private static draw_s_lower;
        private static draw_s_upper;
        private static draw_s_perp;
        draw(draw: Draw): void;
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
        initialize(bA: Body, bB: Body, groundA: XY, groundB: XY, anchorA: XY, anchorB: XY, r: number): void;
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
        private static initVelocityConstraints_s_PA;
        private static initVelocityConstraints_s_PB;
        initVelocityConstraints(data: SolverData): void;
        private static solveVelocityConstraints_s_vpA;
        private static solveVelocityConstraints_s_vpB;
        private static solveVelocityConstraints_s_PA;
        private static solveVelocityConstraints_s_PB;
        solveVelocityConstraints(data: SolverData): void;
        private static solvePositionConstraints_s_PA;
        private static solvePositionConstraints_s_PB;
        solvePositionConstraints(data: SolverData): boolean;
        getAnchorA<T extends XY>(out: T): T;
        getAnchorB<T extends XY>(out: T): T;
        getReactionForce<T extends XY>(inv_dt: number, out: T): T;
        getReactionTorque(inv_dt: number): number;
        getGroundAnchorA(): Vec2;
        getGroundAnchorB(): Vec2;
        private static getCurrentLengthA_s_p;
        getCurrentLengthA(): number;
        private static getCurrentLengthB_s_p;
        getCurrentLengthB(): number;
        dump(log: (format: string, ...args: any[]) => void): void;
        shiftOrigin(newOrigin: Vec2): void;
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
        initialize(bA: Body, bB: Body, anchor1: XY, anchor2: XY): void;
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
        private static initVelocityConstraints_s_P;
        initVelocityConstraints(data: SolverData): void;
        private static solveVelocityConstraints_s_Cdot_v2;
        private static solveVelocityConstraints_s_impulse_v2;
        solveVelocityConstraints(data: SolverData): void;
        private static solvePositionConstraints_s_C_v2;
        private static solvePositionConstraints_s_impulse;
        solvePositionConstraints(data: SolverData): boolean;
        getAnchorA<T extends XY>(out: T): T;
        getAnchorB<T extends XY>(out: T): T;
        getReactionForce<T extends XY>(inv_dt: number, out: T): T;
        getReactionTorque(inv_dt: number): number;
        getJointAngle(): number;
        getJointSpeed(): number;
        isMotorEnabled(): boolean;
        setEnableMotor(flag: boolean): void;
        getMotorTorque(inv_dt: number): number;
        getMotorSpeed(): number;
        setMaxMotorTorque(torque: number): void;
        getMaxMotorTorque(): number;
        isLimitEnabled(): boolean;
        setEnableLimit(flag: boolean): void;
        getLowerLimit(): number;
        getUpperLimit(): number;
        setLimits(lower: number, upper: number): void;
        setMotorSpeed(speed: number): void;
        dump(log: (format: string, ...args: any[]) => void): void;
        private static draw_s_pA;
        private static draw_s_pB;
        private static draw_s_c1;
        private static draw_s_c2;
        private static draw_s_c3;
        private static draw_s_c4;
        private static draw_s_c5;
        private static draw_s_color_;
        private static draw_s_r;
        private static draw_s_rlo;
        private static draw_s_rhi;
        draw(draw: Draw): void;
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
        initialize(bA: Body, bB: Body, anchor1: XY, anchor2: XY): void;
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
        private static initVelocityConstraints_s_P;
        initVelocityConstraints(data: SolverData): void;
        private static solveVelocityConstraints_s_Cdot1;
        private static solveVelocityConstraints_s_impulse1;
        private static solveVelocityConstraints_s_impulse;
        private static solveVelocityConstraints_s_P;
        solveVelocityConstraints(data: SolverData): void;
        private static solvePositionConstraints_s_C1;
        private static solvePositionConstraints_s_P;
        private static solvePositionConstraints_s_impulse;
        solvePositionConstraints(data: SolverData): boolean;
        getAnchorA<T extends XY>(out: T): T;
        getAnchorB<T extends XY>(out: T): T;
        getReactionForce<T extends XY>(inv_dt: number, out: T): T;
        getReactionTorque(inv_dt: number): number;
        setStiffness(stiffness: number): void;
        getStiffness(): number;
        setDamping(damping: number): void;
        getDamping(): number;
        dump(log: (format: string, ...args: any[]) => void): void;
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
        initialize(bA: Body, bB: Body, anchor: Vec2, axis: Vec2): void;
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
        getMotorSpeed(): number;
        getMaxMotorTorque(): number;
        setSpringFrequencyHz(hz: number): void;
        getSpringFrequencyHz(): number;
        setSpringDampingRatio(ratio: number): void;
        getSpringDampingRatio(): number;
        private static initVelocityConstraints_s_d;
        private static initVelocityConstraints_s_P;
        initVelocityConstraints(data: SolverData): void;
        private static solveVelocityConstraints_s_P;
        solveVelocityConstraints(data: SolverData): void;
        private static solvePositionConstraints_s_d;
        private static solvePositionConstraints_s_P;
        solvePositionConstraints(data: SolverData): boolean;
        getDefinition(def: WheelJointDef): WheelJointDef;
        getAnchorA<T extends XY>(out: T): T;
        getAnchorB<T extends XY>(out: T): T;
        getReactionForce<T extends XY>(inv_dt: number, out: T): T;
        getReactionTorque(inv_dt: number): number;
        getJointTranslation(): number;
        getJointLinearSpeed(): number;
        getJointAngle(): number;
        getJointAngularSpeed(): number;
        getPrismaticJointTranslation(): number;
        getPrismaticJointSpeed(): number;
        getRevoluteJointAngle(): number;
        getRevoluteJointSpeed(): number;
        isMotorEnabled(): boolean;
        setEnableMotor(flag: boolean): void;
        setMotorSpeed(speed: number): void;
        setMaxMotorTorque(force: number): void;
        getMotorTorque(inv_dt: number): number;
        isLimitEnabled(): boolean;
        setEnableLimit(flag: boolean): void;
        getLowerLimit(): number;
        getUpperLimit(): number;
        setLimits(lower: number, upper: number): void;
        dump(log: (format: string, ...args: any[]) => void): void;
        private static draw_s_pA;
        private static draw_s_pB;
        private static draw_s_axis;
        private static draw_s_c1;
        private static draw_s_c2;
        private static draw_s_c3;
        private static draw_s_c4;
        private static draw_s_c5;
        private static draw_s_lower;
        private static draw_s_upper;
        private static draw_s_perp;
        draw(draw: Draw): void;
    }
}
declare namespace b2 {
    class World {
        readonly contactManager: ContactManager;
        bodyList: Body;
        jointList: Joint;
        particleSystemList: ParticleSystem;
        bodyCount: number;
        jointCount: number;
        readonly gravity: Vec2;
        allowSleep: boolean;
        destructionListener: DestructionListener;
        debugDrawInstance: Draw;
        inv_dt0: number;
        newContacts: boolean;
        locked: boolean;
        clearForcesFlag: boolean;
        warmStarting: boolean;
        continuousPhysics: boolean;
        subStepping: boolean;
        stepComplete: boolean;
        readonly profile: Profile;
        readonly island: Island;
        readonly s_stack: Array<Body>;
        controllerList: Controller;
        controllerCount: number;
        constructor(gravity: XY);
        setDestructionListener(listener: DestructionListener): void;
        setContactFilter(filter: ContactFilter): void;
        setContactListener(listener: ContactListener): void;
        setDebugDraw(debugDraw: Draw): void;
        addBody(b: Body): void;
        createBody(def?: IBodyDef): Body;
        destroyBody(b: Body): void;
        private static _createJoint;
        private static _destroyJoint;
        createJoint(def: IAreaJointDef): AreaJoint;
        createJoint(def: IDistanceJointDef): DistanceJoint;
        createJoint(def: IFrictionJointDef): FrictionJoint;
        createJoint(def: IGearJointDef): GearJoint;
        createJoint(def: IMotorJointDef): MotorJoint;
        createJoint(def: IMouseJointDef): MouseJoint;
        createJoint(def: IPrismaticJointDef): PrismaticJoint;
        createJoint(def: IPulleyJointDef): PulleyJoint;
        createJoint(def: IRevoluteJointDef): RevoluteJoint;
        createJoint(def: IWeldJointDef): WeldJoint;
        createJoint(def: IWheelJointDef): WheelJoint;
        addJoint(j: Joint): void;
        destroyJoint(j: Joint): void;
        createParticleSystem(def: ParticleSystemDef): ParticleSystem;
        destroyParticleSystem(p: ParticleSystem): void;
        calculateReasonableParticleIterations(timeStep: number): number;
        private static step_s_step;
        private static step_s_stepTimer;
        private static step_s_timer;
        step(dt: number, velocityIterations: number, positionIterations: number, particleIterations?: number): void;
        clearForces(): void;
        drawParticleSystem(system: ParticleSystem): void;
        private static debugdraw_s_color;
        private static debugdraw_s_vs;
        private static debugdraw_s_xf;
        debugDraw(): void;
        queryAABB(callback: QueryCallback, aabb: AABB): void;
        queryAABB(aabb: AABB, fn: QueryCallbackFunction): void;
        private _queryAABB;
        queryAllAABB(aabb: AABB, out?: Fixture[]): Fixture[];
        queryPointAABB(callback: QueryCallback, point: XY): void;
        queryPointAABB(point: XY, fn: QueryCallbackFunction): void;
        private _queryPointAABB;
        queryAllPointAABB(point: XY, out?: Fixture[]): Fixture[];
        queryFixtureShape(callback: QueryCallback, shape: Shape, index: number, transform: Transform): void;
        queryFixtureShape(shape: Shape, index: number, transform: Transform, fn: QueryCallbackFunction): void;
        private static queryFixtureShape_s_aabb;
        private _queryFixtureShape;
        queryAllFixtureShape(shape: Shape, index: number, transform: Transform, out?: Fixture[]): Fixture[];
        queryFixturePoint(callback: QueryCallback, point: XY): void;
        queryFixturePoint(point: XY, fn: QueryCallbackFunction): void;
        private _queryFixturePoint;
        queryAllFixturePoint(point: XY, out?: Fixture[]): Fixture[];
        rayCast(callback: RayCastCallback, point1: XY, point2: XY): void;
        rayCast(point1: XY, point2: XY, fn: RayCastCallbackFunction): void;
        private static rayCast_s_input;
        private static rayCast_s_output;
        private static rayCast_s_point;
        private _rayCast;
        rayCastOne(point1: XY, point2: XY): Fixture;
        rayCastAll(point1: XY, point2: XY, out?: Fixture[]): Fixture[];
        getParticleSystemList(): ParticleSystem;
        getContactList(): Contact;
        setAllowSleeping(flag: boolean): void;
        getAllowSleeping(): boolean;
        setWarmStarting(flag: boolean): void;
        getWarmStarting(): boolean;
        setContinuousPhysics(flag: boolean): void;
        getContinuousPhysics(): boolean;
        setSubStepping(flag: boolean): void;
        getSubStepping(): boolean;
        getProxyCount(): number;
        getBodyCount(): number;
        getJointCount(): number;
        getContactCount(): number;
        getTreeHeight(): number;
        getTreeBalance(): number;
        getTreeQuality(): number;
        setGravity(gravity: XY, wake?: boolean): void;
        isLocked(): boolean;
        setAutoClearForces(flag: boolean): void;
        getAutoClearForces(): boolean;
        shiftOrigin(newOrigin: XY): void;
        getContactManager(): ContactManager;
        getProfile(): Profile;
        dump(log: (format: string, ...args: any[]) => void): void;
        drawShape(fixture: Fixture, color: Color): void;
        solve(step: TimeStep): void;
        private static solveTOI_s_subStep;
        private static solveTOI_s_backup;
        private static solveTOI_s_backup1;
        private static solveTOI_s_backup2;
        private static solveTOI_s_toi_input;
        private static solveTOI_s_toi_output;
        solveTOI(step: TimeStep): void;
        addController(controller: Controller): Controller;
        removeController(controller: Controller): Controller;
    }
}
declare namespace b2 {
    enum ParticleGroupFlag {
        SolidParticleGroup = 1,
        RigidParticleGroup = 2,
        ParticleGroupCanBeEmpty = 4,
        ParticleGroupWillBeDestroyed = 8,
        ParticleGroupNeedsUpdateDepth = 16,
        ParticleGroupInternalMask = 24
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
        group?: ParticleGroup;
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
        group: ParticleGroup;
    }
    class ParticleGroup {
        readonly system: ParticleSystem;
        firstIndex: number;
        lastIndex: number;
        groupFlags: ParticleGroupFlag;
        strength: number;
        prev: ParticleGroup;
        next: ParticleGroup;
        timestamp: number;
        mass: number;
        inertia: number;
        readonly center: Vec2;
        readonly linearVelocity: Vec2;
        angularVelocity: number;
        readonly transform: Transform;
        userData: any;
        constructor(system: ParticleSystem);
        getNext(): ParticleGroup;
        getParticleSystem(): ParticleSystem;
        getParticleCount(): number;
        getBufferIndex(): number;
        containsParticle(index: number): boolean;
        getAllParticleFlags(): ParticleFlag;
        getGroupFlags(): ParticleGroupFlag;
        setGroupFlags(flags: number): void;
        getMass(): number;
        getInertia(): number;
        getCenter(): Vec2;
        getLinearVelocity(): Vec2;
        getAngularVelocity(): number;
        getTransform(): Transform;
        getPosition(): Vec2;
        getAngle(): number;
        getLinearVelocityFromWorldPoint<T extends XY>(worldPoint: XY, out: T): T;
        static readonly GetLinearVelocityFromWorldPoint_s_t0: Vec2;
        getUserData(): void;
        setUserData(data: any): void;
        applyForce(force: XY): void;
        applyLinearImpulse(impulse: XY): void;
        destroyParticles(callDestructionListener: boolean): void;
        updateStatistics(): void;
    }
}
declare namespace b2 {
    class GrowableBuffer<T> {
        data: T[];
        count: number;
        capacity: number;
        allocator: () => T;
        constructor(allocator: () => T);
        append(): number;
        reserve(newCapacity: number): void;
        grow(): void;
        free(): void;
        shorten(newEnd: number): void;
        getData(): T[];
        getCount(): number;
        setCount(newCount: number): void;
        getCapacity(): number;
        removeIf(pred: (t: T) => boolean): void;
        unique(pred: (a: T, b: T) => boolean): void;
    }
    type ParticleIndex = number;
    class FixtureParticleQueryCallback extends QueryCallback {
        system: ParticleSystem;
        constructor(system: ParticleSystem);
        shouldQueryParticleSystem(system: ParticleSystem): boolean;
        reportFixture(fixture: Fixture): boolean;
        reportParticle(system: ParticleSystem, index: number): boolean;
        reportFixtureAndParticle(fixture: Fixture, childIndex: number, index: number): void;
    }
    class ParticleContact {
        indexA: number;
        indexB: number;
        weight: number;
        normal: Vec2;
        flags: ParticleFlag;
        setIndices(a: number, b: number): void;
        setWeight(w: number): void;
        setNormal(n: Vec2): void;
        setFlags(f: ParticleFlag): void;
        getIndexA(): number;
        getIndexB(): number;
        getWeight(): number;
        getNormal(): Vec2;
        getFlags(): ParticleFlag;
        isEqual(rhs: ParticleContact): boolean;
        isNotEqual(rhs: ParticleContact): boolean;
        approximatelyEqual(rhs: ParticleContact): boolean;
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
        clone(): ParticleSystemDef;
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
        handleIndexBuffer: ParticleSysteUserOverridableBuffer<ParticleHandle>;
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
        groupBuffer: Array<ParticleGroup>;
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
        groupList: ParticleGroup;
        def: ParticleSystemDef;
        world: World;
        prev: ParticleSystem;
        next: ParticleSystem;
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
        drop(): void;
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
        createParticle(def: IParticleDef): number;
        /**
         * Retrieve a handle to the particle at the specified index.
         *
         * Please see #ParticleHandle for why you might want a handle.
         */
        getParticleHandleFromIndex(index: number): ParticleHandle;
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
        destroyParticle(index: number, callDestructionListener?: boolean): void;
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
        destroyOldestParticle(index: number, callDestructionListener?: boolean): void;
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
        destroyParticlesInShape(shape: Shape, xf: Transform, callDestructionListener?: boolean): number;
        static readonly DestroyParticlesInShape_s_aabb: AABB;
        /**
         * Create a particle group whose properties have been defined.
         *
         * No reference to the definition is retained.
         *
         * warning: This function is locked during callbacks.
         */
        createParticleGroup(groupDef: IParticleGroupDef): ParticleGroup;
        static readonly createParticleGroup_s_transform: Transform;
        /**
         * Join two particle groups.
         *
         * warning: This function is locked during callbacks.
         *
         * @param groupA the first group. Expands to encompass the second group.
         * @param groupB the second group. It is destroyed.
         */
        joinParticleGroups(groupA: ParticleGroup, groupB: ParticleGroup): void;
        /**
         * Split particle group into multiple disconnected groups.
         *
         * warning: This function is locked during callbacks.
         *
         * @param group the group to be split.
         */
        splitParticleGroup(group: ParticleGroup): void;
        /**
         * Get the world particle group list. With the returned group,
         * use ParticleGroup::GetNext to get the next group in the
         * world list.
         *
         * A null group indicates the end of the list.
         *
         * @return the head of the world particle group list.
         */
        getParticleGroupList(): ParticleGroup;
        /**
         * Get the number of particle groups.
         */
        getParticleGroupCount(): number;
        /**
         * Get the number of particles.
         */
        getParticleCount(): number;
        /**
         * Get the maximum number of particles.
         */
        getMaxParticleCount(): number;
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
        setMaxParticleCount(count: number): void;
        /**
         * Get all existing particle flags.
         */
        getAllParticleFlags(): ParticleFlag;
        /**
         * Get all existing particle group flags.
         */
        getAllGroupFlags(): ParticleGroupFlag;
        /**
         * Pause or unpause the particle system. When paused,
         * World::Step() skips over this particle system. All
         * ParticleSystem function calls still work.
         *
         * @param paused paused is true to pause, false to un-pause.
         */
        setPaused(paused: boolean): void;
        /**
         * Initially, true, then, the last value passed into
         * SetPaused().
         *
         * @return true if the particle system is being updated in World::Step().
         */
        getPaused(): boolean;
        /**
         * Change the particle density.
         *
         * Particle density affects the mass of the particles, which in
         * turn affects how the particles interact with Bodies. Note
         * that the density does not affect how the particles interact
         * with each other.
         */
        setDensity(density: number): void;
        /**
         * Get the particle density.
         */
        getDensity(): number;
        /**
         * Change the particle gravity scale. Adjusts the effect of the
         * global gravity vector on particles.
         */
        setGravityScale(gravityScale: number): void;
        /**
         * Get the particle gravity scale.
         */
        getGravityScale(): number;
        /**
         * Damping is used to reduce the velocity of particles. The
         * damping parameter can be larger than 1.0f but the damping
         * effect becomes sensitive to the time step when the damping
         * parameter is large.
         */
        setDamping(damping: number): void;
        /**
         * Get damping for particles
         */
        getDamping(): number;
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
        setStaticPressureIterations(iterations: number): void;
        /**
         * Get the number of iterations for static pressure of
         * particles.
         */
        getStaticPressureIterations(): number;
        /**
         * Change the particle radius.
         *
         * You should set this only once, on world start.
         * If you change the radius during execution, existing particles
         * may explode, shrink, or behave unexpectedly.
         */
        setRadius(radius: number): void;
        /**
         * Get the particle radius.
         */
        getRadius(): number;
        /**
         * Get the position of each particle
         *
         * Array is length GetParticleCount()
         *
         * @return the pointer to the head of the particle positions array.
         */
        getPositionBuffer(): Vec2[];
        /**
         * Get the velocity of each particle
         *
         * Array is length GetParticleCount()
         *
         * @return the pointer to the head of the particle velocities array.
         */
        getVelocityBuffer(): Vec2[];
        /**
         * Get the color of each particle
         *
         * Array is length GetParticleCount()
         *
         * @return the pointer to the head of the particle colors array.
         */
        getColorBuffer(): Color[];
        /**
         * Get the particle-group of each particle.
         *
         * Array is length GetParticleCount()
         *
         * @return the pointer to the head of the particle group array.
         */
        getGroupBuffer(): Array<ParticleGroup>;
        /**
         * Get the weight of each particle
         *
         * Array is length GetParticleCount()
         *
         * @return the pointer to the head of the particle positions array.
         */
        getWeightBuffer(): number[];
        /**
         * Get the user-specified data of each particle.
         *
         * Array is length GetParticleCount()
         *
         * @return the pointer to the head of the particle user-data array.
         */
        getUserDataBuffer<T>(): T[];
        /**
         * Get the flags for each particle. See the ParticleFlag enum.
         *
         * Array is length GetParticleCount()
         *
         * @return the pointer to the head of the particle-flags array.
         */
        getFlagsBuffer(): ParticleFlag[];
        /**
         * Set flags for a particle. See the ParticleFlag enum.
         */
        setParticleFlags(index: number, newFlags: ParticleFlag): void;
        /**
         * Get flags for a particle. See the ParticleFlag enum.
         */
        getParticleFlags(index: number): ParticleFlag;
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
        setFlagsBuffer(buffer: ParticleFlag[]): void;
        setPositionBuffer(buffer: Vec2[] | Float32Array): void;
        setVelocityBuffer(buffer: TypedVec2[] | Float32Array): void;
        setColorBuffer(buffer: Color[] | Float32Array): void;
        setUserDataBuffer<T>(buffer: T[]): void;
        /**
         * Get contacts between particles
         * Contact data can be used for many reasons, for example to
         * trigger rendering or audio effects.
         */
        getContacts(): ParticleContact[];
        getContactCount(): number;
        /**
         * Get contacts between particles and bodies
         *
         * Contact data can be used for many reasons, for example to
         * trigger rendering or audio effects.
         */
        getBodyContacts(): ParticleBodyContact[];
        getBodyContactCount(): number;
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
        getPairs(): ParticlePair[];
        getPairCount(): number;
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
        getTriads(): ParticleTriad[];
        getTriadCount(): number;
        /**
         * Set an optional threshold for the maximum number of
         * consecutive particle iterations that a particle may contact
         * multiple bodies before it is considered a candidate for being
         * "stuck". Setting to zero or less disables.
         */
        setStuckThreshold(steps: number): void;
        /**
         * Get potentially stuck particles from the last step; the user
         * must decide if they are stuck or not, and if so, delete or
         * move them
         */
        getStuckCandidates(): number[];
        /**
         * Get the number of stuck particle candidates from the last
         * step.
         */
        getStuckCandidateCount(): number;
        /**
         * Compute the kinetic energy that can be lost by damping force
         */
        computeCollisionEnergy(): number;
        static readonly computeCollisionEnergy_s_v: Vec2;
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
        setStrictContactCheck(enabled: boolean): void;
        /**
         * Get the status of the strict contact check.
         */
        getStrictContactCheck(): boolean;
        /**
         * Set the lifetime (in seconds) of a particle relative to the
         * current time.  A lifetime of less than or equal to 0.0f
         * results in the particle living forever until it's manually
         * destroyed by the application.
         */
        setParticleLifetime(index: number, lifetime: number): void;
        /**
         * Get the lifetime (in seconds) of a particle relative to the
         * current time.  A value > 0.0f is returned if the particle is
         * scheduled to be destroyed in the future, values <= 0.0f
         * indicate the particle has an infinite lifetime.
         */
        getParticleLifetime(index: number): number;
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
        setDestructionByAge(enable: boolean): void;
        /**
         * Get whether the oldest particle will be destroyed in
         * CreateParticle() when the maximum number of particles are
         * present in the system.
         */
        getDestructionByAge(): boolean;
        /**
         * Get the array of particle expiration times indexed by
         * particle index.
         *
         * GetParticleCount() items are in the returned array.
         */
        getExpirationTimeBuffer(): number[];
        /**
         * Convert a expiration time value in returned by
         * GetExpirationTimeBuffer() to a time in seconds relative to
         * the current simulation time.
         */
        expirationTimeToLifetime(expirationTime: number): number;
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
        getIndexByExpirationTimeBuffer(): number[];
        /**
         * Apply an impulse to one particle. This immediately modifies
         * the velocity. Similar to Body::ApplyLinearImpulse.
         *
         * @param index the particle that will be modified.
         * @param impulse impulse the world impulse vector, usually in N-seconds or kg-m/s.
         */
        particleApplyLinearImpulse(index: number, impulse: XY): void;
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
        applyLinearImpulse(firstIndex: number, lastIndex: number, impulse: XY): void;
        static isSignificantForce(force: XY): boolean;
        /**
         * Apply a force to the center of a particle.
         *
         * @param index the particle that will be modified.
         * @param force the world force vector, usually in Newtons (N).
         */
        particleApplyForce(index: number, force: XY): void;
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
        applyForce(firstIndex: number, lastIndex: number, force: XY): void;
        /**
         * Get the next particle-system in the world's particle-system
         * list.
         */
        getNext(): ParticleSystem;
        /**
         * Query the particle system for all particles that potentially
         * overlap the provided AABB.
         * QueryCallback::ShouldQueryParticleSystem is ignored.
         *
         * @param callback a user implemented callback class.
         * @param aabb the query box.
         */
        queryAABB(callback: QueryCallback, aabb: AABB): void;
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
        queryShapeAABB(callback: QueryCallback, shape: Shape, xf: Transform, childIndex?: number): void;
        static readonly queryShapeAABB_s_aabb: AABB;
        queryPointAABB(callback: QueryCallback, point: XY, slop?: number): void;
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
        rayCast(callback: RayCastCallback, point1: XY, point2: XY): void;
        static readonly rayCast_s_aabb: AABB;
        static readonly rayCast_s_p: Vec2;
        static readonly rayCast_s_v: Vec2;
        static readonly rayCast_s_n: Vec2;
        static readonly rayCast_s_point: Vec2;
        /**
         * Compute the axis-aligned bounding box for all particles
         * contained within this particle system.
         * @param aabb Returns the axis-aligned bounding box of the system.
         */
        computeAABB(aabb: AABB): void;
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
        freeBuffer<T>(b: T[], capacity: number): void;
        freeUserOverridableBuffer<T>(b: ParticleSysteUserOverridableBuffer<T>): void;
        /**
         * Reallocate a buffer
         */
        reallocateBuffer3<T>(oldBuffer: T[], oldCapacity: number, newCapacity: number): T[];
        /**
         * Reallocate a buffer
         */
        reallocateBuffer5<T>(buffer: T[], userSuppliedCapacity: number, oldCapacity: number, newCapacity: number, deferred: boolean): T[];
        /**
         * Reallocate a buffer
         */
        reallocateBuffer4<T>(buffer: ParticleSysteUserOverridableBuffer<any>, oldCapacity: number, newCapacity: number, deferred: boolean): T[];
        requestBuffer<T>(buffer: T[]): T[];
        /**
         * Reallocate the handle / index map and schedule the allocation
         * of a new pool for handle allocation.
         */
        reallocateHandleBuffers(newCapacity: number): void;
        reallocateInternalAllocatedBuffers(capacity: number): void;
        createParticleForGroup(groupDef: IParticleGroupDef, xf: Transform, p: XY): void;
        createParticlesStrokeShapeForGroup(shape: Shape, groupDef: IParticleGroupDef, xf: Transform): void;
        static readonly createParticlesStrokeShapeForGroup_s_edge: EdgeShape;
        static readonly createParticlesStrokeShapeForGroup_s_d: Vec2;
        static readonly createParticlesStrokeShapeForGroup_s_p: Vec2;
        createParticlesFillShapeForGroup(shape: Shape, groupDef: IParticleGroupDef, xf: Transform): void;
        static readonly createParticlesFillShapeForGroup_s_aabb: AABB;
        static readonly createParticlesFillShapeForGroup_s_p: Vec2;
        createParticlesWithShapeForGroup(shape: Shape, groupDef: IParticleGroupDef, xf: Transform): void;
        createParticlesWithShapesForGroup(shapes: Shape[], shapeCount: number, groupDef: IParticleGroupDef, xf: Transform): void;
        cloneParticle(oldIndex: number, group: ParticleGroup): number;
        destroyParticlesInGroup(group: ParticleGroup, callDestructionListener?: boolean): void;
        destroyParticleGroup(group: ParticleGroup): void;
        static particleCanBeConnected(flags: ParticleFlag, group: ParticleGroup): boolean;
        updatePairsAndTriads(firstIndex: number, lastIndex: number, filter: ParticleSysteConnectionFilter): void;
        private static updatePairsAndTriads_s_dab;
        private static updatePairsAndTriads_s_dbc;
        private static updatePairsAndTriads_s_dca;
        updatePairsAndTriadsWithReactiveParticles(): void;
        static comparePairIndices(a: ParticlePair, b: ParticlePair): boolean;
        static matchPairIndices(a: ParticlePair, b: ParticlePair): boolean;
        static compareTriadIndices(a: ParticleTriad, b: ParticleTriad): boolean;
        static matchTriadIndices(a: ParticleTriad, b: ParticleTriad): boolean;
        static initializeParticleLists(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[]): void;
        mergeParticleListsInContact(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[]): void;
        static mergeParticleLists(listA: ParticleSysteParticleListNode, listB: ParticleSysteParticleListNode): void;
        static findLongestParticleList(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[]): ParticleSysteParticleListNode;
        mergeZombieParticleListNodes(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[], survivingList: ParticleSysteParticleListNode): void;
        static mergeParticleListAndNode(list: ParticleSysteParticleListNode, node: ParticleSysteParticleListNode): void;
        createParticleGroupsFromParticleList(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[], survivingList: ParticleSysteParticleListNode): void;
        updatePairsAndTriadsWithParticleList(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[]): void;
        computeDepth(): void;
        getInsideBoundsEnumerator(aabb: AABB): ParticleSysteInsideBoundsEnumerator;
        updateAllParticleFlags(): void;
        updateAllGroupFlags(): void;
        addContact(a: number, b: number, contacts: GrowableBuffer<ParticleContact>): void;
        static readonly addContact_s_d: Vec2;
        findContacts_Reference(contacts: GrowableBuffer<ParticleContact>): void;
        findContacts(contacts: GrowableBuffer<ParticleContact>): void;
        updateProxies_Reference(proxies: GrowableBuffer<ParticleSysteProxy>): void;
        updateProxies(proxies: GrowableBuffer<ParticleSysteProxy>): void;
        sortProxies(proxies: GrowableBuffer<ParticleSysteProxy>): void;
        filterContacts(contacts: GrowableBuffer<ParticleContact>): void;
        notifyContactListenerPreContact(particlePairs: ParticlePairSet): void;
        notifyContactListenerPostContact(particlePairs: ParticlePairSet): void;
        static particleContactIsZombie(contact: ParticleContact): boolean;
        updateContacts(exceptZombie: boolean): void;
        notifyBodyContactListenerPreContact(fixtureSet: ParticleSysteFixtureParticleSet): void;
        notifyBodyContactListenerPostContact(fixtureSet: ParticleSysteFixtureParticleSet): void;
        updateBodyContacts(): void;
        static readonly updateBodyContacts_s_aabb: AABB;
        updateBodyContacts_callback: ParticleSysteUpdateBodyContactsCallback;
        solve(step: TimeStep): void;
        static readonly solve_s_subStep: TimeStep;
        solveCollision(step: TimeStep): void;
        static readonly solveCollision_s_aabb: AABB;
        solveCollision_callback: ParticleSysteSolveCollisionCallback;
        limitVelocity(step: TimeStep): void;
        solveGravity(step: TimeStep): void;
        static readonly SolveGravity_s_gravity: Vec2;
        solveBarrier(step: TimeStep): void;
        static readonly solveBarrier_s_aabb: AABB;
        static readonly solveBarrier_s_va: Vec2;
        static readonly solveBarrier_s_vb: Vec2;
        static readonly solveBarrier_s_pba: Vec2;
        static readonly solveBarrier_s_vba: Vec2;
        static readonly solveBarrier_s_vc: Vec2;
        static readonly solveBarrier_s_pca: Vec2;
        static readonly solveBarrier_s_vca: Vec2;
        static readonly solveBarrier_s_qba: Vec2;
        static readonly solveBarrier_s_qca: Vec2;
        static readonly solveBarrier_s_dv: Vec2;
        static readonly solveBarrier_s_f: Vec2;
        solveStaticPressure(step: TimeStep): void;
        computeWeight(): void;
        solvePressure(step: TimeStep): void;
        static readonly solvePressure_s_f: Vec2;
        solveDamping(step: TimeStep): void;
        static readonly solveDamping_s_v: Vec2;
        static readonly solveDamping_s_f: Vec2;
        solveRigidDamping(): void;
        static readonly solveRigidDamping_s_t0: Vec2;
        static readonly solveRigidDamping_s_t1: Vec2;
        static readonly solveRigidDamping_s_p: Vec2;
        static readonly solveRigidDamping_s_v: Vec2;
        solveExtraDamping(): void;
        static readonly solveExtraDamping_s_v: Vec2;
        static readonly solveExtraDamping_s_f: Vec2;
        solveWall(): void;
        solveRigid(step: TimeStep): void;
        static readonly solveRigid_s_position: Vec2;
        static readonly solveRigid_s_rotation: Rot;
        static readonly solveRigid_s_transform: Transform;
        static readonly solveRigid_s_velocityTransform: Transform;
        solveElastic(step: TimeStep): void;
        static readonly solveElastic_s_pa: Vec2;
        static readonly solveElastic_s_pb: Vec2;
        static readonly solveElastic_s_pc: Vec2;
        static readonly solveElastic_s_r: Rot;
        static readonly solveElastic_s_t0: Vec2;
        solveSpring(step: TimeStep): void;
        static readonly solveSpring_s_pa: Vec2;
        static readonly solveSpring_s_pb: Vec2;
        static readonly solveSpring_s_d: Vec2;
        static readonly solveSpring_s_f: Vec2;
        solveTensile(step: TimeStep): void;
        static readonly solveTensile_s_weightedNormal: Vec2;
        static readonly solveTensile_s_s: Vec2;
        static readonly solveTensile_s_f: Vec2;
        solveViscous(): void;
        static readonly solveViscous_s_v: Vec2;
        static readonly solveViscous_s_f: Vec2;
        solveRepulsive(step: TimeStep): void;
        static readonly solveRepulsive_s_f: Vec2;
        solvePowder(step: TimeStep): void;
        static readonly solvePowder_s_f: Vec2;
        solveSolid(step: TimeStep): void;
        static readonly SolveSolid_s_f: Vec2;
        solveForce(step: TimeStep): void;
        solveColorMixing(): void;
        solveZombie(): void;
        /**
         * Destroy all particles which have outlived their lifetimes set
         * by SetParticleLifetime().
         */
        solveLifetimes(step: TimeStep): void;
        rotateBuffer(start: number, mid: number, end: number): void;
        getCriticalVelocity(step: TimeStep): number;
        getCriticalVelocitySquared(step: TimeStep): number;
        getCriticalPressure(step: TimeStep): number;
        getParticleStride(): number;
        getParticleMass(): number;
        getParticleInvMass(): number;
        /**
         * Get the world's contact filter if any particles with the
         * contactFilterParticle flag are present in the system.
         */
        getFixtureContactFilter(): ContactFilter;
        /**
         * Get the world's contact filter if any particles with the
         * particleContactFilterParticle flag are present in the
         * system.
         */
        getParticleContactFilter(): ContactFilter;
        /**
         * Get the world's contact listener if any particles with the
         * fixtureContactListenerParticle flag are present in the
         * system.
         */
        getFixtureContactListener(): ContactListener;
        /**
         * Get the world's contact listener if any particles with the
         * particleContactListenerParticle flag are present in the
         * system.
         */
        getParticleContactListener(): ContactListener;
        setUserOverridableBuffer<T>(buffer: ParticleSysteUserOverridableBuffer<T>, data: T[]): void;
        setGroupFlags(group: ParticleGroup, newFlags: ParticleGroupFlag): void;
        static bodyContactCompare(lhs: ParticleBodyContact, rhs: ParticleBodyContact): boolean;
        removeSpuriousBodyContacts(): void;
        private static removeSpuriousBodyContacts_s_n;
        private static removeSpuriousBodyContacts_s_pos;
        private static removeSpuriousBodyContacts_s_normal;
        detectStuckParticle(particle: number): void;
        /**
         * Determine whether a particle index is valid.
         */
        validateParticleIndex(index: number): boolean;
        /**
         * Get the time elapsed in
         * ParticleSystemDef::lifetimeGranularity.
         */
        getQuantizedTimeElapsed(): number;
        /**
         * Convert a lifetime in seconds to an expiration time.
         */
        lifetimeToExpirationTime(lifetime: number): number;
        forceCanBeApplied(flags: ParticleFlag): boolean;
        prepareForceBuffer(): void;
        isRigidGroup(group: ParticleGroup): boolean;
        getLinearVelocity(group: ParticleGroup, particleIndex: number, point: Vec2, out: Vec2): Vec2;
        initDampingParameter(invMass: number[], invInertia: number[], tangentDistance: number[], mass: number, inertia: number, center: Vec2, point: Vec2, normal: Vec2): void;
        initDampingParameterWithRigidGroupOrParticle(invMass: number[], invInertia: number[], tangentDistance: number[], isRigidGroup: boolean, group: ParticleGroup, particleIndex: number, point: Vec2, normal: Vec2): void;
        computeDampingImpulse(invMassA: number, invInertiaA: number, tangentDistanceA: number, invMassB: number, invInertiaB: number, tangentDistanceB: number, normalVelocity: number): number;
        applyDamping(invMass: number, invInertia: number, tangentDistance: number, isRigidGroup: boolean, group: ParticleGroup, particleIndex: number, impulse: number, normal: Vec2): void;
    }
    class ParticleSysteUserOverridableBuffer<T> {
        _data: T[];
        data: T[];
        userSuppliedCapacity: number;
    }
    class ParticleSysteProxy {
        index: number;
        tag: number;
        static compareProxyProxy(a: ParticleSysteProxy, b: ParticleSysteProxy): boolean;
        static compareTagProxy(a: number, b: ParticleSysteProxy): boolean;
        static compareProxyTag(a: ParticleSysteProxy, b: number): boolean;
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
        getNext(): number;
    }
    class ParticleSysteParticleListNode {
        /**
         * The head of the list.
         */
        list: ParticleSysteParticleListNode;
        /**
         * The next node in the list.
         */
        next: ParticleSysteParticleListNode;
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
        allocate(itemSize: number, count: number): number;
        clear(): void;
        getCount(): number;
        invalidate(itemIndex: number): void;
        getValidBuffer(): boolean[];
        getBuffer(): T[];
        setCount(count: number): void;
    }
    class ParticleSysteFixtureParticle {
        first: Fixture;
        second: number;
        constructor(fixture: Fixture, particle: number);
    }
    class ParticleSysteFixtureParticleSet extends ParticleSysteFixedSetAllocator<ParticleSysteFixtureParticle> {
        initialize(bodyContactBuffer: GrowableBuffer<ParticleBodyContact>, flagsBuffer: ParticleSysteUserOverridableBuffer<ParticleFlag>): void;
        find(pair: ParticleSysteFixtureParticle): number;
    }
    class ParticleSysteParticlePair {
        first: number;
        second: number;
        constructor(particleA: number, particleB: number);
    }
    class ParticlePairSet extends ParticleSysteFixedSetAllocator<ParticleSysteParticlePair> {
        initialize(contactBuffer: GrowableBuffer<ParticleContact>, flagsBuffer: ParticleSysteUserOverridableBuffer<ParticleFlag>): void;
        find(pair: ParticleSysteParticlePair): number;
    }
    class ParticleSysteConnectionFilter {
        /**
         * Is the particle necessary for connection?
         * A pair or a triad should contain at least one 'necessary'
         * particle.
         */
        isNecessary(index: number): boolean;
        /**
         * An additional condition for creating a pair.
         */
        shouldCreatePair(a: number, b: number): boolean;
        /**
         * An additional condition for creating a triad.
         */
        shouldCreateTriad(a: number, b: number, c: number): boolean;
    }
    class ParticleSysteDestroyParticlesInShapeCallback extends QueryCallback {
        system: ParticleSystem;
        shape: Shape;
        xf: Transform;
        callDestructionListener: boolean;
        destroyed: number;
        constructor(system: ParticleSystem, shape: Shape, xf: Transform, callDestructionListener: boolean);
        reportFixture(fixture: Fixture): boolean;
        reportParticle(particleSystem: ParticleSystem, index: number): boolean;
        isDestroyed(): number;
    }
    class ParticleSysteJoinParticleGroupsFilter extends ParticleSysteConnectionFilter {
        threshold: number;
        constructor(threshold: number);
        /**
         * An additional condition for creating a pair.
         */
        shouldCreatePair(a: number, b: number): boolean;
        /**
         * An additional condition for creating a triad.
         */
        shouldCreateTriad(a: number, b: number, c: number): boolean;
    }
    class ParticleSysteCompositeShape extends Shape {
        constructor(shapes: Shape[], shapeCount?: number);
        shapes: Shape[];
        shapeCount: number;
        clone(): Shape;
        getChildCount(): number;
        /**
         * @see Shape::TestPoint
         */
        testPoint(xf: Transform, p: XY): boolean;
        /**
         * @see Shape::ComputeDistance
         */
        computeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number;
        /**
         * Implement Shape.
         */
        rayCast(output: RayCastOutput, input: RayCastInput, xf: Transform, childIndex: number): boolean;
        /**
         * @see Shape::ComputeAABB
         */
        computeAABB(aabb: AABB, xf: Transform, childIndex: number): void;
        /**
         * @see Shape::ComputeMass
         */
        computeMass(massData: MassData, density: number): void;
        setupDistanceProxy(proxy: DistanceProxy, index: number): void;
        computeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number;
        dump(log: (format: string, ...args: any[]) => void): void;
    }
    class ParticleSysteReactiveFilter extends ParticleSysteConnectionFilter {
        flagsBuffer: ParticleSysteUserOverridableBuffer<ParticleFlag>;
        constructor(flagsBuffer: ParticleSysteUserOverridableBuffer<ParticleFlag>);
        isNecessary(index: number): boolean;
    }
    class ParticleSysteUpdateBodyContactsCallback extends FixtureParticleQueryCallback {
        contactFilter: ContactFilter;
        constructor(system: ParticleSystem, contactFilter?: ContactFilter);
        shouldCollideFixtureParticle(fixture: Fixture, particleSystem: ParticleSystem, particleIndex: number): boolean;
        reportFixtureAndParticle(fixture: Fixture, childIndex: number, a: number): void;
        static readonly ReportFixtureAndParticle_s_n: Vec2;
        static readonly ReportFixtureAndParticle_s_rp: Vec2;
    }
    class ParticleSysteSolveCollisionCallback extends FixtureParticleQueryCallback {
        step: TimeStep;
        constructor(system: ParticleSystem, step: TimeStep);
        reportFixtureAndParticle(fixture: Fixture, childIndex: number, a: number): void;
        static readonly ReportFixtureAndParticle_s_p1: Vec2;
        static readonly ReportFixtureAndParticle_s_output: RayCastOutput;
        static readonly ReportFixtureAndParticle_s_input: RayCastInput;
        static readonly ReportFixtureAndParticle_s_p: Vec2;
        static readonly ReportFixtureAndParticle_s_v: Vec2;
        static readonly ReportFixtureAndParticle_s_f: Vec2;
        reportParticle(system: ParticleSystem, index: number): boolean;
    }
}
declare namespace b2 {
    class StackQueue<T> {
        readonly buffer: Array<T>;
        front: number;
        back: number;
        readonly capacity: number;
        constructor(capacity: number);
        push(item: T): void;
        pop(): void;
        empty(): boolean;
        getFront(): T;
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
        addGenerator(center: Vec2, tag: number, necessary: boolean): void;
        /**
         * Generate the Voronoi diagram. It is rasterized with a given
         * interval in the same range as the necessary generators exist.
         *
         * @param radius the interval of the diagram.
         * @param margin margin for which the range of the diagram is extended.
         */
        generate(radius: number, margin: number): void;
        /**
         * Enumerate all nodes that contain at least one necessary
         * generator.
         */
        getNodes(callback: VoronoiDiagraNodeCallback): void;
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
        copy(other: RopeTuning): this;
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
        create(def: RopeDef): void;
        setTuning(tuning: RopeTuning): void;
        step(dt: number, iterations: number, position: Vec2): void;
        reset(position: Vec2): void;
        draw(draw: Draw): void;
        private solveStretchPBD;
        private solveStretchXPBD;
        private solveBendPBDAngle;
        private solveBendXPBDAngle;
        private solveBendPBDDistance;
        private solveBendPBDHeight;
        private solveBendPBDTriangle;
        private applyBendForces;
    }
}
