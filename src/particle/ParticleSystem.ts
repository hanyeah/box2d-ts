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
  function std_iter_swap<T>(array: T[], a: number, b: number): void {
    const tmp: T = array[a];
    array[a] = array[b];
    array[b] = tmp;
  }

  function default_compare<T>(a: T, b: T): boolean { return a < b; }

  function std_sort<T>(array: T[], first: number = 0, len: number = array.length - first, cmp: (a: T, b: T) => boolean = default_compare): T[] {
    let left = first;
    const stack: number[] = [];
    let pos = 0;

    for (; ; ) { /* outer loop */
      for (; left + 1 < len; len++) { /* sort left to len-1 */
        const pivot = array[left + Math.floor(Math.random() * (len - left))]; /* pick random pivot */
        stack[pos++] = len; /* sort right part later */
        for (let right = left - 1; ; ) { /* inner loop: partitioning */
          while (cmp(array[++right], pivot)) {} /* look for greater element */
          while (cmp(pivot, array[--len])) {} /* look for smaller element */
          if (right >= len) {
            break;
          } /* partition point found? */
          std_iter_swap(array, right, len); /* the only swap */
        } /* partitioned, continue left part */
      }
      if (pos === 0) {
        break;
      } /* stack empty? */
      left = len; /* left to right is sorted */
      len = stack[--pos]; /* get next range to sort */
    }

    return array;
  }

  function std_stable_sort<T>(array: T[], first: number = 0, len: number = array.length - first, cmp: (a: T, b: T) => boolean = default_compare): T[] {
    return std_sort(array, first, len, cmp);
  }

  function std_remove_if<T>(array: T[], predicate: (value: T) => boolean, length: number = array.length) {
    let l = 0;

    for (let c = 0; c < length; ++c) {
      // if we can be collapsed, keep l where it is.
      if (predicate(array[c])) {
        continue;
      }

      // this node can't be collapsed; push it back as far as we can.
      if (c === l) {
        ++l;
        continue; // quick exit if we're already in the right spot
      }

      // array[l++] = array[c];
      std_iter_swap(array, l++, c);
    }

    return l;
  }

  function std_lower_bound<A, B>(array: A[], first: number, last: number, val: B, cmp: (a: A, b: B) => boolean): number {
    let count = last - first;
    while (count > 0) {
      const step = Math.floor(count / 2);
      let it = first + step;

      if (cmp(array[it], val)) {
        first = ++it;
        count -= step + 1;
      } else {
        count = step;
      }
    }
    return first;
  }

  function std_upper_bound<A, B>(array: B[], first: number, last: number, val: A, cmp: (a: A, b: B) => boolean): number {
    let count = last - first;
    while (count > 0) {
      const step = Math.floor(count / 2);
      let it = first + step;

      if (!cmp(val, array[it])) {
        first = ++it;
        count -= step + 1;
      } else {
        count = step;
      }
    }
    return first;
  }

  function std_rotate<T>(array: T[], first: number, n_first: number, last: number): void {
    let next = n_first;
    while (first !== next) {
      std_iter_swap(array, first++, next++);
      if (next === last) {
        next = n_first;
      } else if (first === n_first) {
        n_first = next;
      }
    }
  }

  function std_unique<T>(array: T[], first: number, last: number, cmp: (a: T, b: T) => boolean): number {
    if (first === last) {
      return last;
    }
    let result = first;
    while (++first !== last) {
      if (!cmp(array[result], array[first])) {
        ///array[++result] = array[first];
        std_iter_swap(array, ++result, first);
      }
    }
    return ++result;
  }

  export class GrowableBuffer<T> {
    public data: T[] = [];
    public count: number = 0;
    public capacity: number = 0;
    public allocator: () => T;

    constructor(allocator: () => T) {
      this.allocator = allocator;
    }

    public Append(): number {
      if (this.count >= this.capacity) {
        this.Grow();
      }
      return this.count++;
    }

    public Reserve(newCapacity: number): void {
      if (this.capacity >= newCapacity) {
        return;
      }

      // DEBUG: Assert(this.capacity === this.data.length);
      for (let i = this.capacity; i < newCapacity; ++i) {
        this.data[i] = this.allocator();
      }
      this.capacity = newCapacity;
    }

    public Grow(): void {
      // Double the capacity.
      const newCapacity = this.capacity ? 2 * this.capacity : minParticleSystemBufferCapacity;
      // DEBUG: Assert(newCapacity > this.capacity);
      this.Reserve(newCapacity);
    }

    public Free(): void {
      if (this.data.length === 0) {
        return;
      }

      this.data = [];
      this.capacity = 0;
      this.count = 0;
    }

    public Shorten(newEnd: number): void {
      // DEBUG: Assert(false);
    }

    public Data(): T[] {
      return this.data;
    }

    public GetCount(): number {
      return this.count;
    }

    public SetCount(newCount: number): void {
      // DEBUG: Assert(0 <= newCount && newCount <= this.capacity);
      this.count = newCount;
    }

    public GetCapacity(): number {
      return this.capacity;
    }

    public RemoveIf(pred: (t: T) => boolean): void {
      // DEBUG: let count = 0;
      // DEBUG: for (let i = 0; i < this.count; ++i) {
      // DEBUG:   if (!pred(this.data[i])) {
      // DEBUG:     count++;
      // DEBUG:   }
      // DEBUG: }

      this.count = std_remove_if(this.data, pred, this.count);

      // DEBUG: Assert(count === this.count);
    }

    public Unique(pred: (a: T, b: T) => boolean): void {
      this.count = std_unique(this.data, 0, this.count, pred);
    }
  }

  export type ParticleIndex = number;

  export class FixtureParticleQueryCallback extends QueryCallback {
    public system: ParticleSystem;
    constructor(system: ParticleSystem) {
      super();
      this.system = system;
    }
    public ShouldQueryParticleSystem(system: ParticleSystem): boolean {
      // Skip reporting particles.
      return false;
    }
    public ReportFixture(fixture: Fixture): boolean {
      if (fixture.IsSensor()) {
        return true;
      }
      const shape = fixture.GetShape();
      const childCount = shape.GetChildCount();
      for (let childIndex = 0; childIndex < childCount; childIndex++) {
        const aabb = fixture.GetAABB(childIndex);
        const enumerator = this.system.GetInsideBoundsEnumerator(aabb);
        let index: number;
        while ((index = enumerator.GetNext()) >= 0) {
          this.ReportFixtureAndParticle(fixture, childIndex, index);
        }
      }
      return true;
    }
    public ReportParticle(system: ParticleSystem, index: number): boolean {
      return false;
    }
    public ReportFixtureAndParticle(fixture: Fixture, childIndex: number, index: number): void {
      // DEBUG: Assert(false); // pure virtual
    }
  }

  export class ParticleContact {
    public indexA: number = 0;
    public indexB: number = 0;
    public weight: number = 0;
    public normal: Vec2 = new Vec2();
    public flags: ParticleFlag = 0;

    public SetIndices(a: number, b: number): void {
      // DEBUG: Assert(a <= maxParticleIndex && b <= maxParticleIndex);
      this.indexA = a;
      this.indexB = b;
    }

    public SetWeight(w: number): void {
      this.weight = w;
    }

    public SetNormal(n: Vec2): void {
      this.normal.Copy(n);
    }

    public SetFlags(f: ParticleFlag): void {
      this.flags = f;
    }

    public GetIndexA(): number {
      return this.indexA;
    }

    public GetIndexB(): number {
      return this.indexB;
    }

    public GetWeight(): number {
      return this.weight;
    }

    public GetNormal(): Vec2 {
      return this.normal;
    }

    public GetFlags(): ParticleFlag {
      return this.flags;
    }

    public IsEqual(rhs: ParticleContact): boolean {
      return this.indexA === rhs.indexA && this.indexB === rhs.indexB && this.flags === rhs.flags && this.weight === rhs.weight && this.normal.x === rhs.normal.x && this.normal.y === rhs.normal.y;
    }

    public IsNotEqual(rhs: ParticleContact): boolean {
      return !this.IsEqual(rhs);
    }

    public ApproximatelyEqual(rhs: ParticleContact): boolean {
      const MAX_WEIGHT_DIFF = 0.01; // Weight 0 ~ 1, so about 1%
      const MAX_NORMAL_DIFF_SQ = 0.01 * 0.01; // Normal length = 1, so 1%
      return this.indexA === rhs.indexA && this.indexB === rhs.indexB && this.flags === rhs.flags && Abs(this.weight - rhs.weight) < MAX_WEIGHT_DIFF && Vec2.DistanceSquaredVV(this.normal, rhs.normal) < MAX_NORMAL_DIFF_SQ;
    }
  }

  export class ParticleBodyContact {
    public index: number = 0; // Index of the particle making contact.
    public body!: Body; // The body making contact.
    public fixture!: Fixture; // The specific fixture making contact
    public weight: number = 0.0; // Weight of the contact. A value between 0.0f and 1.0f.
    public normal: Vec2 = new Vec2(); // The normalized direction from the particle to the body.
    public mass: number = 0.0; // The effective mass used in calculating force.
  }

  export class ParticlePair {
    public indexA: number = 0; // Indices of the respective particles making pair.
    public indexB: number = 0;
    public flags: ParticleFlag = 0; // The logical sum of the particle flags. See the ParticleFlag enum.
    public strength: number = 0.0; // The strength of cohesion among the particles.
    public distance: number = 0.0; // The initial distance of the particles.
  }

  export class ParticleTriad {
    public indexA: number = 0; // Indices of the respective particles making triad.
    public indexB: number = 0;
    public indexC: number = 0;
    public flags: ParticleFlag = 0; // The logical sum of the particle flags. See the ParticleFlag enum.
    public strength: number = 0.0; // The strength of cohesion among the particles.
    public pa: Vec2 = new Vec2(0.0, 0.0); // Values used for calculation.
    public pb: Vec2 = new Vec2(0.0, 0.0);
    public pc: Vec2 = new Vec2(0.0, 0.0);
    public ka: number = 0.0;
    public kb: number = 0.0;
    public kc: number = 0.0;
    public s: number = 0.0;
  }

  export class ParticleSystemDef {
    // Initialize physical coefficients to the maximum values that
    // maintain numerical stability.

    /**
     * Enable strict Particle/Body contact check.
     * See SetStrictContactCheck for details.
     */
    public strictContactCheck: boolean = false;

    /**
     * Set the particle density.
     * See SetDensity for details.
     */
    public density: number = 1.0;

    /**
     * Change the particle gravity scale. Adjusts the effect of the
     * global gravity vector on particles. Default value is 1.0f.
     */
    public gravityScale: number = 1.0;

    /**
     * Particles behave as circles with this radius. In Box2D units.
     */
    public radius: number = 1.0;

    /**
     * Set the maximum number of particles.
     * By default, there is no maximum. The particle buffers can
     * continue to grow while World's block allocator still has
     * memory.
     * See SetMaxParticleCount for details.
     */
    public maxCount: number = 0;

    /**
     * Increases pressure in response to compression
     * Smaller values allow more compression
     */
    public pressureStrength: number = 0.005;

    /**
     * Reduces velocity along the collision normal
     * Smaller value reduces less
     */
    public dampingStrength: number = 1.0;

    /**
     * Restores shape of elastic particle groups
     * Larger values increase elastic particle velocity
     */
    public elasticStrength: number = 0.25;

    /**
     * Restores length of spring particle groups
     * Larger values increase spring particle velocity
     */
    public springStrength: number = 0.25;

    /**
     * Reduces relative velocity of viscous particles
     * Larger values slow down viscous particles more
     */
    public viscousStrength: number = 0.25;

    /**
     * Produces pressure on tensile particles
     * 0~0.2. Larger values increase the amount of surface tension.
     */
    public surfaceTensionPressureStrength: number = 0.2;

    /**
     * Smoothes outline of tensile particles
     * 0~0.2. Larger values result in rounder, smoother,
     * water-drop-like clusters of particles.
     */
    public surfaceTensionNormalStrength: number = 0.2;

    /**
     * Produces additional pressure on repulsive particles
     * Larger values repulse more
     * Negative values mean attraction. The range where particles
     * behave stably is about -0.2 to 2.0.
     */
    public repulsiveStrength: number = 1.0;

    /**
     * Produces repulsion between powder particles
     * Larger values repulse more
     */
    public powderStrength: number = 0.5;

    /**
     * Pushes particles out of solid particle group
     * Larger values repulse more
     */
    public ejectionStrength: number = 0.5;

    /**
     * Produces static pressure
     * Larger values increase the pressure on neighboring partilces
     * For a description of static pressure, see
     * http://en.wikipedia.org/wiki/Static_pressure#Static_pressure_in_fluid_dynamics
     */
    public staticPressureStrength: number = 0.2;

    /**
     * Reduces instability in static pressure calculation
     * Larger values make stabilize static pressure with fewer
     * iterations
     */
    public staticPressureRelaxation: number = 0.2;

    /**
     * Computes static pressure more precisely
     * See SetStaticPressureIterations for details
     */
    public staticPressureIterations: number = 8;

    /**
     * Determines how fast colors are mixed
     * 1.0f ==> mixed immediately
     * 0.5f ==> mixed half way each simulation step (see
     * World::Step())
     */
    public colorMixingStrength: number = 0.5;

    /**
     * Whether to destroy particles by age when no more particles
     * can be created.  See #ParticleSystem::SetDestructionByAge()
     * for more information.
     */
    public destroyByAge: boolean = true;

    /**
     * Granularity of particle lifetimes in seconds.  By default
     * this is set to (1.0f / 60.0f) seconds.  ParticleSystem uses
     * a 32-bit signed value to track particle lifetimes so the
     * maximum lifetime of a particle is (2^32 - 1) / (1.0f /
     * lifetimeGranularity) seconds. With the value set to 1/60 the
     * maximum lifetime or age of a particle is 2.27 years.
     */
    public lifetimeGranularity: number = 1.0 / 60.0;

    public Copy(def: ParticleSystemDef): ParticleSystemDef {
      this.strictContactCheck = def.strictContactCheck;
      this.density = def.density;
      this.gravityScale = def.gravityScale;
      this.radius = def.radius;
      this.maxCount = def.maxCount;
      this.pressureStrength = def.pressureStrength;
      this.dampingStrength = def.dampingStrength;
      this.elasticStrength = def.elasticStrength;
      this.springStrength = def.springStrength;
      this.viscousStrength = def.viscousStrength;
      this.surfaceTensionPressureStrength = def.surfaceTensionPressureStrength;
      this.surfaceTensionNormalStrength = def.surfaceTensionNormalStrength;
      this.repulsiveStrength = def.repulsiveStrength;
      this.powderStrength = def.powderStrength;
      this.ejectionStrength = def.ejectionStrength;
      this.staticPressureStrength = def.staticPressureStrength;
      this.staticPressureRelaxation = def.staticPressureRelaxation;
      this.staticPressureIterations = def.staticPressureIterations;
      this.colorMixingStrength = def.colorMixingStrength;
      this.destroyByAge = def.destroyByAge;
      this.lifetimeGranularity = def.lifetimeGranularity;
      return this;
    }

    public Clone(): ParticleSystemDef {
      return new ParticleSystemDef().Copy(this);
    }
  }

  export class ParticleSystem {
    public paused: boolean = false;
    public timestamp: number = 0;
    public allParticleFlags: ParticleFlag = 0;
    public needsUpdateAllParticleFlags: boolean = false;
    public allGroupFlags: ParticleGroupFlag = 0;
    public needsUpdateAllGroupFlags: boolean = false;
    public hasForce: boolean = false;
    public iterationIndex: number = 0;
    public inverseDensity: number = 0.0;
    public particleDiameter: number = 0.0;
    public inverseDiameter: number = 0.0;
    public squaredDiameter: number = 0.0;
    public count: number = 0;
    public internalAllocatedCapacity: number = 0;
    /**
     * Allocator for ParticleHandle instances.
     */
    ///handleAllocator: any = null;
    /**
     * Maps particle indicies to handles.
     */
    public handleIndexBuffer: ParticleSysteUserOverridableBuffer<ParticleHandle | null> = new ParticleSysteUserOverridableBuffer<ParticleHandle | null>();
    public flagsBuffer: ParticleSysteUserOverridableBuffer<ParticleFlag> = new ParticleSysteUserOverridableBuffer<ParticleFlag>();
    public positionBuffer: ParticleSysteUserOverridableBuffer<Vec2> = new ParticleSysteUserOverridableBuffer<Vec2>();
    public velocityBuffer: ParticleSysteUserOverridableBuffer<Vec2> = new ParticleSysteUserOverridableBuffer<Vec2>();
    public forceBuffer: Vec2[] = [];
    /**
     * this.weightBuffer is populated in ComputeWeight and used in
     * ComputeDepth(), SolveStaticPressure() and SolvePressure().
     */
    public weightBuffer: number[] = [];
    /**
     * When any particles have the flag staticPressureParticle,
     * this.staticPressureBuffer is first allocated and used in
     * SolveStaticPressure() and SolvePressure().  It will be
     * reallocated on subsequent CreateParticle() calls.
     */
    public staticPressureBuffer: number[] = [];
    /**
     * this.accumulationBuffer is used in many functions as a temporary
     * buffer for scalar values.
     */
    public accumulationBuffer: number[] = [];
    /**
     * When any particles have the flag tensileParticle,
     * this.accumulation2Buffer is first allocated and used in
     * SolveTensile() as a temporary buffer for vector values.  It
     * will be reallocated on subsequent CreateParticle() calls.
     */
    public accumulation2Buffer: Vec2[] = [];
    /**
     * When any particle groups have the flag solidParticleGroup,
     * this.depthBuffer is first allocated and populated in
     * ComputeDepth() and used in SolveSolid(). It will be
     * reallocated on subsequent CreateParticle() calls.
     */
    public depthBuffer: number[] = [];
    public colorBuffer: ParticleSysteUserOverridableBuffer<Color> = new ParticleSysteUserOverridableBuffer<Color>();
    public groupBuffer: Array<ParticleGroup | null> = [];
    public userDataBuffer: ParticleSysteUserOverridableBuffer<any> = new ParticleSysteUserOverridableBuffer();
    /**
     * Stuck particle detection parameters and record keeping
     */
    public stuckThreshold: number = 0;
    public lastBodyContactStepBuffer: ParticleSysteUserOverridableBuffer<number> = new ParticleSysteUserOverridableBuffer<number>();
    public bodyContactCountBuffer: ParticleSysteUserOverridableBuffer<number> = new ParticleSysteUserOverridableBuffer<number>();
    public consecutiveContactStepsBuffer: ParticleSysteUserOverridableBuffer<number> = new ParticleSysteUserOverridableBuffer<number>();
    public stuckParticleBuffer: GrowableBuffer<number> = new GrowableBuffer<number>(() => 0);
    public proxyBuffer: GrowableBuffer<ParticleSysteProxy> = new GrowableBuffer<ParticleSysteProxy>(() => new ParticleSysteProxy());
    public contactBuffer: GrowableBuffer<ParticleContact> = new GrowableBuffer<ParticleContact>(() => new ParticleContact());
    public bodyContactBuffer: GrowableBuffer<ParticleBodyContact> = new GrowableBuffer<ParticleBodyContact>(() => new ParticleBodyContact());
    public pairBuffer: GrowableBuffer<ParticlePair> = new GrowableBuffer<ParticlePair>(() => new ParticlePair());
    public triadBuffer: GrowableBuffer<ParticleTriad> = new GrowableBuffer<ParticleTriad>(() => new ParticleTriad());
    /**
     * Time each particle should be destroyed relative to the last
     * time this.timeElapsed was initialized.  Each unit of time
     * corresponds to ParticleSystemDef::lifetimeGranularity
     * seconds.
     */
    public expirationTimeBuffer: ParticleSysteUserOverridableBuffer<number> = new ParticleSysteUserOverridableBuffer<number>();
    /**
     * List of particle indices sorted by expiration time.
     */
    public indexByExpirationTimeBuffer: ParticleSysteUserOverridableBuffer<number> = new ParticleSysteUserOverridableBuffer<number>();
    /**
     * Time elapsed in 32:32 fixed point.  Each non-fractional unit
     * of time corresponds to
     * ParticleSystemDef::lifetimeGranularity seconds.
     */
    public timeElapsed: number = 0;
    /**
     * Whether the expiration time buffer has been modified and
     * needs to be resorted.
     */
    public expirationTimeBufferRequiresSorting: boolean = false;
    public groupCount: number = 0;
    public groupList: ParticleGroup | null = null;
    public def: ParticleSystemDef = new ParticleSystemDef();
    public world: World;
    public prev: ParticleSystem | null = null;
    public next: ParticleSystem | null = null;

    public static readonly xTruncBits: number = 12;
    public static readonly yTruncBits: number = 12;
    public static readonly tagBits: number = 8 * 4; // 8u * sizeof(uint32);
    public static readonly yOffset: number = 1 << (ParticleSystem.yTruncBits - 1);
    public static readonly yShift: number = ParticleSystem.tagBits - ParticleSystem.yTruncBits;
    public static readonly xShift: number = ParticleSystem.tagBits - ParticleSystem.yTruncBits - ParticleSystem.xTruncBits;
    public static readonly xScale: number = 1 << ParticleSystem.xShift;
    public static readonly xOffset: number = ParticleSystem.xScale * (1 << (ParticleSystem.xTruncBits - 1));
    public static readonly yMask: number = ((1 << ParticleSystem.yTruncBits) - 1) << ParticleSystem.yShift;
    public static readonly xMask: number = ~ParticleSystem.yMask;

    public static computeTag(x: number, y: number): number {
      ///return ((uint32)(y + yOffset) << yShift) + (uint32)(xScale * x + xOffset);
      return ((((y + ParticleSystem.yOffset) >>> 0) << ParticleSystem.yShift) + ((ParticleSystem.xScale * x + ParticleSystem.xOffset) >>> 0)) >>> 0;
    }

    public static computeRelativeTag(tag: number, x: number, y: number): number {
      ///return tag + (y << yShift) + (x << xShift);
      return (tag + (y << ParticleSystem.yShift) + (x << ParticleSystem.xShift)) >>> 0;
    }

    constructor(def: ParticleSystemDef, world: World) {
      this.SetStrictContactCheck(def.strictContactCheck);
      this.SetDensity(def.density);
      this.SetGravityScale(def.gravityScale);
      this.SetRadius(def.radius);
      this.SetMaxParticleCount(def.maxCount);
      // DEBUG: Assert(def.lifetimeGranularity > 0.0);
      this.def = def.Clone();
      this.world = world;
      this.SetDestructionByAge(this.def.destroyByAge);
    }

    public Drop(): void {
      while (this.groupList) {
        this.DestroyParticleGroup(this.groupList);
      }

      this.FreeUserOverridableBuffer(this.handleIndexBuffer);
      this.FreeUserOverridableBuffer(this.flagsBuffer);
      this.FreeUserOverridableBuffer(this.lastBodyContactStepBuffer);
      this.FreeUserOverridableBuffer(this.bodyContactCountBuffer);
      this.FreeUserOverridableBuffer(this.consecutiveContactStepsBuffer);
      this.FreeUserOverridableBuffer(this.positionBuffer);
      this.FreeUserOverridableBuffer(this.velocityBuffer);
      this.FreeUserOverridableBuffer(this.colorBuffer);
      this.FreeUserOverridableBuffer(this.userDataBuffer);
      this.FreeUserOverridableBuffer(this.expirationTimeBuffer);
      this.FreeUserOverridableBuffer(this.indexByExpirationTimeBuffer);
      this.FreeBuffer(this.forceBuffer, this.internalAllocatedCapacity);
      this.FreeBuffer(this.weightBuffer, this.internalAllocatedCapacity);
      this.FreeBuffer(this.staticPressureBuffer, this.internalAllocatedCapacity);
      this.FreeBuffer(this.accumulationBuffer, this.internalAllocatedCapacity);
      this.FreeBuffer(this.accumulation2Buffer, this.internalAllocatedCapacity);
      this.FreeBuffer(this.depthBuffer, this.internalAllocatedCapacity);
      this.FreeBuffer(this.groupBuffer, this.internalAllocatedCapacity);
    }

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
    public CreateParticle(def: IParticleDef): number {
      if (this.world.IsLocked()) { throw new Error(); }

      if (this.count >= this.internalAllocatedCapacity) {
        // Double the particle capacity.
        const capacity = this.count ? 2 * this.count : minParticleSystemBufferCapacity;
        this.ReallocateInternalAllocatedBuffers(capacity);
      }
      if (this.count >= this.internalAllocatedCapacity) {
        // If the oldest particle should be destroyed...
        if (this.def.destroyByAge) {
          this.DestroyOldestParticle(0, false);
          // Need to destroy this particle *now* so that it's possible to
          // create a new particle.
          this.SolveZombie();
        } else {
          return invalidParticleIndex;
        }
      }
      const index = this.count++;
      this.flagsBuffer.data[index] = 0;
      if (this.lastBodyContactStepBuffer.data) {
        this.lastBodyContactStepBuffer.data[index] = 0;
      }
      if (this.bodyContactCountBuffer.data) {
        this.bodyContactCountBuffer.data[index] = 0;
      }
      if (this.consecutiveContactStepsBuffer.data) {
        this.consecutiveContactStepsBuffer.data[index] = 0;
      }
      this.positionBuffer.data[index] = (this.positionBuffer.data[index] || new Vec2()).Copy(Maybe(def.position, Vec2.ZERO));
      this.velocityBuffer.data[index] = (this.velocityBuffer.data[index] || new Vec2()).Copy(Maybe(def.velocity, Vec2.ZERO));
      this.weightBuffer[index] = 0;
      this.forceBuffer[index] = (this.forceBuffer[index] || new Vec2()).SetZero();
      if (this.staticPressureBuffer) {
        this.staticPressureBuffer[index] = 0;
      }
      if (this.depthBuffer) {
        this.depthBuffer[index] = 0;
      }
      const color: Color = new Color().Copy(Maybe(def.color, Color.ZERO));
      if (this.colorBuffer.data || !color.IsZero()) {
        this.colorBuffer.data = this.RequestBuffer(this.colorBuffer.data);
        this.colorBuffer.data[index] = (this.colorBuffer.data[index] || new Color()).Copy(color);
      }
      if (this.userDataBuffer.data || def.userData) {
        this.userDataBuffer.data = this.RequestBuffer(this.userDataBuffer.data);
        this.userDataBuffer.data[index] = def.userData;
      }
      if (this.handleIndexBuffer.data) {
        this.handleIndexBuffer.data[index] = null;
      }
      ///Proxy& proxy = proxyBuffer.Append();
      const proxy = this.proxyBuffer.data[this.proxyBuffer.Append()];

      // If particle lifetimes are enabled or the lifetime is set in the particle
      // definition, initialize the lifetime.
      const lifetime = Maybe(def.lifetime, 0.0);
      const finiteLifetime = lifetime > 0.0;
      if (this.expirationTimeBuffer.data || finiteLifetime) {
        this.SetParticleLifetime(index, finiteLifetime ? lifetime :
          this.ExpirationTimeToLifetime(-this.GetQuantizedTimeElapsed()));
        // Add a reference to the newly added particle to the end of the
        // queue.
        this.indexByExpirationTimeBuffer.data[index] = index;
      }

      proxy.index = index;
      const group = Maybe(def.group, null);
      this.groupBuffer[index] = group;
      if (group) {
        if (group.firstIndex < group.lastIndex) {
          // Move particles in the group just before the new particle.
          this.RotateBuffer(group.firstIndex, group.lastIndex, index);
          // DEBUG: Assert(group.lastIndex === index);
          // Update the index range of the group to contain the new particle.
          group.lastIndex = index + 1;
        } else {
          // If the group is empty, reset the index range to contain only the
          // new particle.
          group.firstIndex = index;
          group.lastIndex = index + 1;
        }
      }
      this.SetParticleFlags(index, Maybe(def.flags, 0));
      return index;
    }

    /**
     * Retrieve a handle to the particle at the specified index.
     *
     * Please see #ParticleHandle for why you might want a handle.
     */
    public GetParticleHandleFromIndex(index: number): ParticleHandle {
      // DEBUG: Assert(index >= 0 && index < this.GetParticleCount() && index !== invalidParticleIndex);
      this.handleIndexBuffer.data = this.RequestBuffer(this.handleIndexBuffer.data);
      let handle = this.handleIndexBuffer.data[index];
      if (handle) {
        return handle;
      }
      // Create a handle.
      ///handle = handleAllocator.Allocate();
      handle = new ParticleHandle();
      // DEBUG: Assert(handle !== null);
      handle.SetIndex(index);
      this.handleIndexBuffer.data[index] = handle;
      return handle;
    }

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
    public DestroyParticle(index: number, callDestructionListener: boolean = false): void {
      let flags = ParticleFlag.zombieParticle;
      if (callDestructionListener) {
        flags |= ParticleFlag.destructionListenerParticle;
      }
      this.SetParticleFlags(index, this.flagsBuffer.data[index] | flags);
    }

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
    public DestroyOldestParticle(index: number, callDestructionListener: boolean = false): void {
      const particleCount = this.GetParticleCount();
      // DEBUG: Assert(index >= 0 && index < particleCount);
      // Make sure particle lifetime tracking is enabled.
      // DEBUG: Assert(this.indexByExpirationTimeBuffer.data !== null);
      // Destroy the oldest particle (preferring to destroy finite
      // lifetime particles first) to free a slot in the buffer.
      const oldestFiniteLifetimeParticle =
        this.indexByExpirationTimeBuffer.data[particleCount - (index + 1)];
      const oldestInfiniteLifetimeParticle =
        this.indexByExpirationTimeBuffer.data[index];
      this.DestroyParticle(
        this.expirationTimeBuffer.data[oldestFiniteLifetimeParticle] > 0.0 ?
          oldestFiniteLifetimeParticle : oldestInfiniteLifetimeParticle,
        callDestructionListener);
    }

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
    public DestroyParticlesInShape(shape: Shape, xf: Transform, callDestructionListener: boolean = false): number {
      const s_aabb = ParticleSystem.DestroyParticlesInShape_s_aabb;
      if (this.world.IsLocked()) { throw new Error(); }

      const callback = new ParticleSysteDestroyParticlesInShapeCallback(this, shape, xf, callDestructionListener);

      const aabb = s_aabb;
      shape.ComputeAABB(aabb, xf, 0);
      this.world.QueryAABB(callback, aabb);
      return callback.Destroyed();
    }
    public static readonly DestroyParticlesInShape_s_aabb = new AABB();

    /**
     * Create a particle group whose properties have been defined.
     *
     * No reference to the definition is retained.
     *
     * warning: This function is locked during callbacks.
     */
    public CreateParticleGroup(groupDef: IParticleGroupDef): ParticleGroup {
      const s_transform = ParticleSystem.CreateParticleGroup_s_transform;

      if (this.world.IsLocked()) { throw new Error(); }

      const transform = s_transform;
      transform.SetPositionAngle(Maybe(groupDef.position, Vec2.ZERO), Maybe(groupDef.angle, 0));
      const firstIndex = this.count;
      if (groupDef.shape) {
        this.CreateParticlesWithShapeForGroup(groupDef.shape, groupDef, transform);
      }
      if (groupDef.shapes) {
        this.CreateParticlesWithShapesForGroup(groupDef.shapes, Maybe(groupDef.shapeCount, groupDef.shapes.length), groupDef, transform);
      }
      if (groupDef.positionData) {
        const count = Maybe(groupDef.particleCount, groupDef.positionData.length);
        for (let i = 0; i < count; i++) {
          const p = groupDef.positionData[i];
          this.CreateParticleForGroup(groupDef, transform, p);
        }
      }
      const lastIndex = this.count;

      let group = new ParticleGroup(this);
      group.firstIndex = firstIndex;
      group.lastIndex = lastIndex;
      group.strength = Maybe(groupDef.strength, 1);
      group.userData = groupDef.userData;
      group.transform.Copy(transform);
      group.prev = null;
      group.next = this.groupList;
      if (this.groupList) {
        this.groupList.prev = group;
      }
      this.groupList = group;
      ++this.groupCount;
      for (let i = firstIndex; i < lastIndex; i++) {
        this.groupBuffer[i] = group;
      }
      this.SetGroupFlags(group, Maybe(groupDef.groupFlags, 0));

      // Create pairs and triads between particles in the group.
      const filter = new ParticleSysteConnectionFilter();
      this.UpdateContacts(true);
      this.UpdatePairsAndTriads(firstIndex, lastIndex, filter);

      if (groupDef.group) {
        this.JoinParticleGroups(groupDef.group, group);
        group = groupDef.group;
      }

      return group;
    }
    public static readonly CreateParticleGroup_s_transform = new Transform();

    /**
     * Join two particle groups.
     *
     * warning: This function is locked during callbacks.
     *
     * @param groupA the first group. Expands to encompass the second group.
     * @param groupB the second group. It is destroyed.
     */
    public JoinParticleGroups(groupA: ParticleGroup, groupB: ParticleGroup): void {
      if (this.world.IsLocked()) { throw new Error(); }

      // DEBUG: Assert(groupA !== groupB);
      this.RotateBuffer(groupB.firstIndex, groupB.lastIndex, this.count);
      // DEBUG: Assert(groupB.lastIndex === this.count);
      this.RotateBuffer(groupA.firstIndex, groupA.lastIndex, groupB.firstIndex);
      // DEBUG: Assert(groupA.lastIndex === groupB.firstIndex);

      // Create pairs and triads connecting groupA and groupB.
      const filter = new ParticleSysteJoinParticleGroupsFilter(groupB.firstIndex);
      this.UpdateContacts(true);
      this.UpdatePairsAndTriads(groupA.firstIndex, groupB.lastIndex, filter);

      for (let i = groupB.firstIndex; i < groupB.lastIndex; i++) {
        this.groupBuffer[i] = groupA;
      }
      const groupFlags = groupA.groupFlags | groupB.groupFlags;
      this.SetGroupFlags(groupA, groupFlags);
      groupA.lastIndex = groupB.lastIndex;
      groupB.firstIndex = groupB.lastIndex;
      this.DestroyParticleGroup(groupB);
    }

    /**
     * Split particle group into multiple disconnected groups.
     *
     * warning: This function is locked during callbacks.
     *
     * @param group the group to be split.
     */
    public SplitParticleGroup(group: ParticleGroup): void {
      this.UpdateContacts(true);
      const particleCount = group.GetParticleCount();
      // We create several linked lists. Each list represents a set of connected particles.
      const nodeBuffer: ParticleSysteParticleListNode[] = MakeArray(particleCount, (index: number) => new ParticleSysteParticleListNode());
      ParticleSystem.InitializeParticleLists(group, nodeBuffer);
      this.MergeParticleListsInContact(group, nodeBuffer);
      const survivingList = ParticleSystem.FindLongestParticleList(group, nodeBuffer);
      this.MergeZombieParticleListNodes(group, nodeBuffer, survivingList);
      this.CreateParticleGroupsFromParticleList(group, nodeBuffer, survivingList);
      this.UpdatePairsAndTriadsWithParticleList(group, nodeBuffer);
    }

    /**
     * Get the world particle group list. With the returned group,
     * use ParticleGroup::GetNext to get the next group in the
     * world list.
     *
     * A null group indicates the end of the list.
     *
     * @return the head of the world particle group list.
     */
    public GetParticleGroupList(): ParticleGroup | null {
      return this.groupList;
    }

    /**
     * Get the number of particle groups.
     */
    public GetParticleGroupCount(): number {
      return this.groupCount;
    }

    /**
     * Get the number of particles.
     */
    public GetParticleCount(): number {
      return this.count;
    }

    /**
     * Get the maximum number of particles.
     */
    public GetMaxParticleCount(): number {
      return this.def.maxCount;
    }

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
    public SetMaxParticleCount(count: number): void {
      // DEBUG: Assert(this.count <= count);
      this.def.maxCount = count;
    }

    /**
     * Get all existing particle flags.
     */
    public GetAllParticleFlags(): ParticleFlag {
      return this.allParticleFlags;
    }

    /**
     * Get all existing particle group flags.
     */
    public GetAllGroupFlags(): ParticleGroupFlag {
      return this.allGroupFlags;
    }

    /**
     * Pause or unpause the particle system. When paused,
     * World::Step() skips over this particle system. All
     * ParticleSystem function calls still work.
     *
     * @param paused paused is true to pause, false to un-pause.
     */
    public SetPaused(paused: boolean): void {
      this.paused = paused;
    }

    /**
     * Initially, true, then, the last value passed into
     * SetPaused().
     *
     * @return true if the particle system is being updated in World::Step().
     */
    public GetPaused(): boolean {
      return this.paused;
    }

    /**
     * Change the particle density.
     *
     * Particle density affects the mass of the particles, which in
     * turn affects how the particles interact with Bodies. Note
     * that the density does not affect how the particles interact
     * with each other.
     */
    public SetDensity(density: number): void {
      this.def.density = density;
      this.inverseDensity = 1 / this.def.density;
    }

    /**
     * Get the particle density.
     */
    public GetDensity(): number {
      return this.def.density;
    }

    /**
     * Change the particle gravity scale. Adjusts the effect of the
     * global gravity vector on particles.
     */
    public SetGravityScale(gravityScale: number): void {
      this.def.gravityScale = gravityScale;
    }

    /**
     * Get the particle gravity scale.
     */
    public GetGravityScale(): number {
      return this.def.gravityScale;
    }

    /**
     * Damping is used to reduce the velocity of particles. The
     * damping parameter can be larger than 1.0f but the damping
     * effect becomes sensitive to the time step when the damping
     * parameter is large.
     */
    public SetDamping(damping: number): void {
      this.def.dampingStrength = damping;
    }

    /**
     * Get damping for particles
     */
    public GetDamping(): number {
      return this.def.dampingStrength;
    }

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
    public SetStaticPressureIterations(iterations: number): void {
      this.def.staticPressureIterations = iterations;
    }

    /**
     * Get the number of iterations for static pressure of
     * particles.
     */
    public GetStaticPressureIterations(): number {
      return this.def.staticPressureIterations;
    }

    /**
     * Change the particle radius.
     *
     * You should set this only once, on world start.
     * If you change the radius during execution, existing particles
     * may explode, shrink, or behave unexpectedly.
     */
    public SetRadius(radius: number): void {
      this.particleDiameter = 2 * radius;
      this.squaredDiameter = this.particleDiameter * this.particleDiameter;
      this.inverseDiameter = 1 / this.particleDiameter;
    }

    /**
     * Get the particle radius.
     */
    public GetRadius(): number {
      return this.particleDiameter / 2;
    }

    /**
     * Get the position of each particle
     *
     * Array is length GetParticleCount()
     *
     * @return the pointer to the head of the particle positions array.
     */
    public GetPositionBuffer(): Vec2[] {
      return this.positionBuffer.data;
    }

    /**
     * Get the velocity of each particle
     *
     * Array is length GetParticleCount()
     *
     * @return the pointer to the head of the particle velocities array.
     */
    public GetVelocityBuffer(): Vec2[] {
      return this.velocityBuffer.data;
    }

    /**
     * Get the color of each particle
     *
     * Array is length GetParticleCount()
     *
     * @return the pointer to the head of the particle colors array.
     */
    public GetColorBuffer(): Color[] {
      this.colorBuffer.data = this.RequestBuffer(this.colorBuffer.data);
      return this.colorBuffer.data;
    }

    /**
     * Get the particle-group of each particle.
     *
     * Array is length GetParticleCount()
     *
     * @return the pointer to the head of the particle group array.
     */
    public GetGroupBuffer(): Array<ParticleGroup | null> {
      return this.groupBuffer;
    }

    /**
     * Get the weight of each particle
     *
     * Array is length GetParticleCount()
     *
     * @return the pointer to the head of the particle positions array.
     */
    public GetWeightBuffer(): number[] {
      return this.weightBuffer;
    }

    /**
     * Get the user-specified data of each particle.
     *
     * Array is length GetParticleCount()
     *
     * @return the pointer to the head of the particle user-data array.
     */
    public GetUserDataBuffer<T>(): T[] {
      this.userDataBuffer.data = this.RequestBuffer(this.userDataBuffer.data);
      return this.userDataBuffer.data;
    }

    /**
     * Get the flags for each particle. See the ParticleFlag enum.
     *
     * Array is length GetParticleCount()
     *
     * @return the pointer to the head of the particle-flags array.
     */
    public GetFlagsBuffer(): ParticleFlag[] {
      return this.flagsBuffer.data;
    }

    /**
     * Set flags for a particle. See the ParticleFlag enum.
     */
    public SetParticleFlags(index: number, newFlags: ParticleFlag): void {
      const oldFlags = this.flagsBuffer.data[index];
      if (oldFlags & ~newFlags) {
        // If any flags might be removed
        this.needsUpdateAllParticleFlags = true;
      }
      if (~this.allParticleFlags & newFlags) {
        // If any flags were added
        if (newFlags & ParticleFlag.tensileParticle) {
          this.accumulation2Buffer = this.RequestBuffer(this.accumulation2Buffer);
        }
        if (newFlags & ParticleFlag.colorMixingParticle) {
          this.colorBuffer.data = this.RequestBuffer(this.colorBuffer.data);
        }
        this.allParticleFlags |= newFlags;
      }
      this.flagsBuffer.data[index] = newFlags;
    }

    /**
     * Get flags for a particle. See the ParticleFlag enum.
     */
    public GetParticleFlags(index: number): ParticleFlag {
      return this.flagsBuffer.data[index];
    }

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
    public SetFlagsBuffer(buffer: ParticleFlag[]): void {
      this.SetUserOverridableBuffer(this.flagsBuffer, buffer);
    }

    public SetPositionBuffer(buffer: Vec2[] | Float32Array): void {
      if (buffer instanceof Float32Array) {
        if (buffer.length % 2 !== 0) { throw new Error(); }
        const count: number = buffer.length / 2;
        const array: TypedVec2[] = new Array(count);
        for (let i = 0; i < count; ++i) {
          array[i] = new TypedVec2(buffer.subarray(i * 2, i * 2 + 2));
        }
        buffer = array;
      }
      this.SetUserOverridableBuffer(this.positionBuffer, buffer);
    }

    public SetVelocityBuffer(buffer: TypedVec2[] | Float32Array): void {
      if (buffer instanceof Float32Array) {
        if (buffer.length % 2 !== 0) { throw new Error(); }
        const count: number = buffer.length / 2;
        const array: TypedVec2[] = new Array(count);
        for (let i = 0; i < count; ++i) {
          array[i] = new TypedVec2(buffer.subarray(i * 2, i * 2 + 2));
        }
        buffer = array;
      }
      this.SetUserOverridableBuffer(this.velocityBuffer, buffer);
    }

    public SetColorBuffer(buffer: Color[] | Float32Array): void {
      if (buffer instanceof Float32Array) {
        if (buffer.length % 4 !== 0) { throw new Error(); }
        const count: number = buffer.length / 4;
        const array: Color[] = new Array(count);
        for (let i = 0; i < count; ++i) {
          array[i] = new TypedColor(buffer.subarray(i * 4, i * 4 + 4));
        }
        buffer = array;
      }
      this.SetUserOverridableBuffer(this.colorBuffer, buffer);
    }

    public SetUserDataBuffer<T>(buffer: T[]): void {
      this.SetUserOverridableBuffer(this.userDataBuffer, buffer);
    }

    /**
     * Get contacts between particles
     * Contact data can be used for many reasons, for example to
     * trigger rendering or audio effects.
     */
    public GetContacts(): ParticleContact[] {
      return this.contactBuffer.data;
    }

    public GetContactCount(): number {
      return this.contactBuffer.count;
    }

    /**
     * Get contacts between particles and bodies
     *
     * Contact data can be used for many reasons, for example to
     * trigger rendering or audio effects.
     */
    public GetBodyContacts(): ParticleBodyContact[] {
      return this.bodyContactBuffer.data;
    }

    public GetBodyContactCount(): number {
      return this.bodyContactBuffer.count;
    }

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
    public GetPairs(): ParticlePair[] {
      return this.pairBuffer.data;
    }

    public GetPairCount(): number {
      return this.pairBuffer.count;
    }

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
    public GetTriads(): ParticleTriad[] {
      return this.triadBuffer.data;
    }

    public GetTriadCount(): number {
      return this.triadBuffer.count;
    }

    /**
     * Set an optional threshold for the maximum number of
     * consecutive particle iterations that a particle may contact
     * multiple bodies before it is considered a candidate for being
     * "stuck". Setting to zero or less disables.
     */
    public SetStuckThreshold(steps: number): void {
      this.stuckThreshold = steps;

      if (steps > 0) {
        this.lastBodyContactStepBuffer.data = this.RequestBuffer(this.lastBodyContactStepBuffer.data);
        this.bodyContactCountBuffer.data = this.RequestBuffer(this.bodyContactCountBuffer.data);
        this.consecutiveContactStepsBuffer.data = this.RequestBuffer(this.consecutiveContactStepsBuffer.data);
      }
    }

    /**
     * Get potentially stuck particles from the last step; the user
     * must decide if they are stuck or not, and if so, delete or
     * move them
     */
    public GetStuckCandidates(): number[] {
      ///return stuckParticleBuffer.Data();
      return this.stuckParticleBuffer.Data();
    }

    /**
     * Get the number of stuck particle candidates from the last
     * step.
     */
    public GetStuckCandidateCount(): number {
      ///return stuckParticleBuffer.GetCount();
      return this.stuckParticleBuffer.GetCount();
    }

    /**
     * Compute the kinetic energy that can be lost by damping force
     */
    public ComputeCollisionEnergy(): number {
      const s_v = ParticleSystem.ComputeCollisionEnergy_s_v;
      const vel_data = this.velocityBuffer.data;
      let suv2 = 0;
      for (let k = 0; k < this.contactBuffer.count; k++) {
        const contact = this.contactBuffer.data[k];
        const a = contact.indexA;
        const b = contact.indexB;
        const n = contact.normal;
        ///Vec2 v = velocityBuffer.data[b] - velocityBuffer.data[a];
        const v = Vec2.SubVV(vel_data[b], vel_data[a], s_v);
        const vn = Vec2.DotVV(v, n);
        if (vn < 0) {
          suv2 += vn * vn;
        }
      }
      return 0.5 * this.GetParticleMass() * suv2;
    }
    public static readonly ComputeCollisionEnergy_s_v = new Vec2();

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
    public SetStrictContactCheck(enabled: boolean): void {
      this.def.strictContactCheck = enabled;
    }

    /**
     * Get the status of the strict contact check.
     */
    public GetStrictContactCheck(): boolean {
      return this.def.strictContactCheck;
    }

    /**
     * Set the lifetime (in seconds) of a particle relative to the
     * current time.  A lifetime of less than or equal to 0.0f
     * results in the particle living forever until it's manually
     * destroyed by the application.
     */
    public SetParticleLifetime(index: number, lifetime: number): void {
      // DEBUG: Assert(this.ValidateParticleIndex(index));
      const initializeExpirationTimes = this.indexByExpirationTimeBuffer.data === null;
      this.expirationTimeBuffer.data = this.RequestBuffer(this.expirationTimeBuffer.data);
      this.indexByExpirationTimeBuffer.data = this.RequestBuffer(this.indexByExpirationTimeBuffer.data);

      // Initialize the inverse mapping buffer.
      if (initializeExpirationTimes) {
        const particleCount = this.GetParticleCount();
        for (let i = 0; i < particleCount; ++i) {
          this.indexByExpirationTimeBuffer.data[i] = i;
        }
      }
      ///const int32 quantizedLifetime = (int32)(lifetime / def.lifetimeGranularity);
      const quantizedLifetime = lifetime / this.def.lifetimeGranularity;
      // Use a negative lifetime so that it's possible to track which
      // of the infinite lifetime particles are older.
      const newExpirationTime = quantizedLifetime > 0.0 ? this.GetQuantizedTimeElapsed() + quantizedLifetime : quantizedLifetime;
      if (newExpirationTime !== this.expirationTimeBuffer.data[index]) {
        this.expirationTimeBuffer.data[index] = newExpirationTime;
        this.expirationTimeBufferRequiresSorting = true;
      }
    }

    /**
     * Get the lifetime (in seconds) of a particle relative to the
     * current time.  A value > 0.0f is returned if the particle is
     * scheduled to be destroyed in the future, values <= 0.0f
     * indicate the particle has an infinite lifetime.
     */
    public GetParticleLifetime(index: number): number {
      // DEBUG: Assert(this.ValidateParticleIndex(index));
      return this.ExpirationTimeToLifetime(this.GetExpirationTimeBuffer()[index]);
    }

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
    public SetDestructionByAge(enable: boolean): void {
      if (enable) {
        this.GetExpirationTimeBuffer();
      }
      this.def.destroyByAge = enable;
    }

    /**
     * Get whether the oldest particle will be destroyed in
     * CreateParticle() when the maximum number of particles are
     * present in the system.
     */
    public GetDestructionByAge(): boolean {
      return this.def.destroyByAge;
    }

    /**
     * Get the array of particle expiration times indexed by
     * particle index.
     *
     * GetParticleCount() items are in the returned array.
     */
    public GetExpirationTimeBuffer(): number[] {
      this.expirationTimeBuffer.data = this.RequestBuffer(this.expirationTimeBuffer.data);
      return this.expirationTimeBuffer.data;
    }

    /**
     * Convert a expiration time value in returned by
     * GetExpirationTimeBuffer() to a time in seconds relative to
     * the current simulation time.
     */
    public ExpirationTimeToLifetime(expirationTime: number): number {
      return (expirationTime > 0 ?
        expirationTime - this.GetQuantizedTimeElapsed() :
        expirationTime) * this.def.lifetimeGranularity;
    }

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
    public GetIndexByExpirationTimeBuffer(): number[] {
      // If particles are present, initialize / reinitialize the lifetime buffer.
      if (this.GetParticleCount()) {
        this.SetParticleLifetime(0, this.GetParticleLifetime(0));
      } else {
        this.indexByExpirationTimeBuffer.data = this.RequestBuffer(this.indexByExpirationTimeBuffer.data);
      }
      return this.indexByExpirationTimeBuffer.data;
    }

    /**
     * Apply an impulse to one particle. This immediately modifies
     * the velocity. Similar to Body::ApplyLinearImpulse.
     *
     * @param index the particle that will be modified.
     * @param impulse impulse the world impulse vector, usually in N-seconds or kg-m/s.
     */
    public ParticleApplyLinearImpulse(index: number, impulse: XY): void {
      this.ApplyLinearImpulse(index, index + 1, impulse);
    }

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
    public ApplyLinearImpulse(firstIndex: number, lastIndex: number, impulse: XY): void {
      const vel_data = this.velocityBuffer.data;
      const numParticles = (lastIndex - firstIndex);
      const totalMass = numParticles * this.GetParticleMass();
      ///const Vec2 velocityDelta = impulse / totalMass;
      const velocityDelta = new Vec2().Copy(impulse).SelfMul(1 / totalMass);
      for (let i = firstIndex; i < lastIndex; i++) {
        ///velocityBuffer.data[i] += velocityDelta;
        vel_data[i].SelfAdd(velocityDelta);
      }
    }

    public static IsSignificantForce(force: XY): boolean {
      return force.x !== 0 || force.y !== 0;
    }

    /**
     * Apply a force to the center of a particle.
     *
     * @param index the particle that will be modified.
     * @param force the world force vector, usually in Newtons (N).
     */
    public ParticleApplyForce(index: number, force: XY): void {
      if (ParticleSystem.IsSignificantForce(force) &&
        this.ForceCanBeApplied(this.flagsBuffer.data[index])) {
        this.PrepareForceBuffer();
        ///forceBuffer[index] += force;
        this.forceBuffer[index].SelfAdd(force);
      }
    }

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
    public ApplyForce(firstIndex: number, lastIndex: number, force: XY): void {
      // Ensure we're not trying to apply force to particles that can't move,
      // such as wall particles.
      // DEBUG: let flags = 0;
      // DEBUG: for (let i = firstIndex; i < lastIndex; i++) {
      // DEBUG:   flags |= this.flagsBuffer.data[i];
      // DEBUG: }
      // DEBUG: Assert(this.ForceCanBeApplied(flags));

      // Early out if force does nothing (optimization).
      ///const Vec2 distributedForce = force / (float32)(lastIndex - firstIndex);
      const distributedForce =  new Vec2().Copy(force).SelfMul(1 / (lastIndex - firstIndex));
      if (ParticleSystem.IsSignificantForce(distributedForce)) {
        this.PrepareForceBuffer();

        // Distribute the force over all the particles.
        for (let i = firstIndex; i < lastIndex; i++) {
          ///forceBuffer[i] += distributedForce;
          this.forceBuffer[i].SelfAdd(distributedForce);
        }
      }
    }

    /**
     * Get the next particle-system in the world's particle-system
     * list.
     */
    public GetNext(): ParticleSystem | null {
      return this.next;
    }

    /**
     * Query the particle system for all particles that potentially
     * overlap the provided AABB.
     * QueryCallback::ShouldQueryParticleSystem is ignored.
     *
     * @param callback a user implemented callback class.
     * @param aabb the query box.
     */
    public QueryAABB(callback: QueryCallback, aabb: AABB): void {
      if (this.proxyBuffer.count === 0) {
        return;
      }
      const beginProxy = 0;
      const endProxy = this.proxyBuffer.count;
      const firstProxy = std_lower_bound(this.proxyBuffer.data, beginProxy, endProxy,
        ParticleSystem.computeTag(
          this.inverseDiameter * aabb.lowerBound.x,
          this.inverseDiameter * aabb.lowerBound.y),
        ParticleSysteProxy.CompareProxyTag);
      const lastProxy = std_upper_bound(this.proxyBuffer.data, firstProxy, endProxy,
        ParticleSystem.computeTag(
          this.inverseDiameter * aabb.upperBound.x,
          this.inverseDiameter * aabb.upperBound.y),
        ParticleSysteProxy.CompareTagProxy);
      const pos_data = this.positionBuffer.data;
      for (let k = firstProxy; k < lastProxy; ++k) {
        const proxy = this.proxyBuffer.data[k];
        const i = proxy.index;
        const p = pos_data[i];
        if (aabb.lowerBound.x < p.x && p.x < aabb.upperBound.x &&
          aabb.lowerBound.y < p.y && p.y < aabb.upperBound.y) {
          if (!callback.ReportParticle(this, i)) {
            break;
          }
        }
      }
    }

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
    public QueryShapeAABB(callback: QueryCallback, shape: Shape, xf: Transform, childIndex: number = 0): void {
      const s_aabb = ParticleSystem.QueryShapeAABB_s_aabb;
      const aabb = s_aabb;
      shape.ComputeAABB(aabb, xf, childIndex);
      this.QueryAABB(callback, aabb);
    }
    public static readonly QueryShapeAABB_s_aabb = new AABB();

    public QueryPointAABB(callback: QueryCallback, point: XY, slop: number = linearSlop): void {
      const s_aabb = ParticleSystem.QueryPointAABB_s_aabb;
      const aabb = s_aabb;
      aabb.lowerBound.Set(point.x - slop, point.y - slop);
      aabb.upperBound.Set(point.x + slop, point.y + slop);
      this.QueryAABB(callback, aabb);
    }
    public static readonly QueryPointAABB_s_aabb = new AABB();

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
    public RayCast(callback: RayCastCallback, point1: XY, point2: XY): void {
      const s_aabb = ParticleSystem.RayCast_s_aabb;
      const s_p = ParticleSystem.RayCast_s_p;
      const s_v = ParticleSystem.RayCast_s_v;
      const s_n = ParticleSystem.RayCast_s_n;
      const s_point = ParticleSystem.RayCast_s_point;
      if (this.proxyBuffer.count === 0) {
        return;
      }
      const pos_data = this.positionBuffer.data;
      const aabb = s_aabb;
      Vec2.MinV(point1, point2, aabb.lowerBound);
      Vec2.MaxV(point1, point2, aabb.upperBound);
      let fraction = 1;
      // solving the following equation:
      // ((1-t)*point1+t*point2-position)^2=diameter^2
      // where t is a potential fraction
      ///Vec2 v = point2 - point1;
      const v = Vec2.SubVV(point2, point1, s_v);
      const v2 = Vec2.DotVV(v, v);
      const enumerator = this.GetInsideBoundsEnumerator(aabb);

      let i: number;
      while ((i = enumerator.GetNext()) >= 0) {
        ///Vec2 p = point1 - positionBuffer.data[i];
        const p = Vec2.SubVV(point1, pos_data[i], s_p);
        const pv = Vec2.DotVV(p, v);
        const p2 = Vec2.DotVV(p, p);
        const determinant = pv * pv - v2 * (p2 - this.squaredDiameter);
        if (determinant >= 0) {
          const sqrtDeterminant = Sqrt(determinant);
          // find a solution between 0 and fraction
          let t = (-pv - sqrtDeterminant) / v2;
          if (t > fraction) {
            continue;
          }
          if (t < 0) {
            t = (-pv + sqrtDeterminant) / v2;
            if (t < 0 || t > fraction) {
              continue;
            }
          }
          ///Vec2 n = p + t * v;
          const n = Vec2.AddVMulSV(p, t, v, s_n);
          n.Normalize();
          ///float32 f = callback.ReportParticle(this, i, point1 + t * v, n, t);
          const f = callback.ReportParticle(this, i, Vec2.AddVMulSV(point1, t, v, s_point), n, t);
          fraction = Min(fraction, f);
          if (fraction <= 0) {
            break;
          }
        }
      }
    }
    public static readonly RayCast_s_aabb = new AABB();
    public static readonly RayCast_s_p = new Vec2();
    public static readonly RayCast_s_v = new Vec2();
    public static readonly RayCast_s_n = new Vec2();
    public static readonly RayCast_s_point = new Vec2();

    /**
     * Compute the axis-aligned bounding box for all particles
     * contained within this particle system.
     * @param aabb Returns the axis-aligned bounding box of the system.
     */
    public ComputeAABB(aabb: AABB): void {
      const particleCount = this.GetParticleCount();
      // DEBUG: Assert(aabb !== null);
      aabb.lowerBound.x = +maxFloat;
      aabb.lowerBound.y = +maxFloat;
      aabb.upperBound.x = -maxFloat;
      aabb.upperBound.y = -maxFloat;

      const pos_data = this.positionBuffer.data;
      for (let i = 0; i < particleCount; i++) {
        const p = pos_data[i];
        Vec2.MinV(aabb.lowerBound, p, aabb.lowerBound);
        Vec2.MaxV(aabb.upperBound, p, aabb.upperBound);
      }
      aabb.lowerBound.x -= this.particleDiameter;
      aabb.lowerBound.y -= this.particleDiameter;
      aabb.upperBound.x += this.particleDiameter;
      aabb.upperBound.y += this.particleDiameter;
    }

    /**
     * All particle types that require creating pairs
     */
    public static readonly k_pairFlags: number = ParticleFlag.springParticle;

    /**
     * All particle types that require creating triads
     */
    public static readonly k_triadFlags = ParticleFlag.elasticParticle;

    /**
     * All particle types that do not produce dynamic pressure
     */
    public static readonly k_noPressureFlags = ParticleFlag.powderParticle | ParticleFlag.tensileParticle;

    /**
     * All particle types that apply extra damping force with bodies
     */
    public static readonly k_extraDampingFlags = ParticleFlag.staticPressureParticle;

    public static readonly k_barrierWallFlags = ParticleFlag.barrierParticle | ParticleFlag.wallParticle;

    public FreeBuffer<T>(b: T[] | null, capacity: number): void {
      if (b === null) {
        return;
      }
      b.length = 0;
    }

    public FreeUserOverridableBuffer<T>(b: ParticleSysteUserOverridableBuffer<T>): void {
      if (b.userSuppliedCapacity === 0) {
        this.FreeBuffer(b.data, this.internalAllocatedCapacity);
      }
    }

    /**
     * Reallocate a buffer
     */
    public ReallocateBuffer3<T>(oldBuffer: T[] | null, oldCapacity: number, newCapacity: number): T[] {
      // Assert(newCapacity > oldCapacity);
      if (newCapacity <= oldCapacity) { throw new Error(); }
      const newBuffer = (oldBuffer) ? oldBuffer.slice() : [];
      newBuffer.length = newCapacity;
      return newBuffer;
    }

    /**
     * Reallocate a buffer
     */
    public ReallocateBuffer5<T>(buffer: T[] | null, userSuppliedCapacity: number, oldCapacity: number, newCapacity: number, deferred: boolean): T[] {
      // Assert(newCapacity > oldCapacity);
      if (newCapacity <= oldCapacity) { throw new Error(); }
      // A 'deferred' buffer is reallocated only if it is not NULL.
      // If 'userSuppliedCapacity' is not zero, buffer is user supplied and must
      // be kept.
      // Assert(!userSuppliedCapacity || newCapacity <= userSuppliedCapacity);
      if (!(!userSuppliedCapacity || newCapacity <= userSuppliedCapacity)) { throw new Error(); }
      if ((!deferred || buffer) && !userSuppliedCapacity) {
        buffer = this.ReallocateBuffer3(buffer, oldCapacity, newCapacity);
      }
      return buffer as any; // TODO: fix this
    }

    /**
     * Reallocate a buffer
     */
    public ReallocateBuffer4<T>(buffer: ParticleSysteUserOverridableBuffer<any>, oldCapacity: number, newCapacity: number, deferred: boolean): T[] {
      // DEBUG: Assert(newCapacity > oldCapacity);
      return this.ReallocateBuffer5(buffer.data, buffer.userSuppliedCapacity, oldCapacity, newCapacity, deferred);
    }

    public RequestBuffer<T>(buffer: T[] | null): T[] {
      if (!buffer) {
        if (this.internalAllocatedCapacity === 0) {
          this.ReallocateInternalAllocatedBuffers(minParticleSystemBufferCapacity);
        }

        buffer = [];
        buffer.length = this.internalAllocatedCapacity;
      }
      return buffer;
    }

    /**
     * Reallocate the handle / index map and schedule the allocation
     * of a new pool for handle allocation.
     */
    public ReallocateHandleBuffers(newCapacity: number): void {
      // DEBUG: Assert(newCapacity > this.internalAllocatedCapacity);
      // Reallocate a new handle / index map buffer, copying old handle pointers
      // is fine since they're kept around.
      this.handleIndexBuffer.data = this.ReallocateBuffer4(this.handleIndexBuffer, this.internalAllocatedCapacity, newCapacity, true);
      // Set the size of the next handle allocation.
      ///this.handleAllocator.SetItemsPerSlab(newCapacity - this.internalAllocatedCapacity);
    }

    public ReallocateInternalAllocatedBuffers(capacity: number): void {
      function LimitCapacity(capacity: number, maxCount: number): number {
        return maxCount && capacity > maxCount ? maxCount : capacity;
      }

      // Don't increase capacity beyond the smallest user-supplied buffer size.
      capacity = LimitCapacity(capacity, this.def.maxCount);
      capacity = LimitCapacity(capacity, this.flagsBuffer.userSuppliedCapacity);
      capacity = LimitCapacity(capacity, this.positionBuffer.userSuppliedCapacity);
      capacity = LimitCapacity(capacity, this.velocityBuffer.userSuppliedCapacity);
      capacity = LimitCapacity(capacity, this.colorBuffer.userSuppliedCapacity);
      capacity = LimitCapacity(capacity, this.userDataBuffer.userSuppliedCapacity);
      if (this.internalAllocatedCapacity < capacity) {
        this.ReallocateHandleBuffers(capacity);
        this.flagsBuffer.data = this.ReallocateBuffer4(this.flagsBuffer, this.internalAllocatedCapacity, capacity, false);

        // Conditionally defer these as they are optional if the feature is
        // not enabled.
        const stuck = this.stuckThreshold > 0;
        this.lastBodyContactStepBuffer.data = this.ReallocateBuffer4(this.lastBodyContactStepBuffer, this.internalAllocatedCapacity, capacity, stuck);
        this.bodyContactCountBuffer.data = this.ReallocateBuffer4(this.bodyContactCountBuffer, this.internalAllocatedCapacity, capacity, stuck);
        this.consecutiveContactStepsBuffer.data = this.ReallocateBuffer4(this.consecutiveContactStepsBuffer, this.internalAllocatedCapacity, capacity, stuck);
        this.positionBuffer.data = this.ReallocateBuffer4(this.positionBuffer, this.internalAllocatedCapacity, capacity, false);
        this.velocityBuffer.data = this.ReallocateBuffer4(this.velocityBuffer, this.internalAllocatedCapacity, capacity, false);
        this.forceBuffer = this.ReallocateBuffer5(this.forceBuffer, 0, this.internalAllocatedCapacity, capacity, false);
        this.weightBuffer = this.ReallocateBuffer5(this.weightBuffer, 0, this.internalAllocatedCapacity, capacity, false);
        this.staticPressureBuffer = this.ReallocateBuffer5(this.staticPressureBuffer, 0, this.internalAllocatedCapacity, capacity, true);
        this.accumulationBuffer = this.ReallocateBuffer5(this.accumulationBuffer, 0, this.internalAllocatedCapacity, capacity, false);
        this.accumulation2Buffer = this.ReallocateBuffer5(this.accumulation2Buffer, 0, this.internalAllocatedCapacity, capacity, true);
        this.depthBuffer = this.ReallocateBuffer5(this.depthBuffer, 0, this.internalAllocatedCapacity, capacity, true);
        this.colorBuffer.data = this.ReallocateBuffer4(this.colorBuffer, this.internalAllocatedCapacity, capacity, true);
        this.groupBuffer = this.ReallocateBuffer5(this.groupBuffer, 0, this.internalAllocatedCapacity, capacity, false);
        this.userDataBuffer.data = this.ReallocateBuffer4(this.userDataBuffer, this.internalAllocatedCapacity, capacity, true);
        this.expirationTimeBuffer.data = this.ReallocateBuffer4(this.expirationTimeBuffer, this.internalAllocatedCapacity, capacity, true);
        this.indexByExpirationTimeBuffer.data = this.ReallocateBuffer4(this.indexByExpirationTimeBuffer, this.internalAllocatedCapacity, capacity, false);
        this.internalAllocatedCapacity = capacity;
      }
    }

    public CreateParticleForGroup(groupDef: IParticleGroupDef, xf: Transform, p: XY): void {
      const particleDef = new ParticleDef();
      particleDef.flags = Maybe(groupDef.flags, 0);
      ///particleDef.position = Mul(xf, p);
      Transform.MulXV(xf, p, particleDef.position);
      ///particleDef.velocity =
      ///  groupDef.linearVelocity +
      ///  Cross(groupDef.angularVelocity,
      ///      particleDef.position - groupDef.position);
      Vec2.AddVV(
        Maybe(groupDef.linearVelocity, Vec2.ZERO),
        Vec2.CrossSV(
          Maybe(groupDef.angularVelocity, 0),
          Vec2.SubVV(
            particleDef.position,
            Maybe(groupDef.position, Vec2.ZERO),
            Vec2.s_t0,
          ),
          Vec2.s_t0,
        ),
        particleDef.velocity,
      );
      particleDef.color.Copy(Maybe(groupDef.color, Color.ZERO));
      particleDef.lifetime = Maybe(groupDef.lifetime, 0);
      particleDef.userData = groupDef.userData;
      this.CreateParticle(particleDef);
    }

    public CreateParticlesStrokeShapeForGroup(shape: Shape, groupDef: IParticleGroupDef, xf: Transform): void {
      const s_edge = ParticleSystem.CreateParticlesStrokeShapeForGroup_s_edge;
      const s_d = ParticleSystem.CreateParticlesStrokeShapeForGroup_s_d;
      const s_p = ParticleSystem.CreateParticlesStrokeShapeForGroup_s_p;
      let stride = Maybe(groupDef.stride, 0);
      if (stride === 0) {
        stride = this.GetParticleStride();
      }
      let positionOnEdge = 0;
      const childCount = shape.GetChildCount();
      for (let childIndex = 0; childIndex < childCount; childIndex++) {
        let edge: EdgeShape | null = null;
        if (shape.GetType() === ShapeType.e_edgeShape) {
          edge = shape as EdgeShape;
        } else {
          // DEBUG: Assert(shape.GetType() === ShapeType.e_chainShape);
          edge = s_edge;
          (shape as ChainShape).GetChildEdge(edge, childIndex);
        }
        const d = Vec2.SubVV(edge.vertex2, edge.vertex1, s_d);
        const edgeLength = d.Length();

        while (positionOnEdge < edgeLength) {
          ///Vec2 p = edge.vertex1 + positionOnEdge / edgeLength * d;
          const p = Vec2.AddVMulSV(edge.vertex1, positionOnEdge / edgeLength, d, s_p);
          this.CreateParticleForGroup(groupDef, xf, p);
          positionOnEdge += stride;
        }
        positionOnEdge -= edgeLength;
      }
    }
    public static readonly CreateParticlesStrokeShapeForGroup_s_edge = new EdgeShape();
    public static readonly CreateParticlesStrokeShapeForGroup_s_d = new Vec2();
    public static readonly CreateParticlesStrokeShapeForGroup_s_p = new Vec2();

    public CreateParticlesFillShapeForGroup(shape: Shape, groupDef: IParticleGroupDef, xf: Transform): void {
      const s_aabb = ParticleSystem.CreateParticlesFillShapeForGroup_s_aabb;
      const s_p = ParticleSystem.CreateParticlesFillShapeForGroup_s_p;
      let stride = Maybe(groupDef.stride, 0);
      if (stride === 0) {
        stride = this.GetParticleStride();
      }
      ///Transform identity;
      /// identity.SetIdentity();
      const identity = Transform.IDENTITY;
      const aabb = s_aabb;
      // DEBUG: Assert(shape.GetChildCount() === 1);
      shape.ComputeAABB(aabb, identity, 0);
      for (let y = Math.floor(aabb.lowerBound.y / stride) * stride; y < aabb.upperBound.y; y += stride) {
        for (let x = Math.floor(aabb.lowerBound.x / stride) * stride; x < aabb.upperBound.x; x += stride) {
          const p = s_p.Set(x, y);
          if (shape.TestPoint(identity, p)) {
            this.CreateParticleForGroup(groupDef, xf, p);
          }
        }
      }
    }
    public static readonly CreateParticlesFillShapeForGroup_s_aabb = new AABB();
    public static readonly CreateParticlesFillShapeForGroup_s_p = new Vec2();

    public CreateParticlesWithShapeForGroup(shape: Shape, groupDef: IParticleGroupDef, xf: Transform): void {
      switch (shape.GetType()) {
        case ShapeType.e_edgeShape:
        case ShapeType.e_chainShape:
          this.CreateParticlesStrokeShapeForGroup(shape, groupDef, xf);
          break;
        case ShapeType.e_polygonShape:
        case ShapeType.e_circleShape:
          this.CreateParticlesFillShapeForGroup(shape, groupDef, xf);
          break;
        default:
          // DEBUG: Assert(false);
          break;
      }
    }

    public CreateParticlesWithShapesForGroup(shapes: Shape[], shapeCount: number, groupDef: IParticleGroupDef, xf: Transform): void {
      const compositeShape = new ParticleSysteCompositeShape(shapes, shapeCount);
      this.CreateParticlesFillShapeForGroup(compositeShape, groupDef, xf);
    }

    public CloneParticle(oldIndex: number, group: ParticleGroup): number {
      const def = new ParticleDef();
      def.flags = this.flagsBuffer.data[oldIndex];
      def.position.Copy(this.positionBuffer.data[oldIndex]);
      def.velocity.Copy(this.velocityBuffer.data[oldIndex]);
      if (this.colorBuffer.data) {
        def.color.Copy(this.colorBuffer.data[oldIndex]);
      }
      if (this.userDataBuffer.data) {
        def.userData = this.userDataBuffer.data[oldIndex];
      }
      def.group = group;
      const newIndex = this.CreateParticle(def);
      if (this.handleIndexBuffer.data) {
        const handle = this.handleIndexBuffer.data[oldIndex];
        if (handle) { handle.SetIndex(newIndex); }
        this.handleIndexBuffer.data[newIndex] = handle;
        this.handleIndexBuffer.data[oldIndex] = null;
      }
      if (this.lastBodyContactStepBuffer.data) {
        this.lastBodyContactStepBuffer.data[newIndex] =
          this.lastBodyContactStepBuffer.data[oldIndex];
      }
      if (this.bodyContactCountBuffer.data) {
        this.bodyContactCountBuffer.data[newIndex] =
          this.bodyContactCountBuffer.data[oldIndex];
      }
      if (this.consecutiveContactStepsBuffer.data) {
        this.consecutiveContactStepsBuffer.data[newIndex] =
          this.consecutiveContactStepsBuffer.data[oldIndex];
      }
      if (this.hasForce) {
        this.forceBuffer[newIndex].Copy(this.forceBuffer[oldIndex]);
      }
      if (this.staticPressureBuffer) {
        this.staticPressureBuffer[newIndex] = this.staticPressureBuffer[oldIndex];
      }
      if (this.depthBuffer) {
        this.depthBuffer[newIndex] = this.depthBuffer[oldIndex];
      }
      if (this.expirationTimeBuffer.data) {
        this.expirationTimeBuffer.data[newIndex] =
          this.expirationTimeBuffer.data[oldIndex];
      }
      return newIndex;
    }

    public DestroyParticlesInGroup(group: ParticleGroup, callDestructionListener: boolean = false): void {
      for (let i = group.firstIndex; i < group.lastIndex; i++) {
        this.DestroyParticle(i, callDestructionListener);
      }
    }

    public DestroyParticleGroup(group: ParticleGroup): void {
      // DEBUG: Assert(this.groupCount > 0);
      // DEBUG: Assert(group !== null);

      if (this.world.destructionListener) {
        this.world.destructionListener.SayGoodbyeParticleGroup(group);
      }

      this.SetGroupFlags(group, 0);
      for (let i = group.firstIndex; i < group.lastIndex; i++) {
        this.groupBuffer[i] = null;
      }

      if (group.prev) {
        group.prev.next = group.next;
      }
      if (group.next) {
        group.next.prev = group.prev;
      }
      if (group === this.groupList) {
        this.groupList = group.next;
      }

      --this.groupCount;
    }

    public static ParticleCanBeConnected(flags: ParticleFlag, group: ParticleGroup | null): boolean {
      return ((flags & (ParticleFlag.wallParticle | ParticleFlag.springParticle | ParticleFlag.elasticParticle)) !== 0) ||
        ((group !== null) && ((group.GetGroupFlags() & ParticleGroupFlag.rigidParticleGroup) !== 0));
    }

    public UpdatePairsAndTriads(firstIndex: number, lastIndex: number, filter: ParticleSysteConnectionFilter): void {
      const s_dab = ParticleSystem.UpdatePairsAndTriads_s_dab;
      const s_dbc = ParticleSystem.UpdatePairsAndTriads_s_dbc;
      const s_dca = ParticleSystem.UpdatePairsAndTriads_s_dca;
      const pos_data = this.positionBuffer.data;
      // Create pairs or triads.
      // All particles in each pair/triad should satisfy the following:
      // * firstIndex <= index < lastIndex
      // * don't have zombieParticle
      // * ParticleCanBeConnected returns true
      // * ShouldCreatePair/ShouldCreateTriad returns true
      // Any particles in each pair/triad should satisfy the following:
      // * filter.IsNeeded returns true
      // * have one of k_pairFlags/k_triadsFlags
      // DEBUG: Assert(firstIndex <= lastIndex);
      let particleFlags = 0;
      for (let i = firstIndex; i < lastIndex; i++) {
        particleFlags |= this.flagsBuffer.data[i];
      }
      if (particleFlags & ParticleSystem.k_pairFlags) {
        for (let k = 0; k < this.contactBuffer.count; k++) {
          const contact = this.contactBuffer.data[k];
          const a = contact.indexA;
          const b = contact.indexB;
          const af = this.flagsBuffer.data[a];
          const bf = this.flagsBuffer.data[b];
          const groupA = this.groupBuffer[a];
          const groupB = this.groupBuffer[b];
          if (a >= firstIndex && a < lastIndex &&
            b >= firstIndex && b < lastIndex &&
            !((af | bf) & ParticleFlag.zombieParticle) &&
            ((af | bf) & ParticleSystem.k_pairFlags) &&
            (filter.IsNecessary(a) || filter.IsNecessary(b)) &&
            ParticleSystem.ParticleCanBeConnected(af, groupA) &&
            ParticleSystem.ParticleCanBeConnected(bf, groupB) &&
            filter.ShouldCreatePair(a, b)) {
            ///ParticlePair& pair = pairBuffer.Append();
            const pair = this.pairBuffer.data[this.pairBuffer.Append()];
            pair.indexA = a;
            pair.indexB = b;
            pair.flags = contact.flags;
            pair.strength = Min(
              groupA ? groupA.strength : 1,
              groupB ? groupB.strength : 1);
            ///pair.distance = Distance(pos_data[a], pos_data[b]); // TODO: this was wrong!
            pair.distance = Vec2.DistanceVV(pos_data[a], pos_data[b]);
          }
          ///std::stable_sort(pairBuffer.Begin(), pairBuffer.End(), ComparePairIndices);
          std_stable_sort(this.pairBuffer.data, 0, this.pairBuffer.count, ParticleSystem.ComparePairIndices);
          ///pairBuffer.Unique(MatchPairIndices);
          this.pairBuffer.Unique(ParticleSystem.MatchPairIndices);
        }
      }
      if (particleFlags & ParticleSystem.k_triadFlags) {
        const diagram = new VoronoiDiagram(lastIndex - firstIndex);
        ///let necessary_count = 0;
        for (let i = firstIndex; i < lastIndex; i++) {
          const flags = this.flagsBuffer.data[i];
          const group = this.groupBuffer[i];
          if (!(flags & ParticleFlag.zombieParticle) &&
            ParticleSystem.ParticleCanBeConnected(flags, group)) {
            ///if (filter.IsNecessary(i)) {
            ///++necessary_count;
            ///}
            diagram.AddGenerator(pos_data[i], i, filter.IsNecessary(i));
          }
        }
        ///if (necessary_count === 0) {
        /////debugger;
        ///for (let i = firstIndex; i < lastIndex; i++) {
        ///  filter.IsNecessary(i);
        ///}
        ///}
        const stride = this.GetParticleStride();
        diagram.Generate(stride / 2, stride * 2);
        const system = this;
        const callback = /*UpdateTriadsCallback*/(a: number, b: number, c: number): void => {
          const af = system.flagsBuffer.data[a];
          const bf = system.flagsBuffer.data[b];
          const cf = system.flagsBuffer.data[c];
          if (((af | bf | cf) & ParticleSystem.k_triadFlags) &&
            filter.ShouldCreateTriad(a, b, c)) {
            const pa = pos_data[a];
            const pb = pos_data[b];
            const pc = pos_data[c];
            const dab = Vec2.SubVV(pa, pb, s_dab);
            const dbc = Vec2.SubVV(pb, pc, s_dbc);
            const dca = Vec2.SubVV(pc, pa, s_dca);
            const maxDistanceSquared = maxTriadDistanceSquared * system.squaredDiameter;
            if (Vec2.DotVV(dab, dab) > maxDistanceSquared ||
              Vec2.DotVV(dbc, dbc) > maxDistanceSquared ||
              Vec2.DotVV(dca, dca) > maxDistanceSquared) {
              return;
            }
            const groupA = system.groupBuffer[a];
            const groupB = system.groupBuffer[b];
            const groupC = system.groupBuffer[c];
            ///ParticleTriad& triad = system.triadBuffer.Append();
            const triad = system.triadBuffer.data[system.triadBuffer.Append()];
            triad.indexA = a;
            triad.indexB = b;
            triad.indexC = c;
            triad.flags = af | bf | cf;
            triad.strength = Min(Min(
              groupA ? groupA.strength : 1,
              groupB ? groupB.strength : 1),
              groupC ? groupC.strength : 1);
            ///let midPoint = Vec2.MulSV(1.0 / 3.0, Vec2.AddVV(pa, Vec2.AddVV(pb, pc, new Vec2()), new Vec2()), new Vec2());
            const midPoint_x = (pa.x + pb.x + pc.x) / 3.0;
            const midPoint_y = (pa.y + pb.y + pc.y) / 3.0;
            ///triad.pa = Vec2.SubVV(pa, midPoint, new Vec2());
            triad.pa.x = pa.x - midPoint_x;
            triad.pa.y = pa.y - midPoint_y;
            ///triad.pb = Vec2.SubVV(pb, midPoint, new Vec2());
            triad.pb.x = pb.x - midPoint_x;
            triad.pb.y = pb.y - midPoint_y;
            ///triad.pc = Vec2.SubVV(pc, midPoint, new Vec2());
            triad.pc.x = pc.x - midPoint_x;
            triad.pc.y = pc.y - midPoint_y;
            triad.ka = -Vec2.DotVV(dca, dab);
            triad.kb = -Vec2.DotVV(dab, dbc);
            triad.kc = -Vec2.DotVV(dbc, dca);
            triad.s = Vec2.CrossVV(pa, pb) + Vec2.CrossVV(pb, pc) + Vec2.CrossVV(pc, pa);
          }
        };
        diagram.GetNodes(callback);
        ///std::stable_sort(triadBuffer.Begin(), triadBuffer.End(), CompareTriadIndices);
        std_stable_sort(this.triadBuffer.data, 0, this.triadBuffer.count, ParticleSystem.CompareTriadIndices);
        ///triadBuffer.Unique(MatchTriadIndices);
        this.triadBuffer.Unique(ParticleSystem.MatchTriadIndices);
      }
    }
    private static UpdatePairsAndTriads_s_dab = new Vec2();
    private static UpdatePairsAndTriads_s_dbc = new Vec2();
    private static UpdatePairsAndTriads_s_dca = new Vec2();

    public UpdatePairsAndTriadsWithReactiveParticles(): void {
      const filter = new ParticleSysteReactiveFilter(this.flagsBuffer);
      this.UpdatePairsAndTriads(0, this.count, filter);

      for (let i = 0; i < this.count; i++) {
        this.flagsBuffer.data[i] &= ~ParticleFlag.reactiveParticle;
      }
      this.allParticleFlags &= ~ParticleFlag.reactiveParticle;
    }

    public static ComparePairIndices(a: ParticlePair, b: ParticlePair): boolean {
      const diffA = a.indexA - b.indexA;
      if (diffA !== 0) { return diffA < 0; }
      return a.indexB < b.indexB;
    }

    public static MatchPairIndices(a: ParticlePair, b: ParticlePair): boolean {
      return a.indexA === b.indexA && a.indexB === b.indexB;
    }

    public static CompareTriadIndices(a: ParticleTriad, b: ParticleTriad): boolean {
      const diffA = a.indexA - b.indexA;
      if (diffA !== 0) { return diffA < 0; }
      const diffB = a.indexB - b.indexB;
      if (diffB !== 0) { return diffB < 0; }
      return a.indexC < b.indexC;
    }

    public static MatchTriadIndices(a: ParticleTriad, b: ParticleTriad): boolean {
      return a.indexA === b.indexA && a.indexB === b.indexB && a.indexC === b.indexC;
    }

    public static InitializeParticleLists(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[]): void {
      const bufferIndex = group.GetBufferIndex();
      const particleCount = group.GetParticleCount();
      for (let i = 0; i < particleCount; i++) {
        const node: ParticleSysteParticleListNode = nodeBuffer[i];
        node.list = node;
        node.next = null;
        node.count = 1;
        node.index = i + bufferIndex;
      }
    }

    public MergeParticleListsInContact(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[]): void {
      const bufferIndex = group.GetBufferIndex();
      for (let k = 0; k < this.contactBuffer.count; k++) {
        /*const ParticleContact&*/
        const contact = this.contactBuffer.data[k];
        const a = contact.indexA;
        const b = contact.indexB;
        if (!group.ContainsParticle(a) || !group.ContainsParticle(b)) {
          continue;
        }
        let listA: ParticleSysteParticleListNode = nodeBuffer[a - bufferIndex].list;
        let listB: ParticleSysteParticleListNode = nodeBuffer[b - bufferIndex].list;
        if (listA === listB) {
          continue;
        }
        // To minimize the cost of insertion, make sure listA is longer than
        // listB.
        if (listA.count < listB.count) {
          const _tmp = listA;
          listA = listB;
          listB = _tmp; ///Swap(listA, listB);
        }
        // DEBUG: Assert(listA.count >= listB.count);
        ParticleSystem.MergeParticleLists(listA, listB);
      }
    }

    public static MergeParticleLists(listA: ParticleSysteParticleListNode, listB: ParticleSysteParticleListNode): void {
      // Insert listB between index 0 and 1 of listA
      // Example:
      //     listA => a1 => a2 => a3 => null
      //     listB => b1 =>  => null
      // to
      //     listA => listB => b1 =>  => a1 => a2 => a3 => null
      // DEBUG: Assert(listA !== listB);
      for (let b: ParticleSysteParticleListNode = listB; ; ) {
        b.list = listA;
        const nextB: ParticleSysteParticleListNode | null = b.next;
        if (nextB) {
          b = nextB;
        } else {
          b.next = listA.next;
          break;
        }
      }
      listA.next = listB;
      listA.count += listB.count;
      listB.count = 0;
    }

    public static FindLongestParticleList(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[]): ParticleSysteParticleListNode {
      const particleCount = group.GetParticleCount();
      let result: ParticleSysteParticleListNode = nodeBuffer[0];
      for (let i = 0; i < particleCount; i++) {
        const node: ParticleSysteParticleListNode = nodeBuffer[i];
        if (result.count < node.count) {
          result = node;
        }
      }
      return result;
    }

    public MergeZombieParticleListNodes(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[], survivingList: ParticleSysteParticleListNode): void {
      const particleCount = group.GetParticleCount();
      for (let i = 0; i < particleCount; i++) {
        const node: ParticleSysteParticleListNode = nodeBuffer[i];
        if (node !== survivingList &&
          (this.flagsBuffer.data[node.index] & ParticleFlag.zombieParticle)) {
          ParticleSystem.MergeParticleListAndNode(survivingList, node);
        }
      }
    }

    public static MergeParticleListAndNode(list: ParticleSysteParticleListNode, node: ParticleSysteParticleListNode): void {
      // Insert node between index 0 and 1 of list
      // Example:
      //     list => a1 => a2 => a3 => null
      //     node => null
      // to
      //     list => node => a1 => a2 => a3 => null
      // DEBUG: Assert(node !== list);
      // DEBUG: Assert(node.list === node);
      // DEBUG: Assert(node.count === 1);
      node.list = list;
      node.next = list.next;
      list.next = node;
      list.count++;
      node.count = 0;
    }

    public CreateParticleGroupsFromParticleList(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[], survivingList: ParticleSysteParticleListNode): void {
      const particleCount = group.GetParticleCount();
      const def = new ParticleGroupDef();
      def.groupFlags = group.GetGroupFlags();
      def.userData = group.GetUserData();
      for (let i = 0; i < particleCount; i++) {
        const list: ParticleSysteParticleListNode = nodeBuffer[i];
        if (!list.count || list === survivingList) {
          continue;
        }
        // DEBUG: Assert(list.list === list);
        const newGroup: ParticleGroup = this.CreateParticleGroup(def);
        for (let node: ParticleSysteParticleListNode | null = list; node; node = node.next) {
          const oldIndex = node.index;
          // DEBUG: const flags = this.flagsBuffer.data[oldIndex];
          // DEBUG: Assert(!(flags & ParticleFlag.zombieParticle));
          const newIndex = this.CloneParticle(oldIndex, newGroup);
          this.flagsBuffer.data[oldIndex] |= ParticleFlag.zombieParticle;
          node.index = newIndex;
        }
      }
    }

    public UpdatePairsAndTriadsWithParticleList(group: ParticleGroup, nodeBuffer: ParticleSysteParticleListNode[]): void {
      const bufferIndex = group.GetBufferIndex();
      // Update indices in pairs and triads. If an index belongs to the group,
      // replace it with the corresponding value in nodeBuffer.
      // Note that nodeBuffer is allocated only for the group and the index should
      // be shifted by bufferIndex.
      for (let k = 0; k < this.pairBuffer.count; k++) {
        const pair = this.pairBuffer.data[k];
        const a = pair.indexA;
        const b = pair.indexB;
        if (group.ContainsParticle(a)) {
          pair.indexA = nodeBuffer[a - bufferIndex].index;
        }
        if (group.ContainsParticle(b)) {
          pair.indexB = nodeBuffer[b - bufferIndex].index;
        }
      }
      for (let k = 0; k < this.triadBuffer.count; k++) {
        const triad = this.triadBuffer.data[k];
        const a = triad.indexA;
        const b = triad.indexB;
        const c = triad.indexC;
        if (group.ContainsParticle(a)) {
          triad.indexA = nodeBuffer[a - bufferIndex].index;
        }
        if (group.ContainsParticle(b)) {
          triad.indexB = nodeBuffer[b - bufferIndex].index;
        }
        if (group.ContainsParticle(c)) {
          triad.indexC = nodeBuffer[c - bufferIndex].index;
        }
      }
    }

    public ComputeDepth(): void {
      const contactGroups: ParticleContact[] = []; // TODO: static
      let contactGroupsCount = 0;
      for (let k = 0; k < this.contactBuffer.count; k++) {
        const contact = this.contactBuffer.data[k];
        const a = contact.indexA;
        const b = contact.indexB;
        const groupA = this.groupBuffer[a];
        const groupB = this.groupBuffer[b];
        if (groupA && groupA === groupB &&
          (groupA.groupFlags & ParticleGroupFlag.particleGroupNeedsUpdateDepth)) {
          contactGroups[contactGroupsCount++] = contact;
        }
      }
      const groupsToUpdate: ParticleGroup[] = []; // TODO: static
      let groupsToUpdateCount = 0;
      for (let group = this.groupList; group; group = group.GetNext()) {
        if (group.groupFlags & ParticleGroupFlag.particleGroupNeedsUpdateDepth) {
          groupsToUpdate[groupsToUpdateCount++] = group;
          this.SetGroupFlags(group,
            group.groupFlags &
            ~ParticleGroupFlag.particleGroupNeedsUpdateDepth);
          for (let i = group.firstIndex; i < group.lastIndex; i++) {
            this.accumulationBuffer[i] = 0;
          }
        }
      }
      // Compute sum of weight of contacts except between different groups.
      for (let k = 0; k < contactGroupsCount; k++) {
        const contact = contactGroups[k];
        const a = contact.indexA;
        const b = contact.indexB;
        const w = contact.weight;
        this.accumulationBuffer[a] += w;
        this.accumulationBuffer[b] += w;
      }

      // DEBUG: Assert(this.depthBuffer !== null);
      for (let i = 0; i < groupsToUpdateCount; i++) {
        const group = groupsToUpdate[i];
        for (let i = group.firstIndex; i < group.lastIndex; i++) {
          const w = this.accumulationBuffer[i];
          this.depthBuffer[i] = w < 0.8 ? 0 : maxFloat;
        }
      }
      // The number of iterations is equal to particle number from the deepest
      // particle to the nearest surface particle, and in general it is smaller
      // than sqrt of total particle number.
      ///int32 iterationCount = (int32)Sqrt((float)count);
      const iterationCount = Sqrt(this.count) >> 0;
      for (let t = 0; t < iterationCount; t++) {
        let updated = false;
        for (let k = 0; k < contactGroupsCount; k++) {
          const contact = contactGroups[k];
          const a = contact.indexA;
          const b = contact.indexB;
          const r = 1 - contact.weight;
          ///float32& ap0 = depthBuffer[a];
          const ap0 = this.depthBuffer[a];
          ///float32& bp0 = depthBuffer[b];
          const bp0 = this.depthBuffer[b];
          const ap1 = bp0 + r;
          const bp1 = ap0 + r;
          if (ap0 > ap1) {
            ///ap0 = ap1;
            this.depthBuffer[a] = ap1;
            updated = true;
          }
          if (bp0 > bp1) {
            ///bp0 = bp1;
            this.depthBuffer[b] = bp1;
            updated = true;
          }
        }
        if (!updated) {
          break;
        }
      }
      for (let i = 0; i < groupsToUpdateCount; i++) {
        const group = groupsToUpdate[i];
        for (let i = group.firstIndex; i < group.lastIndex; i++) {
          if (this.depthBuffer[i] < maxFloat) {
            this.depthBuffer[i] *= this.particleDiameter;
          } else {
            this.depthBuffer[i] = 0;
          }
        }
      }
    }

    public GetInsideBoundsEnumerator(aabb: AABB): ParticleSysteInsideBoundsEnumerator {
      const lowerTag = ParticleSystem.computeTag(this.inverseDiameter * aabb.lowerBound.x - 1,
        this.inverseDiameter * aabb.lowerBound.y - 1);
      const upperTag = ParticleSystem.computeTag(this.inverseDiameter * aabb.upperBound.x + 1,
        this.inverseDiameter * aabb.upperBound.y + 1);
      ///const Proxy* beginProxy = proxyBuffer.Begin();
      const beginProxy = 0;
      ///const Proxy* endProxy = proxyBuffer.End();
      const endProxy = this.proxyBuffer.count;
      ///const Proxy* firstProxy = std::lower_bound(beginProxy, endProxy, lowerTag);
      const firstProxy = std_lower_bound(this.proxyBuffer.data, beginProxy, endProxy, lowerTag, ParticleSysteProxy.CompareProxyTag);
      ///const Proxy* lastProxy = std::upper_bound(firstProxy, endProxy, upperTag);
      const lastProxy = std_upper_bound(this.proxyBuffer.data, beginProxy, endProxy, upperTag, ParticleSysteProxy.CompareTagProxy);

      // DEBUG: Assert(beginProxy <= firstProxy);
      // DEBUG: Assert(firstProxy <= lastProxy);
      // DEBUG: Assert(lastProxy <= endProxy);

      return new ParticleSysteInsideBoundsEnumerator(this, lowerTag, upperTag, firstProxy, lastProxy);
    }

    public UpdateAllParticleFlags(): void {
      this.allParticleFlags = 0;
      for (let i = 0; i < this.count; i++) {
        this.allParticleFlags |= this.flagsBuffer.data[i];
      }
      this.needsUpdateAllParticleFlags = false;
    }

    public UpdateAllGroupFlags(): void {
      this.allGroupFlags = 0;
      for (let group = this.groupList; group; group = group.GetNext()) {
        this.allGroupFlags |= group.groupFlags;
      }
      this.needsUpdateAllGroupFlags = false;
    }

    public AddContact(a: number, b: number, contacts: GrowableBuffer<ParticleContact>): void {
      // DEBUG: Assert(contacts === this.contactBuffer);
      const flags_data = this.flagsBuffer.data;
      const pos_data = this.positionBuffer.data;
      ///Vec2 d = positionBuffer.data[b] - positionBuffer.data[a];
      const d = Vec2.SubVV(pos_data[b], pos_data[a], ParticleSystem.AddContact_s_d);
      const distBtParticlesSq = Vec2.DotVV(d, d);
      if (0 < distBtParticlesSq && distBtParticlesSq < this.squaredDiameter) {
        const invD = InvSqrt(distBtParticlesSq);
        ///ParticleContact& contact = contacts.Append();
        const contact = this.contactBuffer.data[this.contactBuffer.Append()];
        contact.indexA = a;
        contact.indexB = b;
        contact.flags = flags_data[a] | flags_data[b];
        contact.weight = 1 - distBtParticlesSq * invD * this.inverseDiameter;
        contact.normal.x = invD * d.x;
        contact.normal.y = invD * d.y;
      }
    }
    public static readonly AddContact_s_d = new Vec2();

    public FindContacts_Reference(contacts: GrowableBuffer<ParticleContact>): void {
      // DEBUG: Assert(contacts === this.contactBuffer);
      const beginProxy = 0;
      const endProxy = this.proxyBuffer.count;

      this.contactBuffer.count = 0;
      for (let a = beginProxy, c = beginProxy; a < endProxy; a++) {
        const rightTag = ParticleSystem.computeRelativeTag(this.proxyBuffer.data[a].tag, 1, 0);
        for (let b = a + 1; b < endProxy; b++) {
          if (rightTag < this.proxyBuffer.data[b].tag) { break; }
          this.AddContact(this.proxyBuffer.data[a].index, this.proxyBuffer.data[b].index, this.contactBuffer);
        }
        const bottomLeftTag = ParticleSystem.computeRelativeTag(this.proxyBuffer.data[a].tag, -1, 1);
        for (; c < endProxy; c++) {
          if (bottomLeftTag <= this.proxyBuffer.data[c].tag) { break; }
        }
        const bottomRightTag = ParticleSystem.computeRelativeTag(this.proxyBuffer.data[a].tag, 1, 1);
        for (let b = c; b < endProxy; b++) {
          if (bottomRightTag < this.proxyBuffer.data[b].tag) { break; }
          this.AddContact(this.proxyBuffer.data[a].index, this.proxyBuffer.data[b].index, this.contactBuffer);
        }
      }
    }

    ///void ReorderForFindContact(FindContactInput* reordered, int alignedCount) const;
    ///void GatherChecksOneParticle(const uint32 bound, const int startIndex, const int particleIndex, int* nextUncheckedIndex, GrowableBuffer<FindContactCheck>& checks) const;
    ///void GatherChecks(GrowableBuffer<FindContactCheck>& checks) const;
    ///void FindContacts_Simd(GrowableBuffer<ParticleContact>& contacts) const;

    public FindContacts(contacts: GrowableBuffer<ParticleContact>): void {
      this.FindContacts_Reference(contacts);
    }

    ///static void UpdateProxyTags(const uint32* const tags, GrowableBuffer<Proxy>& proxies);
    ///static bool ProxyBufferHasIndex(int32 index, const Proxy* const a, int count);
    ///static int NumProxiesWithSameTag(const Proxy* const a, const Proxy* const b, int count);
    ///static bool AreProxyBuffersTheSame(const GrowableBuffer<Proxy>& a, const GrowableBuffer<Proxy>& b);

    public UpdateProxies_Reference(proxies: GrowableBuffer<ParticleSysteProxy>): void {
      // DEBUG: Assert(proxies === this.proxyBuffer);
      const pos_data = this.positionBuffer.data;
      const inv_diam = this.inverseDiameter;
      for (let k = 0; k < this.proxyBuffer.count; ++k) {
        const proxy = this.proxyBuffer.data[k];
        const i = proxy.index;
        const p = pos_data[i];
        proxy.tag = ParticleSystem.computeTag(inv_diam * p.x, inv_diam * p.y);
      }
    }

    ///void UpdateProxies_Simd(GrowableBuffer<Proxy>& proxies) const;

    public UpdateProxies(proxies: GrowableBuffer<ParticleSysteProxy>): void {
      this.UpdateProxies_Reference(proxies);
    }

    public SortProxies(proxies: GrowableBuffer<ParticleSysteProxy>): void {
      // DEBUG: Assert(proxies === this.proxyBuffer);

      ///std::sort(proxies.Begin(), proxies.End());
      std_sort(this.proxyBuffer.data, 0, this.proxyBuffer.count, ParticleSysteProxy.CompareProxyProxy);
    }

    public FilterContacts(contacts: GrowableBuffer<ParticleContact>): void {
      // Optionally filter the contact.
      const contactFilter = this.GetParticleContactFilter();
      if (contactFilter === null) {
        return;
      }

      /// contacts.RemoveIf(ParticleContactRemovePredicate(this, contactFilter));
      // DEBUG: Assert(contacts === this.contactBuffer);
      const system = this;
      const predicate = (contact: ParticleContact): boolean => {
        return ((contact.flags & ParticleFlag.particleContactFilterParticle) !== 0) && !contactFilter.ShouldCollideParticleParticle(system, contact.indexA, contact.indexB);
      };
      this.contactBuffer.RemoveIf(predicate);
    }

    public NotifyContactListenerPreContact(particlePairs: ParticlePairSet): void {
      const contactListener = this.GetParticleContactListener();
      if (contactListener === null) {
        return;
      }

      ///particlePairs.Initialize(contactBuffer.Begin(), contactBuffer.GetCount(), GetFlagsBuffer());
      particlePairs.Initialize(this.contactBuffer, this.flagsBuffer);

      throw new Error(); // TODO: notify
    }

    public NotifyContactListenerPostContact(particlePairs: ParticlePairSet): void {
      const contactListener = this.GetParticleContactListener();
      if (contactListener === null) {
        return;
      }

      // Loop through all new contacts, reporting any new ones, and
      // "invalidating" the ones that still exist.
      ///const ParticleContact* const endContact = contactBuffer.End();
      ///for (ParticleContact* contact = contactBuffer.Begin(); contact < endContact; ++contact)
      for (let k = 0; k < this.contactBuffer.count; ++k) {
        const contact = this.contactBuffer.data[k];
        ///ParticlePair pair;
        ///pair.first = contact.GetIndexA();
        ///pair.second = contact.GetIndexB();
        ///const int32 itemIndex = particlePairs.Find(pair);
        const itemIndex = -1; // TODO
        if (itemIndex >= 0) {
          // Already touching, ignore this contact.
          particlePairs.Invalidate(itemIndex);
        } else {
          // Just started touching, inform the listener.
          contactListener.BeginContactParticleParticle(this, contact);
        }
      }

      // Report particles that are no longer touching.
      // That is, any pairs that were not invalidated above.
      ///const int32 pairCount = particlePairs.GetCount();
      ///const ParticlePair* const pairs = particlePairs.GetBuffer();
      ///const int8* const valid = particlePairs.GetValidBuffer();
      ///for (int32 i = 0; i < pairCount; ++i)
      ///{
      ///  if (valid[i])
      ///  {
      ///    contactListener.EndContactParticleParticle(this, pairs[i].first, pairs[i].second);
      ///  }
      ///}

      throw new Error(); // TODO: notify
    }

    public static ParticleContactIsZombie(contact: ParticleContact): boolean {
      return (contact.flags & ParticleFlag.zombieParticle) === ParticleFlag.zombieParticle;
    }

    public UpdateContacts(exceptZombie: boolean): void {
      this.UpdateProxies(this.proxyBuffer);
      this.SortProxies(this.proxyBuffer);

      const particlePairs = new ParticlePairSet(); // TODO: static
      this.NotifyContactListenerPreContact(particlePairs);

      this.FindContacts(this.contactBuffer);
      this.FilterContacts(this.contactBuffer);

      this.NotifyContactListenerPostContact(particlePairs);

      if (exceptZombie) {
        this.contactBuffer.RemoveIf(ParticleSystem.ParticleContactIsZombie);
      }
    }

    public NotifyBodyContactListenerPreContact(fixtureSet: ParticleSysteFixtureParticleSet): void {
      const contactListener = this.GetFixtureContactListener();
      if (contactListener === null) {
        return;
      }

      ///fixtureSet.Initialize(bodyContactBuffer.Begin(), bodyContactBuffer.GetCount(), GetFlagsBuffer());
      fixtureSet.Initialize(this.bodyContactBuffer, this.flagsBuffer);

      throw new Error(); // TODO: notify
    }

    public NotifyBodyContactListenerPostContact(fixtureSet: ParticleSysteFixtureParticleSet): void {
      const contactListener = this.GetFixtureContactListener();
      if (contactListener === null) {
        return;
      }

      // Loop through all new contacts, reporting any new ones, and
      // "invalidating" the ones that still exist.
      ///for (ParticleBodyContact* contact = bodyContactBuffer.Begin(); contact !== bodyContactBuffer.End(); ++contact)
      for (let k = 0; k < this.bodyContactBuffer.count; k++) {
        const contact = this.bodyContactBuffer.data[k];
        // DEBUG: Assert(contact !== null);
        ///FixtureParticle fixtureParticleToFind;
        ///fixtureParticleToFind.first = contact.fixture;
        ///fixtureParticleToFind.second = contact.index;
        ///const int32 index = fixtureSet.Find(fixtureParticleToFind);
        const index = -1; // TODO
        if (index >= 0) {
          // Already touching remove this from the set.
          fixtureSet.Invalidate(index);
        } else {
          // Just started touching, report it!
          contactListener.BeginContactFixtureParticle(this, contact);
        }
      }

      // If the contact listener is enabled, report all fixtures that are no
      // longer in contact with particles.
      ///const FixtureParticle* const fixtureParticles = fixtureSet.GetBuffer();
      ///const int8* const fixtureParticlesValid = fixtureSet.GetValidBuffer();
      ///const int32 fixtureParticleCount = fixtureSet.GetCount();
      ///for (int32 i = 0; i < fixtureParticleCount; ++i)
      ///{
      ///  if (fixtureParticlesValid[i])
      ///  {
      ///    const FixtureParticle* const fixtureParticle = &fixtureParticles[i];
      ///    contactListener.EndContactFixtureParticle(fixtureParticle.first, this, fixtureParticle.second);
      ///  }
      ///}

      throw new Error(); // TODO: notify
    }

    public UpdateBodyContacts(): void {
      const s_aabb = ParticleSystem.UpdateBodyContacts_s_aabb;

      // If the particle contact listener is enabled, generate a set of
      // fixture / particle contacts.
      const fixtureSet = new ParticleSysteFixtureParticleSet(); // TODO: static
      this.NotifyBodyContactListenerPreContact(fixtureSet);

      if (this.stuckThreshold > 0) {
        const particleCount = this.GetParticleCount();
        for (let i = 0; i < particleCount; i++) {
          // Detect stuck particles, see comment in
          // ParticleSystem::DetectStuckParticle()
          this.bodyContactCountBuffer.data[i] = 0;
          if (this.timestamp > (this.lastBodyContactStepBuffer.data[i] + 1)) {
            this.consecutiveContactStepsBuffer.data[i] = 0;
          }
        }
      }
      this.bodyContactBuffer.SetCount(0);
      this.stuckParticleBuffer.SetCount(0);

      const aabb = s_aabb;
      this.ComputeAABB(aabb);

      if (this.UpdateBodyContacts_callback === null) {
        this.UpdateBodyContacts_callback = new ParticleSysteUpdateBodyContactsCallback(this);
      }
      const callback = this.UpdateBodyContacts_callback;
      callback.contactFilter = this.GetFixtureContactFilter();
      this.world.QueryAABB(callback, aabb);

      if (this.def.strictContactCheck) {
        this.RemoveSpuriousBodyContacts();
      }

      this.NotifyBodyContactListenerPostContact(fixtureSet);
    }
    public static readonly UpdateBodyContacts_s_aabb = new AABB();
    public UpdateBodyContacts_callback: ParticleSysteUpdateBodyContactsCallback | null = null;

    public Solve(step: TimeStep): void {
      const s_subStep = ParticleSystem.Solve_s_subStep;
      if (this.count === 0) {
        return;
      }
      // If particle lifetimes are enabled, destroy particles that are too old.
      if (this.expirationTimeBuffer.data) {
        this.SolveLifetimes(step);
      }
      if (this.allParticleFlags & ParticleFlag.zombieParticle) {
        this.SolveZombie();
      }
      if (this.needsUpdateAllParticleFlags) {
        this.UpdateAllParticleFlags();
      }
      if (this.needsUpdateAllGroupFlags) {
        this.UpdateAllGroupFlags();
      }
      if (this.paused) {
        return;
      }
      for (this.iterationIndex = 0; this.iterationIndex < step.particleIterations; this.iterationIndex++) {
        ++this.timestamp;
        const subStep = s_subStep.Copy(step);
        subStep.dt /= step.particleIterations;
        subStep.inv_dt *= step.particleIterations;
        this.UpdateContacts(false);
        this.UpdateBodyContacts();
        this.ComputeWeight();
        if (this.allGroupFlags & ParticleGroupFlag.particleGroupNeedsUpdateDepth) {
          this.ComputeDepth();
        }
        if (this.allParticleFlags & ParticleFlag.reactiveParticle) {
          this.UpdatePairsAndTriadsWithReactiveParticles();
        }
        if (this.hasForce) {
          this.SolveForce(subStep);
        }
        if (this.allParticleFlags & ParticleFlag.viscousParticle) {
          this.SolveViscous();
        }
        if (this.allParticleFlags & ParticleFlag.repulsiveParticle) {
          this.SolveRepulsive(subStep);
        }
        if (this.allParticleFlags & ParticleFlag.powderParticle) {
          this.SolvePowder(subStep);
        }
        if (this.allParticleFlags & ParticleFlag.tensileParticle) {
          this.SolveTensile(subStep);
        }
        if (this.allGroupFlags & ParticleGroupFlag.solidParticleGroup) {
          this.SolveSolid(subStep);
        }
        if (this.allParticleFlags & ParticleFlag.colorMixingParticle) {
          this.SolveColorMixing();
        }
        this.SolveGravity(subStep);
        if (this.allParticleFlags & ParticleFlag.staticPressureParticle) {
          this.SolveStaticPressure(subStep);
        }
        this.SolvePressure(subStep);
        this.SolveDamping(subStep);
        if (this.allParticleFlags & ParticleSystem.k_extraDampingFlags) {
          this.SolveExtraDamping();
        }
        // SolveElastic and SolveSpring refer the current velocities for
        // numerical stability, they should be called as late as possible.
        if (this.allParticleFlags & ParticleFlag.elasticParticle) {
          this.SolveElastic(subStep);
        }
        if (this.allParticleFlags & ParticleFlag.springParticle) {
          this.SolveSpring(subStep);
        }
        this.LimitVelocity(subStep);
        if (this.allGroupFlags & ParticleGroupFlag.rigidParticleGroup) {
          this.SolveRigidDamping();
        }
        if (this.allParticleFlags & ParticleFlag.barrierParticle) {
          this.SolveBarrier(subStep);
        }
        // SolveCollision, SolveRigid and SolveWall should be called after
        // other force functions because they may require particles to have
        // specific velocities.
        this.SolveCollision(subStep);
        if (this.allGroupFlags & ParticleGroupFlag.rigidParticleGroup) {
          this.SolveRigid(subStep);
        }
        if (this.allParticleFlags & ParticleFlag.wallParticle) {
          this.SolveWall();
        }
        // The particle positions can be updated only at the end of substep.
        for (let i = 0; i < this.count; i++) {
          ///positionBuffer.data[i] += subStep.dt * velocityBuffer.data[i];
          this.positionBuffer.data[i].SelfMulAdd(subStep.dt, this.velocityBuffer.data[i]);
        }
      }
    }
    public static readonly Solve_s_subStep = new TimeStep();

    public SolveCollision(step: TimeStep): void {
      const s_aabb = ParticleSystem.SolveCollision_s_aabb;
      const pos_data = this.positionBuffer.data;
      const vel_data = this.velocityBuffer.data;

      // This function detects particles which are crossing boundary of bodies
      // and modifies velocities of them so that they will move just in front of
      // boundary. This function function also applies the reaction force to
      // bodies as precisely as the numerical stability is kept.
      const aabb = s_aabb;
      aabb.lowerBound.x = +maxFloat;
      aabb.lowerBound.y = +maxFloat;
      aabb.upperBound.x = -maxFloat;
      aabb.upperBound.y = -maxFloat;
      for (let i = 0; i < this.count; i++) {
        const v = vel_data[i];
        const p1 = pos_data[i];
        ///let p2 = p1 + step.dt * v;
        const p2_x = p1.x + step.dt * v.x;
        const p2_y = p1.y + step.dt * v.y;
        ///aabb.lowerBound = Min(aabb.lowerBound, Min(p1, p2));
        aabb.lowerBound.x = Min(aabb.lowerBound.x, Min(p1.x, p2_x));
        aabb.lowerBound.y = Min(aabb.lowerBound.y, Min(p1.y, p2_y));
        ///aabb.upperBound = Max(aabb.upperBound, Max(p1, p2));
        aabb.upperBound.x = Max(aabb.upperBound.x, Max(p1.x, p2_x));
        aabb.upperBound.y = Max(aabb.upperBound.y, Max(p1.y, p2_y));
      }
      if (this.SolveCollision_callback === null) {
        this.SolveCollision_callback = new ParticleSysteSolveCollisionCallback(this, step);
      }
      const callback = this.SolveCollision_callback;
      callback.step = step;
      this.world.QueryAABB(callback, aabb);
    }
    public static readonly SolveCollision_s_aabb = new AABB();
    public SolveCollision_callback: ParticleSysteSolveCollisionCallback | null = null;

    public LimitVelocity(step: TimeStep): void {
      const vel_data = this.velocityBuffer.data;
      const criticalVelocitySquared = this.GetCriticalVelocitySquared(step);
      for (let i = 0; i < this.count; i++) {
        const v = vel_data[i];
        const v2 = Vec2.DotVV(v, v);
        if (v2 > criticalVelocitySquared) {
          ///v *= Sqrt(criticalVelocitySquared / v2);
          v.SelfMul(Sqrt(criticalVelocitySquared / v2));
        }
      }
    }

    public SolveGravity(step: TimeStep): void {
      const s_gravity = ParticleSystem.SolveGravity_s_gravity;
      const vel_data = this.velocityBuffer.data;
      ///Vec2 gravity = step.dt * def.gravityScale * world.GetGravity();
      const gravity = Vec2.MulSV(step.dt * this.def.gravityScale, this.world.GetGravity(), s_gravity);
      for (let i = 0; i < this.count; i++) {
        vel_data[i].SelfAdd(gravity);
      }
    }
    public static readonly SolveGravity_s_gravity = new Vec2();

    public SolveBarrier(step: TimeStep): void {
      const s_aabb = ParticleSystem.SolveBarrier_s_aabb;
      const s_va = ParticleSystem.SolveBarrier_s_va;
      const s_vb = ParticleSystem.SolveBarrier_s_vb;
      const s_pba = ParticleSystem.SolveBarrier_s_pba;
      const s_vba = ParticleSystem.SolveBarrier_s_vba;
      const s_vc = ParticleSystem.SolveBarrier_s_vc;
      const s_pca = ParticleSystem.SolveBarrier_s_pca;
      const s_vca = ParticleSystem.SolveBarrier_s_vca;
      const s_qba = ParticleSystem.SolveBarrier_s_qba;
      const s_qca = ParticleSystem.SolveBarrier_s_qca;
      const s_dv = ParticleSystem.SolveBarrier_s_dv;
      const s_f = ParticleSystem.SolveBarrier_s_f;
      const pos_data = this.positionBuffer.data;
      const vel_data = this.velocityBuffer.data;
      // If a particle is passing between paired barrier particles,
      // its velocity will be decelerated to avoid passing.
      for (let i = 0; i < this.count; i++) {
        const flags = this.flagsBuffer.data[i];
        ///if ((flags & ParticleSystem.k_barrierWallFlags) === ParticleSystem.k_barrierWallFlags)
        if ((flags & ParticleSystem.k_barrierWallFlags) !== 0) {
          vel_data[i].SetZero();
        }
      }
      const tmax = barrierCollisionTime * step.dt;
      const mass = this.GetParticleMass();
      for (let k = 0; k < this.pairBuffer.count; k++) {
        const pair = this.pairBuffer.data[k];
        if (pair.flags & ParticleFlag.barrierParticle) {
          const a = pair.indexA;
          const b = pair.indexB;
          const pa = pos_data[a];
          const pb = pos_data[b];
          /// AABB aabb;
          const aabb = s_aabb;
          ///aabb.lowerBound = Min(pa, pb);
          Vec2.MinV(pa, pb, aabb.lowerBound);
          ///aabb.upperBound = Max(pa, pb);
          Vec2.MaxV(pa, pb, aabb.upperBound);
          const aGroup = this.groupBuffer[a];
          const bGroup = this.groupBuffer[b];
          ///Vec2 va = GetLinearVelocity(aGroup, a, pa);
          const va = this.GetLinearVelocity(aGroup, a, pa, s_va);
          ///Vec2 vb = GetLinearVelocity(bGroup, b, pb);
          const vb = this.GetLinearVelocity(bGroup, b, pb, s_vb);
          ///Vec2 pba = pb - pa;
          const pba = Vec2.SubVV(pb, pa, s_pba);
          ///Vec2 vba = vb - va;
          const vba = Vec2.SubVV(vb, va, s_vba);
          ///InsideBoundsEnumerator enumerator = GetInsideBoundsEnumerator(aabb);
          const enumerator = this.GetInsideBoundsEnumerator(aabb);
          let c: number;
          while ((c = enumerator.GetNext()) >= 0) {
            const pc = pos_data[c];
            const cGroup = this.groupBuffer[c];
            if (aGroup !== cGroup && bGroup !== cGroup) {
              ///Vec2 vc = GetLinearVelocity(cGroup, c, pc);
              const vc = this.GetLinearVelocity(cGroup, c, pc, s_vc);
              // Solve the equation below:
              //   (1-s)*(pa+t*va)+s*(pb+t*vb) = pc+t*vc
              // which expresses that the particle c will pass a line
              // connecting the particles a and b at the time of t.
              // if s is between 0 and 1, c will pass between a and b.
              ///Vec2 pca = pc - pa;
              const pca = Vec2.SubVV(pc, pa, s_pca);
              ///Vec2 vca = vc - va;
              const vca = Vec2.SubVV(vc, va, s_vca);
              const e2 = Vec2.CrossVV(vba, vca);
              const e1 = Vec2.CrossVV(pba, vca) - Vec2.CrossVV(pca, vba);
              const e0 = Vec2.CrossVV(pba, pca);
              let s: number, t: number;
              ///Vec2 qba, qca;
              const qba = s_qba,
                qca = s_qca;
              if (e2 === 0) {
                if (e1 === 0) { continue; }
                t = -e0 / e1;
                if (!(t >= 0 && t < tmax)) { continue; }
                ///qba = pba + t * vba;
                Vec2.AddVMulSV(pba, t, vba, qba);
                ///qca = pca + t * vca;
                Vec2.AddVMulSV(pca, t, vca, qca);
                s = Vec2.DotVV(qba, qca) / Vec2.DotVV(qba, qba);
                if (!(s >= 0 && s <= 1)) { continue; }
              } else {
                const det = e1 * e1 - 4 * e0 * e2;
                if (det < 0) { continue; }
                const sqrtDet = Sqrt(det);
                let t1 = (-e1 - sqrtDet) / (2 * e2);
                let t2 = (-e1 + sqrtDet) / (2 * e2);
                ///if (t1 > t2) Swap(t1, t2);
                if (t1 > t2) {
                  const tmp = t1;
                  t1 = t2;
                  t2 = tmp;
                }
                t = t1;
                ///qba = pba + t * vba;
                Vec2.AddVMulSV(pba, t, vba, qba);
                ///qca = pca + t * vca;
                Vec2.AddVMulSV(pca, t, vca, qca);
                ///s = Dot(qba, qca) / Dot(qba, qba);
                s = Vec2.DotVV(qba, qca) / Vec2.DotVV(qba, qba);
                if (!(t >= 0 && t < tmax && s >= 0 && s <= 1)) {
                  t = t2;
                  if (!(t >= 0 && t < tmax)) { continue; }
                  ///qba = pba + t * vba;
                  Vec2.AddVMulSV(pba, t, vba, qba);
                  ///qca = pca + t * vca;
                  Vec2.AddVMulSV(pca, t, vca, qca);
                  ///s = Dot(qba, qca) / Dot(qba, qba);
                  s = Vec2.DotVV(qba, qca) / Vec2.DotVV(qba, qba);
                  if (!(s >= 0 && s <= 1)) { continue; }
                }
              }
              // Apply a force to particle c so that it will have the
              // interpolated velocity at the collision point on line ab.
              ///Vec2 dv = va + s * vba - vc;
              const dv = s_dv;
              dv.x = va.x + s * vba.x - vc.x;
              dv.y = va.y + s * vba.y - vc.y;
              ///Vec2 f = GetParticleMass() * dv;
              const f = Vec2.MulSV(mass, dv, s_f);
              if (cGroup && this.IsRigidGroup(cGroup)) {
                // If c belongs to a rigid group, the force will be
                // distributed in the group.
                const mass = cGroup.GetMass();
                const inertia = cGroup.GetInertia();
                if (mass > 0) {
                  ///cGroup.linearVelocity += 1 / mass * f;
                  cGroup.linearVelocity.SelfMulAdd(1 / mass, f);
                }
                if (inertia > 0) {
                  ///cGroup.angularVelocity += Cross(pc - cGroup.GetCenter(), f) / inertia;
                  cGroup.angularVelocity += Vec2.CrossVV(
                    Vec2.SubVV(pc, cGroup.GetCenter(), Vec2.s_t0),
                    f) / inertia;
                }
              } else {
                ///velocityBuffer.data[c] += dv;
                vel_data[c].SelfAdd(dv);
              }
              // Apply a reversed force to particle c after particle
              // movement so that momentum will be preserved.
              ///ParticleApplyForce(c, -step.inv_dt * f);
              this.ParticleApplyForce(c, f.SelfMul(-step.inv_dt));
            }
          }
        }
      }
    }
    public static readonly SolveBarrier_s_aabb = new AABB();
    public static readonly SolveBarrier_s_va = new Vec2();
    public static readonly SolveBarrier_s_vb = new Vec2();
    public static readonly SolveBarrier_s_pba = new Vec2();
    public static readonly SolveBarrier_s_vba = new Vec2();
    public static readonly SolveBarrier_s_vc = new Vec2();
    public static readonly SolveBarrier_s_pca = new Vec2();
    public static readonly SolveBarrier_s_vca = new Vec2();
    public static readonly SolveBarrier_s_qba = new Vec2();
    public static readonly SolveBarrier_s_qca = new Vec2();
    public static readonly SolveBarrier_s_dv = new Vec2();
    public static readonly SolveBarrier_s_f = new Vec2();

    public SolveStaticPressure(step: TimeStep): void {
      this.staticPressureBuffer = this.RequestBuffer(this.staticPressureBuffer);
      const criticalPressure = this.GetCriticalPressure(step);
      const pressurePerWeight = this.def.staticPressureStrength * criticalPressure;
      const maxPressure = maxParticlePressure * criticalPressure;
      const relaxation = this.def.staticPressureRelaxation;
      /// Compute pressure satisfying the modified Poisson equation:
      ///   Sufor_j((p_i - p_j) * w_ij) + relaxation * p_i =
      ///   pressurePerWeight * (w_i - minParticleWeight)
      /// by iterating the calculation:
      ///   p_i = (Sufor_j(p_j * w_ij) + pressurePerWeight *
      ///         (w_i - minParticleWeight)) / (w_i + relaxation)
      /// where
      ///   p_i and p_j are static pressure of particle i and j
      ///   w_ij is contact weight between particle i and j
      ///   w_i is sum of contact weight of particle i
      for (let t = 0; t < this.def.staticPressureIterations; t++) {
        ///memset(accumulationBuffer, 0, sizeof(*accumulationBuffer) * count);
        for (let i = 0; i < this.count; i++) {
          this.accumulationBuffer[i] = 0;
        }
        for (let k = 0; k < this.contactBuffer.count; k++) {
          const contact = this.contactBuffer.data[k];
          if (contact.flags & ParticleFlag.staticPressureParticle) {
            const a = contact.indexA;
            const b = contact.indexB;
            const w = contact.weight;
            this.accumulationBuffer[a] += w * this.staticPressureBuffer[b]; // a <- b
            this.accumulationBuffer[b] += w * this.staticPressureBuffer[a]; // b <- a
          }
        }
        for (let i = 0; i < this.count; i++) {
          const w = this.weightBuffer[i];
          if (this.flagsBuffer.data[i] & ParticleFlag.staticPressureParticle) {
            const wh = this.accumulationBuffer[i];
            const h =
              (wh + pressurePerWeight * (w - minParticleWeight)) /
              (w + relaxation);
            this.staticPressureBuffer[i] = Clamp(h, 0.0, maxPressure);
          } else {
            this.staticPressureBuffer[i] = 0;
          }
        }
      }
    }

    public ComputeWeight(): void {
      // calculates the sum of contact-weights for each particle
      // that means dimensionless density
      ///memset(weightBuffer, 0, sizeof(*weightBuffer) * count);
      for (let k = 0; k < this.count; k++) {
        this.weightBuffer[k] = 0;
      }
      for (let k = 0; k < this.bodyContactBuffer.count; k++) {
        const contact = this.bodyContactBuffer.data[k];
        const a = contact.index;
        const w = contact.weight;
        this.weightBuffer[a] += w;
      }
      for (let k = 0; k < this.contactBuffer.count; k++) {
        const contact = this.contactBuffer.data[k];
        const a = contact.indexA;
        const b = contact.indexB;
        const w = contact.weight;
        this.weightBuffer[a] += w;
        this.weightBuffer[b] += w;
      }
    }

    public SolvePressure(step: TimeStep): void {
      const s_f = ParticleSystem.SolvePressure_s_f;
      const pos_data = this.positionBuffer.data;
      const vel_data = this.velocityBuffer.data;
      // calculates pressure as a linear function of density
      const criticalPressure = this.GetCriticalPressure(step);
      const pressurePerWeight = this.def.pressureStrength * criticalPressure;
      const maxPressure = maxParticlePressure * criticalPressure;
      for (let i = 0; i < this.count; i++) {
        const w = this.weightBuffer[i];
        const h = pressurePerWeight * Max(0.0, w - minParticleWeight);
        this.accumulationBuffer[i] = Min(h, maxPressure);
      }
      // ignores particles which have their own repulsive force
      if (this.allParticleFlags & ParticleSystem.k_noPressureFlags) {
        for (let i = 0; i < this.count; i++) {
          if (this.flagsBuffer.data[i] & ParticleSystem.k_noPressureFlags) {
            this.accumulationBuffer[i] = 0;
          }
        }
      }
      // static pressure
      if (this.allParticleFlags & ParticleFlag.staticPressureParticle) {
        // DEBUG: Assert(this.staticPressureBuffer !== null);
        for (let i = 0; i < this.count; i++) {
          if (this.flagsBuffer.data[i] & ParticleFlag.staticPressureParticle) {
            this.accumulationBuffer[i] += this.staticPressureBuffer[i];
          }
        }
      }
      // applies pressure between each particles in contact
      const velocityPerPressure = step.dt / (this.def.density * this.particleDiameter);
      const inv_mass = this.GetParticleInvMass();
      for (let k = 0; k < this.bodyContactBuffer.count; k++) {
        const contact = this.bodyContactBuffer.data[k];
        const a = contact.index;
        const b = contact.body;
        const w = contact.weight;
        const m = contact.mass;
        const n = contact.normal;
        const p = pos_data[a];
        const h = this.accumulationBuffer[a] + pressurePerWeight * w;
        ///Vec2 f = velocityPerPressure * w * m * h * n;
        const f = Vec2.MulSV(velocityPerPressure * w * m * h, n, s_f);
        ///velocityBuffer.data[a] -= GetParticleInvMass() * f;
        vel_data[a].SelfMulSub(inv_mass, f);
        b.ApplyLinearImpulse(f, p, true);
      }
      for (let k = 0; k < this.contactBuffer.count; k++) {
        const contact = this.contactBuffer.data[k];
        const a = contact.indexA;
        const b = contact.indexB;
        const w = contact.weight;
        const n = contact.normal;
        const h = this.accumulationBuffer[a] + this.accumulationBuffer[b];
        ///Vec2 f = velocityPerPressure * w * h * n;
        const f = Vec2.MulSV(velocityPerPressure * w * h, n, s_f);
        ///velocityBuffer.data[a] -= f;
        vel_data[a].SelfSub(f);
        ///velocityBuffer.data[b] += f;
        vel_data[b].SelfAdd(f);
      }
    }
    public static readonly SolvePressure_s_f = new Vec2();

    public SolveDamping(step: TimeStep): void {
      const s_v = ParticleSystem.SolveDamping_s_v;
      const s_f = ParticleSystem.SolveDamping_s_f;
      const pos_data = this.positionBuffer.data;
      const vel_data = this.velocityBuffer.data;
      // reduces normal velocity of each contact
      const linearDamping = this.def.dampingStrength;
      const quadraticDamping = 1 / this.GetCriticalVelocity(step);
      const inv_mass = this.GetParticleInvMass();
      for (let k = 0; k < this.bodyContactBuffer.count; k++) {
        const contact = this.bodyContactBuffer.data[k];
        const a = contact.index;
        const b = contact.body;
        const w = contact.weight;
        const m = contact.mass;
        const n = contact.normal;
        const p = pos_data[a];
        ///Vec2 v = b.GetLinearVelocityFromWorldPoint(p) - velocityBuffer.data[a];
        const v = Vec2.SubVV(b.GetLinearVelocityFromWorldPoint(p, Vec2.s_t0), vel_data[a], s_v);
        const vn = Vec2.DotVV(v, n);
        if (vn < 0) {
          const damping = Max(linearDamping * w, Min(-quadraticDamping * vn, 0.5));
          ///Vec2 f = damping * m * vn * n;
          const f = Vec2.MulSV(damping * m * vn, n, s_f);
          ///velocityBuffer.data[a] += GetParticleInvMass() * f;
          vel_data[a].SelfMulAdd(inv_mass, f);
          ///b.ApplyLinearImpulse(-f, p, true);
          b.ApplyLinearImpulse(f.SelfNeg(), p, true);
        }
      }
      for (let k = 0; k < this.contactBuffer.count; k++) {
        const contact = this.contactBuffer.data[k];
        const a = contact.indexA;
        const b = contact.indexB;
        const w = contact.weight;
        const n = contact.normal;
        ///Vec2 v = velocityBuffer.data[b] - velocityBuffer.data[a];
        const v = Vec2.SubVV(vel_data[b], vel_data[a], s_v);
        const vn = Vec2.DotVV(v, n);
        if (vn < 0) {
          ///float32 damping = Max(linearDamping * w, Min(- quadraticDamping * vn, 0.5f));
          const damping = Max(linearDamping * w, Min(-quadraticDamping * vn, 0.5));
          ///Vec2 f = damping * vn * n;
          const f = Vec2.MulSV(damping * vn, n, s_f);
          ///this.velocityBuffer.data[a] += f;
          vel_data[a].SelfAdd(f);
          ///this.velocityBuffer.data[b] -= f;
          vel_data[b].SelfSub(f);
        }
      }
    }
    public static readonly SolveDamping_s_v = new Vec2();
    public static readonly SolveDamping_s_f = new Vec2();

    public SolveRigidDamping(): void {
      const s_t0 = ParticleSystem.SolveRigidDamping_s_t0;
      const s_t1 = ParticleSystem.SolveRigidDamping_s_t1;
      const s_p = ParticleSystem.SolveRigidDamping_s_p;
      const s_v = ParticleSystem.SolveRigidDamping_s_v;
      const invMassA = [0.0],
        invInertiaA = [0.0],
        tangentDistanceA = [0.0]; // TODO: static
      const invMassB = [0.0],
        invInertiaB = [0.0],
        tangentDistanceB = [0.0]; // TODO: static
      // Apply impulse to rigid particle groups colliding with other objects
      // to reduce relative velocity at the colliding point.
      const pos_data = this.positionBuffer.data;
      const damping = this.def.dampingStrength;
      for (let k = 0; k < this.bodyContactBuffer.count; k++) {
        const contact = this.bodyContactBuffer.data[k];
        const a = contact.index;
        const aGroup = this.groupBuffer[a];
        if (aGroup && this.IsRigidGroup(aGroup)) {
          const b = contact.body;
          const n = contact.normal;
          const w = contact.weight;
          const p = pos_data[a];
          ///Vec2 v = b.GetLinearVelocityFromWorldPoint(p) - aGroup.GetLinearVelocityFromWorldPoint(p);
          const v = Vec2.SubVV(b.GetLinearVelocityFromWorldPoint(p, s_t0), aGroup.GetLinearVelocityFromWorldPoint(p, s_t1), s_v);
          const vn = Vec2.DotVV(v, n);
          if (vn < 0) {
            // The group's average velocity at particle position 'p' is pushing
            // the particle into the body.
            ///this.InitDampingParameterWithRigidGroupOrParticle(&invMassA, &invInertiaA, &tangentDistanceA, true, aGroup, a, p, n);
            this.InitDampingParameterWithRigidGroupOrParticle(invMassA, invInertiaA, tangentDistanceA, true, aGroup, a, p, n);
            // Calculate b.I from public functions of Body.
            ///this.InitDampingParameter(&invMassB, &invInertiaB, &tangentDistanceB, b.GetMass(), b.GetInertia() - b.GetMass() * b.GetLocalCenter().LengthSquared(), b.GetWorldCenter(), p, n);
            this.InitDampingParameter(invMassB, invInertiaB, tangentDistanceB, b.GetMass(), b.GetInertia() - b.GetMass() * b.GetLocalCenter().LengthSquared(), b.GetWorldCenter(), p, n);
            ///float32 f = damping * Min(w, 1.0) * this.ComputeDampingImpulse(invMassA, invInertiaA, tangentDistanceA, invMassB, invInertiaB, tangentDistanceB, vn);
            const f = damping * Min(w, 1.0) * this.ComputeDampingImpulse(invMassA[0], invInertiaA[0], tangentDistanceA[0], invMassB[0], invInertiaB[0], tangentDistanceB[0], vn);
            ///this.ApplyDamping(invMassA, invInertiaA, tangentDistanceA, true, aGroup, a, f, n);
            this.ApplyDamping(invMassA[0], invInertiaA[0], tangentDistanceA[0], true, aGroup, a, f, n);
            ///b.ApplyLinearImpulse(-f * n, p, true);
            b.ApplyLinearImpulse(Vec2.MulSV(-f, n, Vec2.s_t0), p, true);
          }
        }
      }
      for (let k = 0; k < this.contactBuffer.count; k++) {
        const contact = this.contactBuffer.data[k];
        const a = contact.indexA;
        const b = contact.indexB;
        const n = contact.normal;
        const w = contact.weight;
        const aGroup = this.groupBuffer[a];
        const bGroup = this.groupBuffer[b];
        const aRigid = this.IsRigidGroup(aGroup);
        const bRigid = this.IsRigidGroup(bGroup);
        if (aGroup !== bGroup && (aRigid || bRigid)) {
          ///Vec2 p = 0.5f * (this.positionBuffer.data[a] + this.positionBuffer.data[b]);
          const p = Vec2.MidVV(pos_data[a], pos_data[b], s_p);
          ///Vec2 v = GetLinearVelocity(bGroup, b, p) - GetLinearVelocity(aGroup, a, p);
          const v = Vec2.SubVV(this.GetLinearVelocity(bGroup, b, p, s_t0), this.GetLinearVelocity(aGroup, a, p, s_t1), s_v);
          const vn = Vec2.DotVV(v, n);
          if (vn < 0) {
            ///this.InitDampingParameterWithRigidGroupOrParticle(&invMassA, &invInertiaA, &tangentDistanceA, aRigid, aGroup, a, p, n);
            this.InitDampingParameterWithRigidGroupOrParticle(invMassA, invInertiaA, tangentDistanceA, aRigid, aGroup, a, p, n);
            ///this.InitDampingParameterWithRigidGroupOrParticle(&invMassB, &invInertiaB, &tangentDistanceB, bRigid, bGroup, b, p, n);
            this.InitDampingParameterWithRigidGroupOrParticle(invMassB, invInertiaB, tangentDistanceB, bRigid, bGroup, b, p, n);
            ///float32 f = damping * w * this.ComputeDampingImpulse(invMassA, invInertiaA, tangentDistanceA, invMassB, invInertiaB, tangentDistanceB, vn);
            const f = damping * w * this.ComputeDampingImpulse(invMassA[0], invInertiaA[0], tangentDistanceA[0], invMassB[0], invInertiaB[0], tangentDistanceB[0], vn);
            ///this.ApplyDamping(invMassA, invInertiaA, tangentDistanceA, aRigid, aGroup, a, f, n);
            this.ApplyDamping(invMassA[0], invInertiaA[0], tangentDistanceA[0], aRigid, aGroup, a, f, n);
            ///this.ApplyDamping(invMassB, invInertiaB, tangentDistanceB, bRigid, bGroup, b, -f, n);
            this.ApplyDamping(invMassB[0], invInertiaB[0], tangentDistanceB[0], bRigid, bGroup, b, -f, n);
          }
        }
      }
    }
    public static readonly SolveRigidDamping_s_t0 = new Vec2();
    public static readonly SolveRigidDamping_s_t1 = new Vec2();
    public static readonly SolveRigidDamping_s_p = new Vec2();
    public static readonly SolveRigidDamping_s_v = new Vec2();

    public SolveExtraDamping(): void {
      const s_v = ParticleSystem.SolveExtraDamping_s_v;
      const s_f = ParticleSystem.SolveExtraDamping_s_f;
      const vel_data = this.velocityBuffer.data;
      // Applies additional damping force between bodies and particles which can
      // produce strong repulsive force. Applying damping force multiple times
      // is effective in suppressing vibration.
      const pos_data = this.positionBuffer.data;
      const inv_mass = this.GetParticleInvMass();
      for (let k = 0; k < this.bodyContactBuffer.count; k++) {
        const contact = this.bodyContactBuffer.data[k];
        const a = contact.index;
        if (this.flagsBuffer.data[a] & ParticleSystem.k_extraDampingFlags) {
          const b = contact.body;
          const m = contact.mass;
          const n = contact.normal;
          const p = pos_data[a];
          ///Vec2 v = b.GetLinearVelocityFromWorldPoint(p) - velocityBuffer.data[a];
          const v = Vec2.SubVV(b.GetLinearVelocityFromWorldPoint(p, Vec2.s_t0), vel_data[a], s_v);
          ///float32 vn = Dot(v, n);
          const vn = Vec2.DotVV(v, n);
          if (vn < 0) {
            ///Vec2 f = 0.5f * m * vn * n;
            const f = Vec2.MulSV(0.5 * m * vn, n, s_f);
            ///velocityBuffer.data[a] += GetParticleInvMass() * f;
            vel_data[a].SelfMulAdd(inv_mass, f);
            ///b.ApplyLinearImpulse(-f, p, true);
            b.ApplyLinearImpulse(f.SelfNeg(), p, true);
          }
        }
      }
    }
    public static readonly SolveExtraDamping_s_v = new Vec2();
    public static readonly SolveExtraDamping_s_f = new Vec2();

    public SolveWall(): void {
      const vel_data = this.velocityBuffer.data;
      for (let i = 0; i < this.count; i++) {
        if (this.flagsBuffer.data[i] & ParticleFlag.wallParticle) {
          vel_data[i].SetZero();
        }
      }
    }

    public SolveRigid(step: TimeStep): void {
      const s_position = ParticleSystem.SolveRigid_s_position;
      const s_rotation = ParticleSystem.SolveRigid_s_rotation;
      const s_transform = ParticleSystem.SolveRigid_s_transform;
      const s_velocityTransform = ParticleSystem.SolveRigid_s_velocityTransform;
      const pos_data = this.positionBuffer.data;
      const vel_data = this.velocityBuffer.data;
      for (let group = this.groupList; group; group = group.GetNext()) {
        if (group.groupFlags & ParticleGroupFlag.rigidParticleGroup) {
          group.UpdateStatistics();
          ///Rot rotation(step.dt * group.angularVelocity);
          const rotation = s_rotation;
          rotation.SetAngle(step.dt * group.angularVelocity);
          ///Transform transform(group.center + step.dt * group.linearVelocity - Mul(rotation, group.center), rotation);
          const position = Vec2.AddVV(
            group.center,
            Vec2.SubVV(
              Vec2.MulSV(step.dt, group.linearVelocity, Vec2.s_t0),
              Rot.MulRV(rotation, group.center, Vec2.s_t1),
              Vec2.s_t0),
            s_position);
          const transform = s_transform;
          transform.SetPositionRotation(position, rotation);
          ///group.transform = Mul(transform, group.transform);
          Transform.MulXX(transform, group.transform, group.transform);
          const velocityTransform = s_velocityTransform;
          velocityTransform.p.x = step.inv_dt * transform.p.x;
          velocityTransform.p.y = step.inv_dt * transform.p.y;
          velocityTransform.q.s = step.inv_dt * transform.q.s;
          velocityTransform.q.c = step.inv_dt * (transform.q.c - 1);
          for (let i = group.firstIndex; i < group.lastIndex; i++) {
            ///velocityBuffer.data[i] = Mul(velocityTransform, positionBuffer.data[i]);
            Transform.MulXV(velocityTransform, pos_data[i], vel_data[i]);
          }
        }
      }
    }
    public static readonly SolveRigid_s_position = new Vec2();
    public static readonly SolveRigid_s_rotation = new Rot();
    public static readonly SolveRigid_s_transform = new Transform();
    public static readonly SolveRigid_s_velocityTransform = new Transform();

    public SolveElastic(step: TimeStep): void {
      const s_pa = ParticleSystem.SolveElastic_s_pa;
      const s_pb = ParticleSystem.SolveElastic_s_pb;
      const s_pc = ParticleSystem.SolveElastic_s_pc;
      const s_r = ParticleSystem.SolveElastic_s_r;
      const s_t0 = ParticleSystem.SolveElastic_s_t0;
      const pos_data = this.positionBuffer.data;
      const vel_data = this.velocityBuffer.data;
      const elasticStrength = step.inv_dt * this.def.elasticStrength;
      for (let k = 0; k < this.triadBuffer.count; k++) {
        const triad = this.triadBuffer.data[k];
        if (triad.flags & ParticleFlag.elasticParticle) {
          const a = triad.indexA;
          const b = triad.indexB;
          const c = triad.indexC;
          const oa = triad.pa;
          const ob = triad.pb;
          const oc = triad.pc;
          ///Vec2 pa = positionBuffer.data[a];
          const pa = s_pa.Copy(pos_data[a]);
          ///Vec2 pb = positionBuffer.data[b];
          const pb = s_pb.Copy(pos_data[b]);
          ///Vec2 pc = positionBuffer.data[c];
          const pc = s_pc.Copy(pos_data[c]);
          const va = vel_data[a];
          const vb = vel_data[b];
          const vc = vel_data[c];
          ///pa += step.dt * va;
          pa.SelfMulAdd(step.dt, va);
          ///pb += step.dt * vb;
          pb.SelfMulAdd(step.dt, vb);
          ///pc += step.dt * vc;
          pc.SelfMulAdd(step.dt, vc);
          ///Vec2 midPoint = (float32) 1 / 3 * (pa + pb + pc);
          const midPoint_x = (pa.x + pb.x + pc.x) / 3.0;
          const midPoint_y = (pa.y + pb.y + pc.y) / 3.0;
          ///pa -= midPoint;
          pa.x -= midPoint_x;
          pa.y -= midPoint_y;
          ///pb -= midPoint;
          pb.x -= midPoint_x;
          pb.y -= midPoint_y;
          ///pc -= midPoint;
          pc.x -= midPoint_x;
          pc.y -= midPoint_y;
          ///Rot r;
          const r = s_r;
          r.s = Vec2.CrossVV(oa, pa) + Vec2.CrossVV(ob, pb) + Vec2.CrossVV(oc, pc);
          r.c = Vec2.DotVV(oa, pa) + Vec2.DotVV(ob, pb) + Vec2.DotVV(oc, pc);
          const r2 = r.s * r.s + r.c * r.c;
          let invR = InvSqrt(r2);
          if (!isFinite(invR)) {
            invR = 1.98177537e+019;
          }
          r.s *= invR;
          r.c *= invR;
          ///r.angle = Math.atan2(r.s, r.c); // TODO: optimize
          const strength = elasticStrength * triad.strength;
          ///va += strength * (Mul(r, oa) - pa);
          Rot.MulRV(r, oa, s_t0);
          Vec2.SubVV(s_t0, pa, s_t0);
          Vec2.MulSV(strength, s_t0, s_t0);
          va.SelfAdd(s_t0);
          ///vb += strength * (Mul(r, ob) - pb);
          Rot.MulRV(r, ob, s_t0);
          Vec2.SubVV(s_t0, pb, s_t0);
          Vec2.MulSV(strength, s_t0, s_t0);
          vb.SelfAdd(s_t0);
          ///vc += strength * (Mul(r, oc) - pc);
          Rot.MulRV(r, oc, s_t0);
          Vec2.SubVV(s_t0, pc, s_t0);
          Vec2.MulSV(strength, s_t0, s_t0);
          vc.SelfAdd(s_t0);
        }
      }
    }
    public static readonly SolveElastic_s_pa = new Vec2();
    public static readonly SolveElastic_s_pb = new Vec2();
    public static readonly SolveElastic_s_pc = new Vec2();
    public static readonly SolveElastic_s_r = new Rot();
    public static readonly SolveElastic_s_t0 = new Vec2();

    public SolveSpring(step: TimeStep): void {
      const s_pa = ParticleSystem.SolveSpring_s_pa;
      const s_pb = ParticleSystem.SolveSpring_s_pb;
      const s_d = ParticleSystem.SolveSpring_s_d;
      const s_f = ParticleSystem.SolveSpring_s_f;
      const pos_data = this.positionBuffer.data;
      const vel_data = this.velocityBuffer.data;
      const springStrength = step.inv_dt * this.def.springStrength;
      for (let k = 0; k < this.pairBuffer.count; k++) {
        const pair = this.pairBuffer.data[k];
        if (pair.flags & ParticleFlag.springParticle) {
          ///int32 a = pair.indexA;
          const a = pair.indexA;
          ///int32 b = pair.indexB;
          const b = pair.indexB;
          ///Vec2 pa = positionBuffer.data[a];
          const pa = s_pa.Copy(pos_data[a]);
          ///Vec2 pb = positionBuffer.data[b];
          const pb = s_pb.Copy(pos_data[b]);
          ///Vec2& va = velocityBuffer.data[a];
          const va = vel_data[a];
          ///Vec2& vb = velocityBuffer.data[b];
          const vb = vel_data[b];
          ///pa += step.dt * va;
          pa.SelfMulAdd(step.dt, va);
          ///pb += step.dt * vb;
          pb.SelfMulAdd(step.dt, vb);
          ///Vec2 d = pb - pa;
          const d = Vec2.SubVV(pb, pa, s_d);
          ///float32 r0 = pair.distance;
          const r0 = pair.distance;
          ///float32 r1 = d.Length();
          const r1 = d.Length();
          ///float32 strength = springStrength * pair.strength;
          const strength = springStrength * pair.strength;
          ///Vec2 f = strength * (r0 - r1) / r1 * d;
          const f = Vec2.MulSV(strength * (r0 - r1) / r1, d, s_f);
          ///va -= f;
          va.SelfSub(f);
          ///vb += f;
          vb.SelfAdd(f);
        }
      }
    }
    public static readonly SolveSpring_s_pa = new Vec2();
    public static readonly SolveSpring_s_pb = new Vec2();
    public static readonly SolveSpring_s_d = new Vec2();
    public static readonly SolveSpring_s_f = new Vec2();

    public SolveTensile(step: TimeStep): void {
      const s_weightedNormal = ParticleSystem.SolveTensile_s_weightedNormal;
      const s_s = ParticleSystem.SolveTensile_s_s;
      const s_f = ParticleSystem.SolveTensile_s_f;
      const vel_data = this.velocityBuffer.data;
      // DEBUG: Assert(this.accumulation2Buffer !== null);
      for (let i = 0; i < this.count; i++) {
        this.accumulation2Buffer[i] = new Vec2();
        this.accumulation2Buffer[i].SetZero();
      }
      for (let k = 0; k < this.contactBuffer.count; k++) {
        const contact = this.contactBuffer.data[k];
        if (contact.flags & ParticleFlag.tensileParticle) {
          const a = contact.indexA;
          const b = contact.indexB;
          const w = contact.weight;
          const n = contact.normal;
          ///Vec2 weightedNormal = (1 - w) * w * n;
          const weightedNormal = Vec2.MulSV((1 - w) * w, n, s_weightedNormal);
          ///accumulation2Buffer[a] -= weightedNormal;
          this.accumulation2Buffer[a].SelfSub(weightedNormal);
          ///accumulation2Buffer[b] += weightedNormal;
          this.accumulation2Buffer[b].SelfAdd(weightedNormal);
        }
      }
      const criticalVelocity = this.GetCriticalVelocity(step);
      const pressureStrength = this.def.surfaceTensionPressureStrength * criticalVelocity;
      const normalStrength = this.def.surfaceTensionNormalStrength * criticalVelocity;
      const maxVelocityVariation = maxParticleForce * criticalVelocity;
      for (let k = 0; k < this.contactBuffer.count; k++) {
        const contact = this.contactBuffer.data[k];
        if (contact.flags & ParticleFlag.tensileParticle) {
          const a = contact.indexA;
          const b = contact.indexB;
          const w = contact.weight;
          const n = contact.normal;
          const h = this.weightBuffer[a] + this.weightBuffer[b];
          ///Vec2 s = accumulation2Buffer[b] - accumulation2Buffer[a];
          const s = Vec2.SubVV(this.accumulation2Buffer[b], this.accumulation2Buffer[a], s_s);
          const fn = Min(
            pressureStrength * (h - 2) + normalStrength * Vec2.DotVV(s, n),
            maxVelocityVariation) * w;
          ///Vec2 f = fn * n;
          const f = Vec2.MulSV(fn, n, s_f);
          ///velocityBuffer.data[a] -= f;
          vel_data[a].SelfSub(f);
          ///velocityBuffer.data[b] += f;
          vel_data[b].SelfAdd(f);
        }
      }
    }
    public static readonly SolveTensile_s_weightedNormal = new Vec2();
    public static readonly SolveTensile_s_s = new Vec2();
    public static readonly SolveTensile_s_f = new Vec2();

    public SolveViscous(): void {
      const s_v = ParticleSystem.SolveViscous_s_v;
      const s_f = ParticleSystem.SolveViscous_s_f;
      const pos_data = this.positionBuffer.data;
      const vel_data = this.velocityBuffer.data;
      const viscousStrength = this.def.viscousStrength;
      const inv_mass = this.GetParticleInvMass();
      for (let k = 0; k < this.bodyContactBuffer.count; k++) {
        const contact = this.bodyContactBuffer.data[k];
        const a = contact.index;
        if (this.flagsBuffer.data[a] & ParticleFlag.viscousParticle) {
          const b = contact.body;
          const w = contact.weight;
          const m = contact.mass;
          const p = pos_data[a];
          ///Vec2 v = b.GetLinearVelocityFromWorldPoint(p) - velocityBuffer.data[a];
          const v = Vec2.SubVV(b.GetLinearVelocityFromWorldPoint(p, Vec2.s_t0), vel_data[a], s_v);
          ///Vec2 f = viscousStrength * m * w * v;
          const f = Vec2.MulSV(viscousStrength * m * w, v, s_f);
          ///velocityBuffer.data[a] += GetParticleInvMass() * f;
          vel_data[a].SelfMulAdd(inv_mass, f);
          ///b.ApplyLinearImpulse(-f, p, true);
          b.ApplyLinearImpulse(f.SelfNeg(), p, true);
        }
      }
      for (let k = 0; k < this.contactBuffer.count; k++) {
        const contact = this.contactBuffer.data[k];
        if (contact.flags & ParticleFlag.viscousParticle) {
          const a = contact.indexA;
          const b = contact.indexB;
          const w = contact.weight;
          ///Vec2 v = velocityBuffer.data[b] - velocityBuffer.data[a];
          const v = Vec2.SubVV(vel_data[b], vel_data[a], s_v);
          ///Vec2 f = viscousStrength * w * v;
          const f = Vec2.MulSV(viscousStrength * w, v, s_f);
          ///velocityBuffer.data[a] += f;
          vel_data[a].SelfAdd(f);
          ///velocityBuffer.data[b] -= f;
          vel_data[b].SelfSub(f);
        }
      }
    }
    public static readonly SolveViscous_s_v = new Vec2();
    public static readonly SolveViscous_s_f = new Vec2();

    public SolveRepulsive(step: TimeStep): void {
      const s_f = ParticleSystem.SolveRepulsive_s_f;
      const vel_data = this.velocityBuffer.data;
      const repulsiveStrength = this.def.repulsiveStrength * this.GetCriticalVelocity(step);
      for (let k = 0; k < this.contactBuffer.count; k++) {
        const contact = this.contactBuffer.data[k];
        if (contact.flags & ParticleFlag.repulsiveParticle) {
          const a = contact.indexA;
          const b = contact.indexB;
          if (this.groupBuffer[a] !== this.groupBuffer[b]) {
            const w = contact.weight;
            const n = contact.normal;
            ///Vec2 f = repulsiveStrength * w * n;
            const f = Vec2.MulSV(repulsiveStrength * w, n, s_f);
            ///velocityBuffer.data[a] -= f;
            vel_data[a].SelfSub(f);
            ///velocityBuffer.data[b] += f;
            vel_data[b].SelfAdd(f);
          }
        }
      }
    }
    public static readonly SolveRepulsive_s_f = new Vec2();

    public SolvePowder(step: TimeStep): void {
      const s_f = ParticleSystem.SolvePowder_s_f;
      const pos_data = this.positionBuffer.data;
      const vel_data = this.velocityBuffer.data;
      const powderStrength = this.def.powderStrength * this.GetCriticalVelocity(step);
      const minWeight = 1.0 - particleStride;
      const inv_mass = this.GetParticleInvMass();
      for (let k = 0; k < this.bodyContactBuffer.count; k++) {
        const contact = this.bodyContactBuffer.data[k];
        const a = contact.index;
        if (this.flagsBuffer.data[a] & ParticleFlag.powderParticle) {
          const w = contact.weight;
          if (w > minWeight) {
            const b = contact.body;
            const m = contact.mass;
            const p = pos_data[a];
            const n = contact.normal;
            const f = Vec2.MulSV(powderStrength * m * (w - minWeight), n, s_f);
            vel_data[a].SelfMulSub(inv_mass, f);
            b.ApplyLinearImpulse(f, p, true);
          }
        }
      }
      for (let k = 0; k < this.contactBuffer.count; k++) {
        const contact = this.contactBuffer.data[k];
        if (contact.flags & ParticleFlag.powderParticle) {
          const w = contact.weight;
          if (w > minWeight) {
            const a = contact.indexA;
            const b = contact.indexB;
            const n = contact.normal;
            const f = Vec2.MulSV(powderStrength * (w - minWeight), n, s_f);
            vel_data[a].SelfSub(f);
            vel_data[b].SelfAdd(f);
          }
        }
      }
    }
    public static readonly SolvePowder_s_f = new Vec2();

    public SolveSolid(step: TimeStep): void {
      const s_f = ParticleSystem.SolveSolid_s_f;
      const vel_data = this.velocityBuffer.data;
      // applies extra repulsive force from solid particle groups
      this.depthBuffer = this.RequestBuffer(this.depthBuffer);
      const ejectionStrength = step.inv_dt * this.def.ejectionStrength;
      for (let k = 0; k < this.contactBuffer.count; k++) {
        const contact = this.contactBuffer.data[k];
        const a = contact.indexA;
        const b = contact.indexB;
        if (this.groupBuffer[a] !== this.groupBuffer[b]) {
          const w = contact.weight;
          const n = contact.normal;
          const h = this.depthBuffer[a] + this.depthBuffer[b];
          const f = Vec2.MulSV(ejectionStrength * h * w, n, s_f);
          vel_data[a].SelfSub(f);
          vel_data[b].SelfAdd(f);
        }
      }
    }
    public static readonly SolveSolid_s_f = new Vec2();

    public SolveForce(step: TimeStep): void {
      const vel_data = this.velocityBuffer.data;
      const velocityPerForce = step.dt * this.GetParticleInvMass();
      for (let i = 0; i < this.count; i++) {
        ///velocityBuffer.data[i] += velocityPerForce * forceBuffer[i];
        vel_data[i].SelfMulAdd(velocityPerForce, this.forceBuffer[i]);
      }
      this.hasForce = false;
    }

    public SolveColorMixing(): void {
      // mixes color between contacting particles
      const colorMixing = 0.5 * this.def.colorMixingStrength;
      if (colorMixing) {
        for (let k = 0; k < this.contactBuffer.count; k++) {
          const contact = this.contactBuffer.data[k];
          const a = contact.indexA;
          const b = contact.indexB;
          if (this.flagsBuffer.data[a] & this.flagsBuffer.data[b] &
            ParticleFlag.colorMixingParticle) {
            const colorA = this.colorBuffer.data[a];
            const colorB = this.colorBuffer.data[b];
            // Use the static method to ensure certain compilers inline
            // this correctly.
            Color.MixColors(colorA, colorB, colorMixing);
          }
        }
      }
    }

    public SolveZombie(): void {
      // removes particles with zombie flag
      let newCount = 0;
      const newIndices: number[] = []; // TODO: static
      for (let i = 0; i < this.count; i++) {
        newIndices[i] = invalidParticleIndex;
      }
      // DEBUG: Assert(newIndices.length === this.count);
      let allParticleFlags = 0;
      for (let i = 0; i < this.count; i++) {
        const flags = this.flagsBuffer.data[i];
        if (flags & ParticleFlag.zombieParticle) {
          const destructionListener = this.world.destructionListener;
          if ((flags & ParticleFlag.destructionListenerParticle) && destructionListener) {
            destructionListener.SayGoodbyeParticle(this, i);
          }
          // Destroy particle handle.
          if (this.handleIndexBuffer.data) {
            const handle = this.handleIndexBuffer.data[i];
            if (handle) {
              handle.SetIndex(invalidParticleIndex);
              this.handleIndexBuffer.data[i] = null;
              ///handleAllocator.Free(handle);
            }
          }
          newIndices[i] = invalidParticleIndex;
        } else {
          newIndices[i] = newCount;
          if (i !== newCount) {
            // Update handle to reference new particle index.
            if (this.handleIndexBuffer.data) {
              const handle = this.handleIndexBuffer.data[i];
              if (handle) { handle.SetIndex(newCount); }
              this.handleIndexBuffer.data[newCount] = handle;
            }
            this.flagsBuffer.data[newCount] = this.flagsBuffer.data[i];
            if (this.lastBodyContactStepBuffer.data) {
              this.lastBodyContactStepBuffer.data[newCount] = this.lastBodyContactStepBuffer.data[i];
            }
            if (this.bodyContactCountBuffer.data) {
              this.bodyContactCountBuffer.data[newCount] = this.bodyContactCountBuffer.data[i];
            }
            if (this.consecutiveContactStepsBuffer.data) {
              this.consecutiveContactStepsBuffer.data[newCount] = this.consecutiveContactStepsBuffer.data[i];
            }
            this.positionBuffer.data[newCount].Copy(this.positionBuffer.data[i]);
            this.velocityBuffer.data[newCount].Copy(this.velocityBuffer.data[i]);
            this.groupBuffer[newCount] = this.groupBuffer[i];
            if (this.hasForce) {
              this.forceBuffer[newCount].Copy(this.forceBuffer[i]);
            }
            if (this.staticPressureBuffer) {
              this.staticPressureBuffer[newCount] = this.staticPressureBuffer[i];
            }
            if (this.depthBuffer) {
              this.depthBuffer[newCount] = this.depthBuffer[i];
            }
            if (this.colorBuffer.data) {
              this.colorBuffer.data[newCount].Copy(this.colorBuffer.data[i]);
            }
            if (this.userDataBuffer.data) {
              this.userDataBuffer.data[newCount] = this.userDataBuffer.data[i];
            }
            if (this.expirationTimeBuffer.data) {
              this.expirationTimeBuffer.data[newCount] = this.expirationTimeBuffer.data[i];
            }
          }
          newCount++;
          allParticleFlags |= flags;
        }
      }

      // predicate functions
      const Test = {
        ///static bool IsProxyInvalid(const Proxy& proxy)
        IsProxyInvalid: (proxy: ParticleSysteProxy) => {
          return proxy.index < 0;
        },
        ///static bool IsContactInvalid(const ParticleContact& contact)
        IsContactInvalid: (contact: ParticleContact) => {
          return contact.indexA < 0 || contact.indexB < 0;
        },
        ///static bool IsBodyContactInvalid(const ParticleBodyContact& contact)
        IsBodyContactInvalid: (contact: ParticleBodyContact) => {
          return contact.index < 0;
        },
        ///static bool IsPairInvalid(const ParticlePair& pair)
        IsPairInvalid: (pair: ParticlePair) => {
          return pair.indexA < 0 || pair.indexB < 0;
        },
        ///static bool IsTriadInvalid(const ParticleTriad& triad)
        IsTriadInvalid: (triad: ParticleTriad) => {
          return triad.indexA < 0 || triad.indexB < 0 || triad.indexC < 0;
        },
      };

      // update proxies
      for (let k = 0; k < this.proxyBuffer.count; k++) {
        const proxy = this.proxyBuffer.data[k];
        proxy.index = newIndices[proxy.index];
      }
      this.proxyBuffer.RemoveIf(Test.IsProxyInvalid);

      // update contacts
      for (let k = 0; k < this.contactBuffer.count; k++) {
        const contact = this.contactBuffer.data[k];
        contact.indexA = newIndices[contact.indexA];
        contact.indexB = newIndices[contact.indexB];
      }
      this.contactBuffer.RemoveIf(Test.IsContactInvalid);

      // update particle-body contacts
      for (let k = 0; k < this.bodyContactBuffer.count; k++) {
        const contact = this.bodyContactBuffer.data[k];
        contact.index = newIndices[contact.index];
      }
      this.bodyContactBuffer.RemoveIf(Test.IsBodyContactInvalid);

      // update pairs
      for (let k = 0; k < this.pairBuffer.count; k++) {
        const pair = this.pairBuffer.data[k];
        pair.indexA = newIndices[pair.indexA];
        pair.indexB = newIndices[pair.indexB];
      }
      this.pairBuffer.RemoveIf(Test.IsPairInvalid);

      // update triads
      for (let k = 0; k < this.triadBuffer.count; k++) {
        const triad = this.triadBuffer.data[k];
        triad.indexA = newIndices[triad.indexA];
        triad.indexB = newIndices[triad.indexB];
        triad.indexC = newIndices[triad.indexC];
      }
      this.triadBuffer.RemoveIf(Test.IsTriadInvalid);

      // Update lifetime indices.
      if (this.indexByExpirationTimeBuffer.data) {
        let writeOffset = 0;
        for (let readOffset = 0; readOffset < this.count; readOffset++) {
          const newIndex = newIndices[this.indexByExpirationTimeBuffer.data[readOffset]];
          if (newIndex !== invalidParticleIndex) {
            this.indexByExpirationTimeBuffer.data[writeOffset++] = newIndex;
          }
        }
      }

      // update groups
      for (let group = this.groupList; group; group = group.GetNext()) {
        let firstIndex = newCount;
        let lastIndex = 0;
        let modified = false;
        for (let i = group.firstIndex; i < group.lastIndex; i++) {
          const j = newIndices[i];
          if (j >= 0) {
            firstIndex = Min(firstIndex, j);
            lastIndex = Max(lastIndex, j + 1);
          } else {
            modified = true;
          }
        }
        if (firstIndex < lastIndex) {
          group.firstIndex = firstIndex;
          group.lastIndex = lastIndex;
          if (modified) {
            if (group.groupFlags & ParticleGroupFlag.solidParticleGroup) {
              this.SetGroupFlags(group, group.groupFlags | ParticleGroupFlag.particleGroupNeedsUpdateDepth);
            }
          }
        } else {
          group.firstIndex = 0;
          group.lastIndex = 0;
          if (!(group.groupFlags & ParticleGroupFlag.particleGroupCanBeEmpty)) {
            this.SetGroupFlags(group, group.groupFlags | ParticleGroupFlag.particleGroupWillBeDestroyed);
          }
        }
      }

      // update particle count
      this.count = newCount;
      this.allParticleFlags = allParticleFlags;
      this.needsUpdateAllParticleFlags = false;

      // destroy bodies with no particles
      for (let group = this.groupList; group; ) {
        const next = group.GetNext();
        if (group.groupFlags & ParticleGroupFlag.particleGroupWillBeDestroyed) {
          this.DestroyParticleGroup(group);
        }
        group = next;
      }
    }

    /**
     * Destroy all particles which have outlived their lifetimes set
     * by SetParticleLifetime().
     */
    public SolveLifetimes(step: TimeStep): void {
      // Update the time elapsed.
      this.timeElapsed = this.LifetimeToExpirationTime(step.dt);
      // Get the floor (non-fractional component) of the elapsed time.
      const quantizedTimeElapsed = this.GetQuantizedTimeElapsed();

      const expirationTimes = this.expirationTimeBuffer.data;
      const expirationTimeIndices = this.indexByExpirationTimeBuffer.data;
      const particleCount = this.GetParticleCount();
      // Sort the lifetime buffer if it's required.
      if (this.expirationTimeBufferRequiresSorting) {
        ///const ExpirationTimeComparator expirationTimeComparator(expirationTimes);
        ///std::sort(expirationTimeIndices, expirationTimeIndices + particleCount, expirationTimeComparator);

        /**
         * Compare the lifetime of particleIndexA and particleIndexB
         * returning true if the lifetime of A is greater than B for
         * particles that will expire.  If either particle's lifetime is
         * infinite (<= 0.0f) this function return true if the lifetime
         * of A is lesser than B. When used with std::sort() this
         * results in an array of particle indicies sorted in reverse
         * order by particle lifetime.
         *
         * For example, the set of lifetimes
         * (1.0, 0.7, 0.3, 0.0, -1.0, 2.0)
         * would be sorted as
         * (0.0, 1.0, -2.0, 1.0, 0.7, 0.3)
         */
        const ExpirationTimeComparator = (particleIndexA: number, particleIndexB: number): boolean => {
          const expirationTimeA = expirationTimes[particleIndexA];
          const expirationTimeB = expirationTimes[particleIndexB];
          const infiniteExpirationTimeA = expirationTimeA <= 0.0;
          const infiniteExpirationTimeB = expirationTimeB <= 0.0;
          return infiniteExpirationTimeA === infiniteExpirationTimeB ?
            expirationTimeA > expirationTimeB : infiniteExpirationTimeA;
        };

        std_sort(expirationTimeIndices, 0, particleCount, ExpirationTimeComparator);

        this.expirationTimeBufferRequiresSorting = false;
      }

      // Destroy particles which have expired.
      for (let i = particleCount - 1; i >= 0; --i) {
        const particleIndex = expirationTimeIndices[i];
        const expirationTime = expirationTimes[particleIndex];
        // If no particles need to be destroyed, skip this.
        if (quantizedTimeElapsed < expirationTime || expirationTime <= 0) {
          break;
        }
        // Destroy this particle.
        this.DestroyParticle(particleIndex);
      }
    }

    public RotateBuffer(start: number, mid: number, end: number): void {
      // move the particles assigned to the given group toward the end of array
      if (start === mid || mid === end) {
        return;
      }
      // DEBUG: Assert(mid >= start && mid <= end);

      function newIndices(i: number): number {
        if (i < start) {
          return i;
        } else if (i < mid) {
          return i + end - mid;
        } else if (i < end) {
          return i + start - mid;
        } else {
          return i;
        }
      }

      ///std::rotate(flagsBuffer.data + start, flagsBuffer.data + mid, flagsBuffer.data + end);
      std_rotate(this.flagsBuffer.data, start, mid, end);
      if (this.lastBodyContactStepBuffer.data) {
        ///std::rotate(lastBodyContactStepBuffer.data + start, lastBodyContactStepBuffer.data + mid, lastBodyContactStepBuffer.data + end);
        std_rotate(this.lastBodyContactStepBuffer.data, start, mid, end);
      }
      if (this.bodyContactCountBuffer.data) {
        ///std::rotate(bodyContactCountBuffer.data + start, bodyContactCountBuffer.data + mid, bodyContactCountBuffer.data + end);
        std_rotate(this.bodyContactCountBuffer.data, start, mid, end);
      }
      if (this.consecutiveContactStepsBuffer.data) {
        ///std::rotate(consecutiveContactStepsBuffer.data + start, consecutiveContactStepsBuffer.data + mid, consecutiveContactStepsBuffer.data + end);
        std_rotate(this.consecutiveContactStepsBuffer.data, start, mid, end);
      }
      ///std::rotate(positionBuffer.data + start, positionBuffer.data + mid, positionBuffer.data + end);
      std_rotate(this.positionBuffer.data, start, mid, end);
      ///std::rotate(velocityBuffer.data + start, velocityBuffer.data + mid, velocityBuffer.data + end);
      std_rotate(this.velocityBuffer.data, start, mid, end);
      ///std::rotate(groupBuffer + start, groupBuffer + mid, groupBuffer + end);
      std_rotate(this.groupBuffer, start, mid, end);
      if (this.hasForce) {
        ///std::rotate(forceBuffer + start, forceBuffer + mid, forceBuffer + end);
        std_rotate(this.forceBuffer, start, mid, end);
      }
      if (this.staticPressureBuffer) {
        ///std::rotate(staticPressureBuffer + start, staticPressureBuffer + mid, staticPressureBuffer + end);
        std_rotate(this.staticPressureBuffer, start, mid, end);
      }
      if (this.depthBuffer) {
        ///std::rotate(depthBuffer + start, depthBuffer + mid, depthBuffer + end);
        std_rotate(this.depthBuffer, start, mid, end);
      }
      if (this.colorBuffer.data) {
        ///std::rotate(colorBuffer.data + start, colorBuffer.data + mid, colorBuffer.data + end);
        std_rotate(this.colorBuffer.data, start, mid, end);
      }
      if (this.userDataBuffer.data) {
        ///std::rotate(userDataBuffer.data + start, userDataBuffer.data + mid, userDataBuffer.data + end);
        std_rotate(this.userDataBuffer.data, start, mid, end);
      }

      // Update handle indices.
      if (this.handleIndexBuffer.data) {
        ///std::rotate(handleIndexBuffer.data + start, handleIndexBuffer.data + mid, handleIndexBuffer.data + end);
        std_rotate(this.handleIndexBuffer.data, start, mid, end);
        for (let i = start; i < end; ++i) {
          const handle = this.handleIndexBuffer.data[i];
          if (handle) { handle.SetIndex(newIndices(handle.GetIndex())); }
        }
      }

      if (this.expirationTimeBuffer.data) {
        ///std::rotate(expirationTimeBuffer.data + start, expirationTimeBuffer.data + mid, expirationTimeBuffer.data + end);
        std_rotate(this.expirationTimeBuffer.data, start, mid, end);
        // Update expiration time buffer indices.
        const particleCount = this.GetParticleCount();
        const indexByExpirationTime = this.indexByExpirationTimeBuffer.data;
        for (let i = 0; i < particleCount; ++i) {
          indexByExpirationTime[i] = newIndices(indexByExpirationTime[i]);
        }
      }

      // update proxies
      for (let k = 0; k < this.proxyBuffer.count; k++) {
        const proxy = this.proxyBuffer.data[k];
        proxy.index = newIndices(proxy.index);
      }

      // update contacts
      for (let k = 0; k < this.contactBuffer.count; k++) {
        const contact = this.contactBuffer.data[k];
        contact.indexA = newIndices(contact.indexA);
        contact.indexB = newIndices(contact.indexB);
      }

      // update particle-body contacts
      for (let k = 0; k < this.bodyContactBuffer.count; k++) {
        const contact = this.bodyContactBuffer.data[k];
        contact.index = newIndices(contact.index);
      }

      // update pairs
      for (let k = 0; k < this.pairBuffer.count; k++) {
        const pair = this.pairBuffer.data[k];
        pair.indexA = newIndices(pair.indexA);
        pair.indexB = newIndices(pair.indexB);
      }

      // update triads
      for (let k = 0; k < this.triadBuffer.count; k++) {
        const triad = this.triadBuffer.data[k];
        triad.indexA = newIndices(triad.indexA);
        triad.indexB = newIndices(triad.indexB);
        triad.indexC = newIndices(triad.indexC);
      }

      // update groups
      for (let group = this.groupList; group; group = group.GetNext()) {
        group.firstIndex = newIndices(group.firstIndex);
        group.lastIndex = newIndices(group.lastIndex - 1) + 1;
      }
    }

    public GetCriticalVelocity(step: TimeStep): number {
      return this.particleDiameter * step.inv_dt;
    }

    public GetCriticalVelocitySquared(step: TimeStep): number {
      const velocity = this.GetCriticalVelocity(step);
      return velocity * velocity;
    }

    public GetCriticalPressure(step: TimeStep): number {
      return this.def.density * this.GetCriticalVelocitySquared(step);
    }

    public GetParticleStride(): number {
      return particleStride * this.particleDiameter;
    }

    public GetParticleMass(): number {
      const stride = this.GetParticleStride();
      return this.def.density * stride * stride;
    }

    public GetParticleInvMass(): number {
      ///return 1.777777 * this.inverseDensity * this.inverseDiameter * this.inverseDiameter;
      // mass = density * stride^2, so we take the inverse of this.
      const inverseStride = this.inverseDiameter * (1.0 / particleStride);
      return this.inverseDensity * inverseStride * inverseStride;
    }

    /**
     * Get the world's contact filter if any particles with the
     * contactFilterParticle flag are present in the system.
     */
    public GetFixtureContactFilter(): ContactFilter | null {
      return (this.allParticleFlags & ParticleFlag.fixtureContactFilterParticle) ?
        this.world.contactManager.contactFilter : null;
    }

    /**
     * Get the world's contact filter if any particles with the
     * particleContactFilterParticle flag are present in the
     * system.
     */
    public GetParticleContactFilter(): ContactFilter | null {
      return (this.allParticleFlags & ParticleFlag.particleContactFilterParticle) ?
        this.world.contactManager.contactFilter : null;
    }

    /**
     * Get the world's contact listener if any particles with the
     * fixtureContactListenerParticle flag are present in the
     * system.
     */
    public GetFixtureContactListener(): ContactListener | null {
      return (this.allParticleFlags & ParticleFlag.fixtureContactListenerParticle) ?
        this.world.contactManager.contactListener : null;
    }

    /**
     * Get the world's contact listener if any particles with the
     * particleContactListenerParticle flag are present in the
     * system.
     */
    public GetParticleContactListener(): ContactListener | null {
      return (this.allParticleFlags & ParticleFlag.particleContactListenerParticle) ?
        this.world.contactManager.contactListener : null;
    }

    public SetUserOverridableBuffer<T>(buffer: ParticleSysteUserOverridableBuffer<T>, data: T[]): void {
      buffer.data = data;
      buffer.userSuppliedCapacity = data.length;
    }

    public SetGroupFlags(group: ParticleGroup, newFlags: ParticleGroupFlag): void {
      const oldFlags = group.groupFlags;
      if ((oldFlags ^ newFlags) & ParticleGroupFlag.solidParticleGroup) {
        // If the solidParticleGroup flag changed schedule depth update.
        newFlags |= ParticleGroupFlag.particleGroupNeedsUpdateDepth;
      }
      if (oldFlags & ~newFlags) {
        // If any flags might be removed
        this.needsUpdateAllGroupFlags = true;
      }
      if (~this.allGroupFlags & newFlags) {
        // If any flags were added
        if (newFlags & ParticleGroupFlag.solidParticleGroup) {
          this.depthBuffer = this.RequestBuffer(this.depthBuffer);
        }
        this.allGroupFlags |= newFlags;
      }
      group.groupFlags = newFlags;
    }

    public static BodyContactCompare(lhs: ParticleBodyContact, rhs: ParticleBodyContact): boolean {
      if (lhs.index === rhs.index) {
        // Subsort by weight, decreasing.
        return lhs.weight > rhs.weight;
      }
      return lhs.index < rhs.index;
    }

    public RemoveSpuriousBodyContacts(): void {
      // At this point we have a list of contact candidates based on AABB
      // overlap.The AABB query that  generated this returns all collidable
      // fixtures overlapping particle bounding boxes.  This breaks down around
      // vertices where two shapes intersect, such as a "ground" surface made
      // of multiple PolygonShapes; it potentially applies a lot of spurious
      // impulses from normals that should not actually contribute.  See the
      // Ramp example in Testbed.
      //
      // To correct for this, we apply this algorithm:
      //   * sort contacts by particle and subsort by weight (nearest to farthest)
      //   * for each contact per particle:
      //      - project a point at the contact distance along the inverse of the
      //        contact normal
      //      - if this intersects the fixture that generated the contact, apply
      //         it, otherwise discard as impossible
      //      - repeat for up to n nearest contacts, currently we get good results
      //        from n=3.
      ///std::sort(bodyContactBuffer.Begin(), bodyContactBuffer.End(), ParticleSystem::BodyContactCompare);
      std_sort(this.bodyContactBuffer.data, 0, this.bodyContactBuffer.count, ParticleSystem.BodyContactCompare);

      ///int32 discarded = 0;
      ///std::remove_if(bodyContactBuffer.Begin(), bodyContactBuffer.End(), ParticleBodyContactRemovePredicate(this, &discarded));
      ///
      ///bodyContactBuffer.SetCount(bodyContactBuffer.GetCount() - discarded);

      const s_n = ParticleSystem.RemoveSpuriousBodyContacts_s_n;
      const s_pos = ParticleSystem.RemoveSpuriousBodyContacts_s_pos;
      const s_normal = ParticleSystem.RemoveSpuriousBodyContacts_s_normal;

      // Max number of contacts processed per particle, from nearest to farthest.
      // This must be at least 2 for correctness with concave shapes; 3 was
      // experimentally arrived at as looking reasonable.
      const k_maxContactsPerPoint = 3;
      const system = this;
      // Index of last particle processed.
      let lastIndex = -1;
      // Number of contacts processed for the current particle.
      let currentContacts = 0;
      // Output the number of discarded contacts.
      // let discarded = 0;
      const ParticleBodyContactRemovePredicate = (contact: ParticleBodyContact): boolean => {
        // This implements the selection criteria described in
        // RemoveSpuriousBodyContacts().
        // This functor is iterating through a list of Body contacts per
        // Particle, ordered from near to far.  For up to the maximum number of
        // contacts we allow per point per step, we verify that the contact
        // normal of the Body that genenerated the contact makes physical sense
        // by projecting a point back along that normal and seeing if it
        // intersects the fixture generating the contact.

        if (contact.index !== lastIndex) {
          currentContacts = 0;
          lastIndex = contact.index;
        }

        if (currentContacts++ > k_maxContactsPerPoint) {
          // ++discarded;
          return true;
        }

        // Project along inverse normal (as returned in the contact) to get the
        // point to check.
        ///Vec2 n = contact.normal;
        const n = s_n.Copy(contact.normal);
        // weight is 1-(inv(diameter) * distance)
        ///n *= system.particleDiameter * (1 - contact.weight);
        n.SelfMul(system.particleDiameter * (1 - contact.weight));
        ///Vec2 pos = system.positionBuffer.data[contact.index] + n;
        const pos = Vec2.AddVV(system.positionBuffer.data[contact.index], n, s_pos);

        // pos is now a point projected back along the contact normal to the
        // contact distance. If the surface makes sense for a contact, pos will
        // now lie on or in the fixture generating
        if (!contact.fixture.TestPoint(pos)) {
          const childCount = contact.fixture.GetShape().GetChildCount();
          for (let childIndex = 0; childIndex < childCount; childIndex++) {
            const normal = s_normal;
            const distance = contact.fixture.ComputeDistance(pos, normal, childIndex);
            if (distance < linearSlop) {
              return false;
            }
          }
          // ++discarded;
          return true;
        }

        return false;
      };
      this.bodyContactBuffer.count = std_remove_if(this.bodyContactBuffer.data, ParticleBodyContactRemovePredicate, this.bodyContactBuffer.count);
    }
    private static RemoveSpuriousBodyContacts_s_n = new Vec2();
    private static RemoveSpuriousBodyContacts_s_pos = new Vec2();
    private static RemoveSpuriousBodyContacts_s_normal = new Vec2();

    public DetectStuckParticle(particle: number): void {
      // Detect stuck particles
      //
      // The basic algorithm is to allow the user to specify an optional
      // threshold where we detect whenever a particle is contacting
      // more than one fixture for more than threshold consecutive
      // steps. This is considered to be "stuck", and these are put
      // in a list the user can query per step, if enabled, to deal with
      // such particles.

      if (this.stuckThreshold <= 0) {
        return;
      }

      // Get the state variables for this particle.
      ///int32 * const consecutiveCount = &consecutiveContactStepsBuffer.data[particle];
      ///int32 * const lastStep = &lastBodyContactStepBuffer.data[particle];
      ///int32 * const bodyCount = &bodyContactCountBuffer.data[particle];

      // This is only called when there is a body contact for this particle.
      ///++(*bodyCount);
      ++this.bodyContactCountBuffer.data[particle];

      // We want to only trigger detection once per step, the first time we
      // contact more than one fixture in a step for a given particle.
      ///if (*bodyCount === 2)
      if (this.bodyContactCountBuffer.data[particle] === 2) {
        ///++(*consecutiveCount);
        ++this.consecutiveContactStepsBuffer.data[particle];
        ///if (*consecutiveCount > stuckThreshold)
        if (this.consecutiveContactStepsBuffer.data[particle] > this.stuckThreshold) {
          ///int32& newStuckParticle = stuckParticleBuffer.Append();
          ///newStuckParticle = particle;
          this.stuckParticleBuffer.data[this.stuckParticleBuffer.Append()] = particle;
        }
      }
      ///*lastStep = timestamp;
      this.lastBodyContactStepBuffer.data[particle] = this.timestamp;
    }

    /**
     * Determine whether a particle index is valid.
     */
    public ValidateParticleIndex(index: number): boolean {
      return index >= 0 && index < this.GetParticleCount() &&
        index !== invalidParticleIndex;
    }

    /**
     * Get the time elapsed in
     * ParticleSystemDef::lifetimeGranularity.
     */
    public GetQuantizedTimeElapsed(): number {
      ///return (int32)(timeElapsed >> 32);
      return Math.floor(this.timeElapsed / 0x100000000);
    }

    /**
     * Convert a lifetime in seconds to an expiration time.
     */
    public LifetimeToExpirationTime(lifetime: number): number {
      ///return timeElapsed + (int64)((lifetime / def.lifetimeGranularity) * (float32)(1LL << 32));
      return this.timeElapsed + Math.floor(((lifetime / this.def.lifetimeGranularity) * 0x100000000));
    }

    public ForceCanBeApplied(flags: ParticleFlag): boolean {
      return !(flags & ParticleFlag.wallParticle);
    }

    public PrepareForceBuffer(): void {
      if (!this.hasForce) {
        ///memset(forceBuffer, 0, sizeof(*forceBuffer) * count);
        for (let i = 0; i < this.count; i++) {
          this.forceBuffer[i].SetZero();
        }
        this.hasForce = true;
      }
    }

    public IsRigidGroup(group: ParticleGroup | null): boolean {
      return (group !== null) && ((group.groupFlags & ParticleGroupFlag.rigidParticleGroup) !== 0);
    }

    public GetLinearVelocity(group: ParticleGroup | null, particleIndex: number, point: Vec2, out: Vec2): Vec2 {
      if (group && this.IsRigidGroup(group)) {
        return group.GetLinearVelocityFromWorldPoint(point, out);
      } else {
        ///return velocityBuffer.data[particleIndex];
        return out.Copy(this.velocityBuffer.data[particleIndex]);
      }
    }

    public InitDampingParameter(invMass: number[], invInertia: number[], tangentDistance: number[], mass: number, inertia: number, center: Vec2, point: Vec2, normal: Vec2): void {
      ///*invMass = mass > 0 ? 1 / mass : 0;
      invMass[0] = mass > 0 ? 1 / mass : 0;
      ///*invInertia = inertia > 0 ? 1 / inertia : 0;
      invInertia[0] = inertia > 0 ? 1 / inertia : 0;
      ///*tangentDistance = Cross(point - center, normal);
      tangentDistance[0] = Vec2.CrossVV(Vec2.SubVV(point, center, Vec2.s_t0), normal);
    }

    public InitDampingParameterWithRigidGroupOrParticle(invMass: number[], invInertia: number[], tangentDistance: number[], isRigidGroup: boolean, group: ParticleGroup | null, particleIndex: number, point: Vec2, normal: Vec2): void {
      if (group && isRigidGroup) {
        this.InitDampingParameter(invMass, invInertia, tangentDistance, group.GetMass(), group.GetInertia(), group.GetCenter(), point, normal);
      } else {
        const flags = this.flagsBuffer.data[particleIndex];
        this.InitDampingParameter(invMass, invInertia, tangentDistance, flags & ParticleFlag.wallParticle ? 0 : this.GetParticleMass(), 0, point, point, normal);
      }
    }

    public ComputeDampingImpulse(invMassA: number, invInertiaA: number, tangentDistanceA: number, invMassB: number, invInertiaB: number, tangentDistanceB: number, normalVelocity: number): number {
      const invMass =
        invMassA + invInertiaA * tangentDistanceA * tangentDistanceA +
        invMassB + invInertiaB * tangentDistanceB * tangentDistanceB;
      return invMass > 0 ? normalVelocity / invMass : 0;
    }

    public ApplyDamping(invMass: number, invInertia: number, tangentDistance: number, isRigidGroup: boolean, group: ParticleGroup | null, particleIndex: number, impulse: number, normal: Vec2): void {
      if (group && isRigidGroup) {
        ///group.linearVelocity += impulse * invMass * normal;
        group.linearVelocity.SelfMulAdd(impulse * invMass, normal);
        ///group.angularVelocity += impulse * tangentDistance * invInertia;
        group.angularVelocity += impulse * tangentDistance * invInertia;
      } else {
        ///velocityBuffer.data[particleIndex] += impulse * invMass * normal;
        this.velocityBuffer.data[particleIndex].SelfMulAdd(impulse * invMass, normal);
      }
    }
  }

  export class ParticleSysteUserOverridableBuffer<T> {
    public _data: T[] | null = null;
    public get data(): T[] { return this._data as T[]; } // HACK: may return null
    public set data(value: T[]) { this._data = value; }
    public userSuppliedCapacity: number = 0;
  }

  export class ParticleSysteProxy {
    public index: number = invalidParticleIndex;
    public tag: number = 0;
    public static CompareProxyProxy(a: ParticleSysteProxy, b: ParticleSysteProxy): boolean {
      return a.tag < b.tag;
    }
    public static CompareTagProxy(a: number, b: ParticleSysteProxy): boolean {
      return a < b.tag;
    }
    public static CompareProxyTag(a: ParticleSysteProxy, b: number): boolean {
      return a.tag < b;
    }
  }

  export class ParticleSysteInsideBoundsEnumerator {
    public system: ParticleSystem;
    public xLower: number;
    public xUpper: number;
    public yLower: number;
    public yUpper: number;
    public first: number;
    public last: number;
    /**
     * InsideBoundsEnumerator enumerates all particles inside the
     * given bounds.
     *
     * Construct an enumerator with bounds of tags and a range of
     * proxies.
     */
    constructor(system: ParticleSystem, lower: number, upper: number, first: number, last: number) {
      this.system = system;
      this.xLower = (lower & ParticleSystem.xMask) >>> 0;
      this.xUpper = (upper & ParticleSystem.xMask) >>> 0;
      this.yLower = (lower & ParticleSystem.yMask) >>> 0;
      this.yUpper = (upper & ParticleSystem.yMask) >>> 0;
      this.first = first;
      this.last = last;
      // DEBUG: Assert(this.first <= this.last);
    }

    /**
     * Get index of the next particle. Returns
     * invalidParticleIndex if there are no more particles.
     */
    public GetNext(): number {
      while (this.first < this.last) {
        const xTag = (this.system.proxyBuffer.data[this.first].tag & ParticleSystem.xMask) >>> 0;
        // #if ASSERT_ENABLED
        // DEBUG: const yTag = (this.system.proxyBuffer.data[this.first].tag & ParticleSysteyMask) >>> 0;
        // DEBUG: Assert(yTag >= this.yLower);
        // DEBUG: Assert(yTag <= this.yUpper);
        // #endif
        if (xTag >= this.xLower && xTag <= this.xUpper) {
          return (this.system.proxyBuffer.data[this.first++]).index;
        }
        this.first++;
      }
      return invalidParticleIndex;
    }
  }

  export class ParticleSysteParticleListNode {
    /**
     * The head of the list.
     */
    public list!: ParticleSysteParticleListNode;
    /**
     * The next node in the list.
     */
    public next: ParticleSysteParticleListNode | null = null;
    /**
     * Number of entries in the list. Valid only for the node at the
     * head of the list.
     */
    public count: number = 0;
    /**
     * Particle index.
     */
    public index: number = 0;
  }

  /**
   * @constructor
   */
  export class ParticleSysteFixedSetAllocator<T> {
    public Allocate(itemSize: number, count: number): number {
      // TODO
      return count;
    }

    public Clear(): void {
      // TODO
    }

    public GetCount(): number {
      // TODO
      return 0;
    }

    public Invalidate(itemIndex: number): void {
      // TODO
    }

    public GetValidBuffer(): boolean[] {
      // TODO
      return [];
    }

    public GetBuffer(): T[] {
      // TODO
      return [];
    }

    public SetCount(count: number): void {
      // TODO
    }
  }

  export class ParticleSysteFixtureParticle {
    public first: Fixture;
    public second: number = invalidParticleIndex;
    constructor(fixture: Fixture, particle: number) {
      this.first = fixture;
      this.second = particle;
    }
  }

  export class ParticleSysteFixtureParticleSet extends ParticleSysteFixedSetAllocator<ParticleSysteFixtureParticle> {
    public Initialize(bodyContactBuffer: GrowableBuffer<ParticleBodyContact>, flagsBuffer: ParticleSysteUserOverridableBuffer<ParticleFlag>): void {
      // TODO
    }
    public Find(pair: ParticleSysteFixtureParticle): number {
      // TODO
      return invalidParticleIndex;
    }
  }

  export class ParticleSysteParticlePair {
    public first: number = invalidParticleIndex;
    public second: number = invalidParticleIndex;
    constructor(particleA: number, particleB: number) {
      this.first = particleA;
      this.second = particleB;
    }
  }

  export class ParticlePairSet extends ParticleSysteFixedSetAllocator<ParticleSysteParticlePair> {
    public Initialize(contactBuffer: GrowableBuffer<ParticleContact>, flagsBuffer: ParticleSysteUserOverridableBuffer<ParticleFlag>): void {
      // TODO
    }

    public Find(pair: ParticleSysteParticlePair): number {
      // TODO
      return invalidParticleIndex;
    }
  }

  export class ParticleSysteConnectionFilter {
    /**
     * Is the particle necessary for connection?
     * A pair or a triad should contain at least one 'necessary'
     * particle.
     */
    public IsNecessary(index: number): boolean {
      return true;
    }

    /**
     * An additional condition for creating a pair.
     */
    public ShouldCreatePair(a: number, b: number): boolean {
      return true;
    }

    /**
     * An additional condition for creating a triad.
     */
    public ShouldCreateTriad(a: number, b: number, c: number): boolean {
      return true;
    }
  }

  export class ParticleSysteDestroyParticlesInShapeCallback extends QueryCallback {
    public system: ParticleSystem;
    public shape: Shape;
    public xf: Transform;
    public callDestructionListener: boolean = false;
    public destroyed: number = 0;

    constructor(system: ParticleSystem, shape: Shape, xf: Transform, callDestructionListener: boolean) {
      super();
      this.system = system;
      this.shape = shape;
      this.xf = xf;
      this.callDestructionListener = callDestructionListener;
      this.destroyed = 0;
    }

    public ReportFixture(fixture: Fixture): boolean {
      return false;
    }

    public ReportParticle(particleSystem: ParticleSystem, index: number): boolean {
      if (particleSystem !== this.system) {
        return false;
      }
      // DEBUG: Assert(index >= 0 && index < this.system.count);
      if (this.shape.TestPoint(this.xf, this.system.positionBuffer.data[index])) {
        this.system.DestroyParticle(index, this.callDestructionListener);
        this.destroyed++;
      }
      return true;
    }

    public Destroyed(): number {
      return this.destroyed;
    }
  }

  export class ParticleSysteJoinParticleGroupsFilter extends ParticleSysteConnectionFilter {
    public threshold: number = 0;

    constructor(threshold: number) {
      super();
      this.threshold = threshold;
    }

    /**
     * An additional condition for creating a pair.
     */
    public ShouldCreatePair(a: number, b: number): boolean {
      return (a < this.threshold && this.threshold <= b) ||
        (b < this.threshold && this.threshold <= a);
    }

    /**
     * An additional condition for creating a triad.
     */
    public ShouldCreateTriad(a: number, b: number, c: number): boolean {
      return (a < this.threshold || b < this.threshold || c < this.threshold) &&
        (this.threshold <= a || this.threshold <= b || this.threshold <= c);
    }
  }

  export class ParticleSysteCompositeShape extends Shape {
    constructor(shapes: Shape[], shapeCount: number = shapes.length) {
      super(ShapeType.e_unknown, 0);
      this.shapes = shapes;
      this.shapeCount = shapeCount;
    }

    public shapes: Shape[];
    public shapeCount: number = 0;

    public Clone(): Shape {
      // DEBUG: Assert(false);
      throw new Error();
    }

    public GetChildCount(): number {
      return 1;
    }

    /**
     * @see Shape::TestPoint
     */
    public TestPoint(xf: Transform, p: XY): boolean {
      for (let i = 0; i < this.shapeCount; i++) {
        if (this.shapes[i].TestPoint(xf, p)) {
          return true;
        }
      }
      return false;
    }

    /**
     * @see Shape::ComputeDistance
     */
    public ComputeDistance(xf: Transform, p: Vec2, normal: Vec2, childIndex: number): number {
      // DEBUG: Assert(false);
      return 0;
    }

    /**
     * Implement Shape.
     */
    public RayCast(output: RayCastOutput, input: RayCastInput, xf: Transform, childIndex: number): boolean {
      // DEBUG: Assert(false);
      return false;
    }

    /**
     * @see Shape::ComputeAABB
     */
    public ComputeAABB(aabb: AABB, xf: Transform, childIndex: number): void {
      const s_subaabb = new AABB();
      aabb.lowerBound.x = +maxFloat;
      aabb.lowerBound.y = +maxFloat;
      aabb.upperBound.x = -maxFloat;
      aabb.upperBound.y = -maxFloat;
      // DEBUG: Assert(childIndex === 0);
      for (let i = 0; i < this.shapeCount; i++) {
        const childCount = this.shapes[i].GetChildCount();
        for (let j = 0; j < childCount; j++) {
          const subaabb = s_subaabb;
          this.shapes[i].ComputeAABB(subaabb, xf, j);
          aabb.Combine1(subaabb);
        }
      }
    }

    /**
     * @see Shape::ComputeMass
     */
    public ComputeMass(massData: MassData, density: number): void {
      // DEBUG: Assert(false);
    }

    public SetupDistanceProxy(proxy: DistanceProxy, index: number): void {
      // DEBUG: Assert(false);
    }

    public ComputeSubmergedArea(normal: Vec2, offset: number, xf: Transform, c: Vec2): number {
      // DEBUG: Assert(false);
      return 0;
    }

    public Dump(log: (format: string, ...args: any[]) => void): void {
      // DEBUG: Assert(false);
    }
  }

  export class ParticleSysteReactiveFilter extends ParticleSysteConnectionFilter {
    public flagsBuffer: ParticleSysteUserOverridableBuffer<ParticleFlag>;
    constructor(flagsBuffer: ParticleSysteUserOverridableBuffer<ParticleFlag>) {
      super();
      this.flagsBuffer = flagsBuffer;
    }
    public IsNecessary(index: number): boolean {
      return (this.flagsBuffer.data[index] & ParticleFlag.reactiveParticle) !== 0;
    }
  }

  export class ParticleSysteUpdateBodyContactsCallback extends FixtureParticleQueryCallback {
    public contactFilter: ContactFilter | null = null;
    constructor(system: ParticleSystem, contactFilter: ContactFilter | null = null) {
      super(system); // base class constructor
      this.contactFilter = contactFilter;
    }

    public ShouldCollideFixtureParticle(fixture: Fixture, particleSystem: ParticleSystem, particleIndex: number): boolean {
      // Call the contact filter if it's set, to determine whether to
      // filter this contact.  Returns true if contact calculations should
      // be performed, false otherwise.
      if (this.contactFilter) {
        const flags = this.system.GetFlagsBuffer();
        if (flags[particleIndex] & ParticleFlag.fixtureContactFilterParticle) {
          return this.contactFilter.ShouldCollideFixtureParticle(fixture, this.system, particleIndex);
        }
      }
      return true;
    }

    public ReportFixtureAndParticle(fixture: Fixture, childIndex: number, a: number): void {
      const s_n = ParticleSysteUpdateBodyContactsCallback.ReportFixtureAndParticle_s_n;
      const s_rp = ParticleSysteUpdateBodyContactsCallback.ReportFixtureAndParticle_s_rp;
      const ap = this.system.positionBuffer.data[a];
      const n = s_n;
      const d = fixture.ComputeDistance(ap, n, childIndex);
      if (d < this.system.particleDiameter && this.ShouldCollideFixtureParticle(fixture, this.system, a)) {
        const b = fixture.GetBody();
        const bp = b.GetWorldCenter();
        const bm = b.GetMass();
        const bI = b.GetInertia() - bm * b.GetLocalCenter().LengthSquared();
        const invBm = bm > 0 ? 1 / bm : 0;
        const invBI = bI > 0 ? 1 / bI : 0;
        const invAm =
          this.system.flagsBuffer.data[a] &
          ParticleFlag.wallParticle ? 0 : this.system.GetParticleInvMass();
        ///Vec2 rp = ap - bp;
        const rp = Vec2.SubVV(ap, bp, s_rp);
        const rpn = Vec2.CrossVV(rp, n);
        const invM = invAm + invBm + invBI * rpn * rpn;

        ///ParticleBodyContact& contact = system.bodyContactBuffer.Append();
        const contact = this.system.bodyContactBuffer.data[this.system.bodyContactBuffer.Append()];
        contact.index = a;
        contact.body = b;
        contact.fixture = fixture;
        contact.weight = 1 - d * this.system.inverseDiameter;
        ///contact.normal = -n;
        contact.normal.Copy(n.SelfNeg());
        contact.mass = invM > 0 ? 1 / invM : 0;
        this.system.DetectStuckParticle(a);
      }
    }
    public static readonly ReportFixtureAndParticle_s_n = new Vec2();
    public static readonly ReportFixtureAndParticle_s_rp = new Vec2();
  }

  export class ParticleSysteSolveCollisionCallback extends FixtureParticleQueryCallback {
    public step: TimeStep;
    constructor(system: ParticleSystem, step: TimeStep) {
      super(system); // base class constructor
      this.step = step;
    }

    public ReportFixtureAndParticle(fixture: Fixture, childIndex: number, a: number): void {
      const s_p1 = ParticleSysteSolveCollisionCallback.ReportFixtureAndParticle_s_p1;
      const s_output = ParticleSysteSolveCollisionCallback.ReportFixtureAndParticle_s_output;
      const s_input = ParticleSysteSolveCollisionCallback.ReportFixtureAndParticle_s_input;
      const s_p = ParticleSysteSolveCollisionCallback.ReportFixtureAndParticle_s_p;
      const s_v = ParticleSysteSolveCollisionCallback.ReportFixtureAndParticle_s_v;
      const s_f = ParticleSysteSolveCollisionCallback.ReportFixtureAndParticle_s_f;

      const body = fixture.GetBody();
      const ap = this.system.positionBuffer.data[a];
      const av = this.system.velocityBuffer.data[a];
      const output = s_output;
      const input = s_input;
      if (this.system.iterationIndex === 0) {
        // Put 'ap' in the local space of the previous frame
        ///Vec2 p1 = MulT(body.xf0, ap);
        const p1 = Transform.MulTXV(body.xf0, ap, s_p1);
        if (fixture.GetShape().GetType() === ShapeType.e_circleShape) {
          // Make relative to the center of the circle
          ///p1 -= body.GetLocalCenter();
          p1.SelfSub(body.GetLocalCenter());
          // Re-apply rotation about the center of the circle
          ///p1 = Mul(body.xf0.q, p1);
          Rot.MulRV(body.xf0.q, p1, p1);
          // Subtract rotation of the current frame
          ///p1 = MulT(body.xf.q, p1);
          Rot.MulTRV(body.xf.q, p1, p1);
          // Return to local space
          ///p1 += body.GetLocalCenter();
          p1.SelfAdd(body.GetLocalCenter());
        }
        // Return to global space and apply rotation of current frame
        ///input.p1 = Mul(body.xf, p1);
        Transform.MulXV(body.xf, p1, input.p1);
      } else {
        ///input.p1 = ap;
        input.p1.Copy(ap);
      }
      ///input.p2 = ap + step.dt * av;
      Vec2.AddVMulSV(ap, this.step.dt, av, input.p2);
      input.maxFraction = 1;
      if (fixture.RayCast(output, input, childIndex)) {
        const n = output.normal;
        ///Vec2 p = (1 - output.fraction) * input.p1 + output.fraction * input.p2 + linearSlop * n;
        const p = s_p;
        p.x = (1 - output.fraction) * input.p1.x + output.fraction * input.p2.x + linearSlop * n.x;
        p.y = (1 - output.fraction) * input.p1.y + output.fraction * input.p2.y + linearSlop * n.y;
        ///Vec2 v = step.inv_dt * (p - ap);
        const v = s_v;
        v.x = this.step.inv_dt * (p.x - ap.x);
        v.y = this.step.inv_dt * (p.y - ap.y);
        ///system.velocityBuffer.data[a] = v;
        this.system.velocityBuffer.data[a].Copy(v);
        ///Vec2 f = step.inv_dt * system.GetParticleMass() * (av - v);
        const f = s_f;
        f.x = this.step.inv_dt * this.system.GetParticleMass() * (av.x - v.x);
        f.y = this.step.inv_dt * this.system.GetParticleMass() * (av.y - v.y);
        this.system.ParticleApplyForce(a, f);
      }
    }
    public static readonly ReportFixtureAndParticle_s_p1 = new Vec2();
    public static readonly ReportFixtureAndParticle_s_output = new RayCastOutput();
    public static readonly ReportFixtureAndParticle_s_input = new RayCastInput();
    public static readonly ReportFixtureAndParticle_s_p = new Vec2();
    public static readonly ReportFixtureAndParticle_s_v = new Vec2();
    public static readonly ReportFixtureAndParticle_s_f = new Vec2();

    public ReportParticle(system: ParticleSystem, index: number): boolean {
      return false;
    }
  }

}
