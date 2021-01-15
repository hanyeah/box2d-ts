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

  /**
   * The particle type. Can be combined with the | operator.
   */
  export enum ParticleFlag {
    /// Water particle.
    WaterParticle = 0,
    /// Removed after next simulation step.
    ZombieParticle = 1 << 1,
    /// Zero velocity.
    WallParticle = 1 << 2,
    /// With restitution from stretching.
    SpringParticle = 1 << 3,
    /// With restitution from deformation.
    ElasticParticle = 1 << 4,
    /// With viscosity.
    ViscousParticle = 1 << 5,
    /// Without isotropic pressure.
    PowderParticle = 1 << 6,
    /// With surface tension.
    TensileParticle = 1 << 7,
    /// Mix color between contacting particles.
    ColorMixingParticle = 1 << 8,
    /// Call DestructionListener on destruction.
    DestructionListenerParticle = 1 << 9,
    /// Prevents other particles from leaking.
    BarrierParticle = 1 << 10,
    /// Less compressibility.
    StaticPressureParticle = 1 << 11,
    /// Makes pairs or triads with other particles.
    ReactiveParticle = 1 << 12,
    /// With high repulsive force.
    RepulsiveParticle = 1 << 13,
    /// Call ContactListener when this particle is about to interact with
    /// a rigid body or stops interacting with a rigid body.
    /// This results in an expensive operation compared to using
    /// fixtureContactFilterParticle to detect collisions between
    /// particles.
    FixtureContactListenerParticle = 1 << 14,
    /// Call ContactListener when this particle is about to interact with
    /// another particle or stops interacting with another particle.
    /// This results in an expensive operation compared to using
    /// particleContactFilterParticle to detect collisions between
    /// particles.
    ParticleContactListenerParticle = 1 << 15,
    /// Call ContactFilter when this particle interacts with rigid bodies.
    FixtureContactFilterParticle = 1 << 16,
    /// Call ContactFilter when this particle interacts with other
    /// particles.
    ParticleContactFilterParticle = 1 << 17,
  }

  export interface IParticleDef {
    flags?: ParticleFlag;
    position?: XY;
    velocity?: XY;
    color?: RGBA;
    lifetime?: number;
    userData?: any;
    group?: ParticleGroup;
  }

  export class ParticleDef implements IParticleDef {
    public flags: ParticleFlag = 0;
    public readonly position: Vec2 = new Vec2();
    public readonly velocity: Vec2 = new Vec2();
    public readonly color: Color = new Color(0, 0, 0, 0);
    public lifetime: number = 0.0;
    public userData: any = null;
    public group: ParticleGroup = null;
  }

  export function CalculateParticleIterations(gravity: number, radius: number, timeStep: number): number {
    // In some situations you may want more particle iterations than this,
    // but to avoid excessive cycle cost, don't recommend more than this.
    const MAX_RECOMMENDED_PARTICLE_ITERATIONS = 8;
    const RADIUS_THRESHOLD = 0.01;
    const iterations = Math.ceil(Math.sqrt(gravity / (RADIUS_THRESHOLD * radius)) * timeStep);
    return clamp(iterations, 1, MAX_RECOMMENDED_PARTICLE_ITERATIONS);
  }

  export class ParticleHandle {
    public index: number = invalidParticleIndex;
    public getIndex(): number { return this.index; }
    public setIndex(index: number): void { this.index = index; }
  }

}

