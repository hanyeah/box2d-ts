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
  export enum ParticleGroupFlag {
    /// Prevents overlapping or leaking.
    SolidParticleGroup = 1 << 0,
    /// Keeps its shape.
    RigidParticleGroup = 1 << 1,
    /// Won't be destroyed if it gets empty.
    ParticleGroupCanBeEmpty = 1 << 2,
    /// Will be destroyed on next simulation step.
    ParticleGroupWillBeDestroyed = 1 << 3,
    /// Updates depth data on next simulation step.
    ParticleGroupNeedsUpdateDepth = 1 << 4,

    ParticleGroupInternalMask = ParticleGroupWillBeDestroyed | ParticleGroupNeedsUpdateDepth,
  }

  export interface IParticleGroupDef {
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

  export class ParticleGroupDef implements IParticleGroupDef {
    public flags: ParticleFlag = 0;
    public groupFlags: ParticleGroupFlag = 0;
    public readonly position: Vec2 = new Vec2();
    public angle: number = 0.0;
    public readonly linearVelocity: Vec2 = new Vec2();
    public angularVelocity: number = 0.0;
    public readonly color: Color = new Color();
    public strength: number = 1.0;
    public shape?: Shape;
    public shapes?: Shape[];
    public shapeCount: number = 0;
    public stride: number = 0;
    public particleCount: number = 0;
    public positionData?: Vec2[];
    public lifetime: number = 0;
    public userData: any = null;
    public group: ParticleGroup = null;
  }

  export class ParticleGroup {

    public readonly system: ParticleSystem;
    public firstIndex: number = 0;
    public lastIndex: number = 0;
    public groupFlags: ParticleGroupFlag = 0;
    public strength: number = 1.0;
    public prev: ParticleGroup = null;
    public next: ParticleGroup = null;
    public timestamp: number = -1;
    public mass: number = 0.0;
    public inertia: number = 0.0;
    public readonly center: Vec2 = new Vec2();
    public readonly linearVelocity: Vec2 = new Vec2();
    public angularVelocity: number = 0.0;
    public readonly transform: Transform = new Transform();
    ///transform.SetIdentity();
    public userData: any = null;

    constructor(system: ParticleSystem) {
      this.system = system;
    }

    public getNext(): ParticleGroup {
      return this.next;
    }

    public getParticleSystem(): ParticleSystem {
      return this.system;
    }

    public getParticleCount(): number {
      return this.lastIndex - this.firstIndex;
    }

    public getBufferIndex(): number {
      return this.firstIndex;
    }

    public containsParticle(index: number): boolean {
      return this.firstIndex <= index && index < this.lastIndex;
    }

    public getAllParticleFlags(): ParticleFlag {
      if (!this.system.flagsBuffer.data) { throw new Error(); }
      let flags = 0;
      for (let i = this.firstIndex; i < this.lastIndex; i++) {
        flags |= this.system.flagsBuffer.data[i];
      }
      return flags;
    }

    public getGroupFlags(): ParticleGroupFlag {
      return this.groupFlags;
    }

    public setGroupFlags(flags: number): void {
      // DEBUG: Assert((flags & ParticleGroupFlag.particleGroupInternalMask) === 0);
      flags |= this.groupFlags & ParticleGroupFlag.ParticleGroupInternalMask;
      this.system.setGroupFlags(this, flags);
    }

    public getMass(): number {
      this.updateStatistics();
      return this.mass;
    }

    public getInertia(): number {
      this.updateStatistics();
      return this.inertia;
    }

    public getCenter(): Vec2 {
      this.updateStatistics();
      return this.center;
    }

    public getLinearVelocity(): Vec2 {
      this.updateStatistics();
      return this.linearVelocity;
    }

    public getAngularVelocity(): number {
      this.updateStatistics();
      return this.angularVelocity;
    }

    public getTransform(): Transform {
      return this.transform;
    }

    public getPosition(): Vec2 {
      return this.transform.p;
    }

    public getAngle(): number {
      return this.transform.q.getAngle();
    }

    public getLinearVelocityFromWorldPoint<T extends XY>(worldPoint: XY, out: T): T {
      const s_t0 = ParticleGroup.GetLinearVelocityFromWorldPoint_s_t0;
      this.updateStatistics();
      ///  return linearVelocity + Cross(angularVelocity, worldPoint - center);
      return Vec2.AddVCrossSV(this.linearVelocity, this.angularVelocity, Vec2.SubVV(worldPoint, this.center, s_t0), out);
    }
    public static readonly GetLinearVelocityFromWorldPoint_s_t0 = new Vec2();

    public getUserData(): void {
      return this.userData;
    }

    public setUserData(data: any): void {
      this.userData = data;
    }

    public applyForce(force: XY): void {
      this.system.applyForce(this.firstIndex, this.lastIndex, force);
    }

    public applyLinearImpulse(impulse: XY): void {
      this.system.applyLinearImpulse(this.firstIndex, this.lastIndex, impulse);
    }

    public destroyParticles(callDestructionListener: boolean): void {
      if (this.system.world.isLocked()) { throw new Error(); }

      for (let i = this.firstIndex; i < this.lastIndex; i++) {
        this.system.destroyParticle(i, callDestructionListener);
      }
    }

    public updateStatistics(): void {
      if (!this.system.positionBuffer.data) { throw new Error(); }
      if (!this.system.velocityBuffer.data) { throw new Error(); }
      const p = new Vec2();
      const v = new Vec2();
      if (this.timestamp !== this.system.timestamp) {
        const m = this.system.getParticleMass();
        ///  this.mass = 0;
        this.mass = m * (this.lastIndex - this.firstIndex);
        this.center.setZero();
        this.linearVelocity.setZero();
        for (let i = this.firstIndex; i < this.lastIndex; i++) {
          ///  this.mass += m;
          ///  this.center += m * this.system.positionBuffer.data[i];
          this.center.selfMulAdd(m, this.system.positionBuffer.data[i]);
          ///  this.linearVelocity += m * this.system.velocityBuffer.data[i];
          this.linearVelocity.selfMulAdd(m, this.system.velocityBuffer.data[i]);
        }
        if (this.mass > 0) {
          const inv_mass = 1 / this.mass;
          ///this.center *= 1 / this.mass;
          this.center.selfMul(inv_mass);
          ///this.linearVelocity *= 1 / this.mass;
          this.linearVelocity.selfMul(inv_mass);
        }
        this.inertia = 0;
        this.angularVelocity = 0;
        for (let i = this.firstIndex; i < this.lastIndex; i++) {
          ///Vec2 p = this.system.positionBuffer.data[i] - this.center;
          Vec2.SubVV(this.system.positionBuffer.data[i], this.center, p);
          ///Vec2 v = this.system.velocityBuffer.data[i] - this.linearVelocity;
          Vec2.SubVV(this.system.velocityBuffer.data[i], this.linearVelocity, v);
          this.inertia += m * Vec2.DotVV(p, p);
          this.angularVelocity += m * Vec2.CrossVV(p, v);
        }
        if (this.inertia > 0) {
          this.angularVelocity *= 1 / this.inertia;
        }
        this.timestamp = this.system.timestamp;
      }
    }
  }

}

