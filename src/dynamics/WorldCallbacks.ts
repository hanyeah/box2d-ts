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
/// Joints and fixtures are destroyed when their associated
/// body is destroyed. Implement this listener so that you
/// may nullify references to these joints and shapes.
namespace b2 {
  export class DestructionListener {
    /// Called when any joint is about to be destroyed due
    /// to the destruction of one of its attached bodies.
    public sayGoodbyeJoint(joint: Joint): void {}

    /// Called when any fixture is about to be destroyed due
    /// to the destruction of its parent body.
    public sayGoodbyeFixture(fixture: Fixture): void {}

    // #if ENABLE_PARTICLE
    /// Called when any particle group is about to be destroyed.
    public sayGoodbyeParticleGroup(group: ParticleGroup): void {}

    /// Called when a particle is about to be destroyed.
    /// The index can be used in conjunction with
    /// ParticleSystem::GetUserDataBuffer() or
    /// ParticleSystem::GetParticleHandleFromIndex() to determine which
    /// particle has been destroyed.
    public sayGoodbyeParticle(system: ParticleSystem, index: number): void {}
    // #endif
  }

/// Implement this class to provide collision filtering. In other words, you can implement
/// this class if you want finer control over contact creation.
  export class ContactFilter {
    /// Return true if contact calculations should be performed between these two shapes.
    /// @warning for performance reasons this is only called when the AABBs begin to overlap.
    public ShouldCollide(fixtureA: Fixture, fixtureB: Fixture): boolean {
      const bodyA: Body = fixtureA.getBody();
      const bodyB: Body = fixtureB.getBody();

      // At least one body should be dynamic or kinematic.
      if (bodyB.getType() === BodyType.StaticBody && bodyA.getType() === BodyType.StaticBody) {
        return false;
      }

      // Does a joint prevent collision?
      if (!bodyB.shouldCollideConnected(bodyA)) {
        return false;
      }

      const filter1: Filter = fixtureA.getFilterData();
      const filter2: Filter = fixtureB.getFilterData();

      if (filter1.groupIndex === filter2.groupIndex && filter1.groupIndex !== 0) {
        return (filter1.groupIndex > 0);
      }

      const collide: boolean = (((filter1.maskBits & filter2.categoryBits) !== 0) && ((filter1.categoryBits & filter2.maskBits) !== 0));
      return collide;
    }

    // #if ENABLE_PARTICLE
    public shouldCollideFixtureParticle(fixture: Fixture, system: ParticleSystem, index: number): boolean {
      return true;
    }

    public shouldCollideParticleParticle(system: ParticleSystem, indexA: number, indexB: number): boolean {
      return true;
    }
    // #endif

    public static readonly defaultFilter: ContactFilter = new ContactFilter();
  }

/// Contact impulses for reporting. Impulses are used instead of forces because
/// sub-step forces may approach infinity for rigid body collisions. These
/// match up one-to-one with the contact points in Manifold.
  export class ContactImpulse {
    public normalImpulses: number[] = MakeNumberArray(maxManifoldPoints);
    public tangentImpulses: number[] = MakeNumberArray(maxManifoldPoints);
    public count: number = 0;
  }

/// Implement this class to get contact information. You can use these results for
/// things like sounds and game logic. You can also get contact results by
/// traversing the contact lists after the time step. However, you might miss
/// some contacts because continuous physics leads to sub-stepping.
/// Additionally you may receive multiple callbacks for the same contact in a
/// single time step.
/// You should strive to make your callbacks efficient because there may be
/// many callbacks per time step.
/// @warning You cannot create/destroy Box2D entities inside these callbacks.
  export class ContactListener {
    /// Called when two fixtures begin to touch.
    public beginContact(contact: Contact): void {}

    /// Called when two fixtures cease to touch.
    public endContact(contact: Contact): void {}

    // #if ENABLE_PARTICLE
    public beginContactFixtureParticle(system: ParticleSystem, contact: ParticleBodyContact): void {}
    public endContactFixtureParticle(system: ParticleSystem, contact: ParticleBodyContact): void {}
    public beginContactParticleParticle(system: ParticleSystem, contact: ParticleContact): void {}
    public endContactParticleParticle(system: ParticleSystem, contact: ParticleContact): void {}
    // #endif

    /// This is called after a contact is updated. This allows you to inspect a
    /// contact before it goes to the solver. If you are careful, you can modify the
    /// contact manifold (e.g. disable contact).
    /// A copy of the old manifold is provided so that you can detect changes.
    /// Note: this is called only for awake bodies.
    /// Note: this is called even when the number of contact points is zero.
    /// Note: this is not called for sensors.
    /// Note: if you set the number of contact points to zero, you will not
    /// get an EndContact callback. However, you may get a BeginContact callback
    /// the next step.
    public preSolve(contact: Contact, oldManifold: Manifold): void {}

    /// This lets you inspect a contact after the solver is finished. This is useful
    /// for inspecting impulses.
    /// Note: the contact manifold does not include time of impact impulses, which can be
    /// arbitrarily large if the sub-step is small. Hence the impulse is provided explicitly
    /// in a separate data structure.
    /// Note: this is only called for contacts that are touching, solid, and awake.
    public postSolve(contact: Contact, impulse: ContactImpulse): void {}

    public static readonly defaultListener: ContactListener = new ContactListener();
  }

/// Callback class for AABB queries.
/// See World::Query
  export class QueryCallback {
    /// Called for each fixture found in the query AABB.
    /// @return false to terminate the query.
    public reportFixture(fixture: Fixture): boolean {
      return true;
    }

    // #if ENABLE_PARTICLE
    public reportParticle(system: ParticleSystem, index: number): boolean {
      return false;
    }
    public shouldQueryParticleSystem(system: ParticleSystem): boolean {
      return true;
    }
    // #endif
  }

  export type QueryCallbackFunction = (fixture: Fixture) => boolean;

/// Callback class for ray casts.
/// See World::RayCast
  export class RayCastCallback {
    /// Called for each fixture found in the query. You control how the ray cast
    /// proceeds by returning a float:
    /// return -1: ignore this fixture and continue
    /// return 0: terminate the ray cast
    /// return fraction: clip the ray to this point
    /// return 1: don't clip the ray and continue
    /// @param fixture the fixture hit by the ray
    /// @param point the point of initial intersection
    /// @param normal the normal vector at the point of intersection
    /// @return -1 to filter, 0 to terminate, fraction to clip the ray for
    /// closest hit, 1 to continue
    public reportFixture(fixture: Fixture, point: Vec2, normal: Vec2, fraction: number): number {
      return fraction;
    }

    // #if ENABLE_PARTICLE
    public reportParticle(system: ParticleSystem, index: number, point: Vec2, normal: Vec2, fraction: number): number {
      return 0;
    }
    public shouldQueryParticleSystem(system: ParticleSystem): boolean {
      return true;
    }
    // #endif
  }

  export type RayCastCallbackFunction = (fixture: Fixture, point: Vec2, normal: Vec2, fraction: number) => number;

}
