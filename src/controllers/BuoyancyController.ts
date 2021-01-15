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

  /**
   * Calculates buoyancy forces for fluids in the form of a half
   * plane.
   */
  export class BuoyancyController extends Controller {
    /**
     * The outer surface normal
     */
    public readonly normal = new Vec2(0, 1);
    /**
     * The height of the fluid surface along the normal
     */
    public offset = 0;
    /**
     * The fluid density
     */
    public density = 0;
    /**
     * Fluid velocity, for drag calculations
     */
    public readonly velocity = new Vec2(0, 0);
    /**
     * Linear drag co-efficient
     */
    public linearDrag = 0;
    /**
     * Angular drag co-efficient
     */
    public angularDrag = 0;
    /**
     * If false, bodies are assumed to be uniformly dense, otherwise
     * use the shapes densities
     */
    public useDensity = false; //False by default to prevent a gotcha
    /**
     * If true, gravity is taken from the world instead of the
     */
    public useWorldGravity = true;
    /**
     * Gravity vector, if the world's gravity is not used
     */
    public readonly gravity = new Vec2(0, 0);

    public step(step: TimeStep) {
      if (!this.bodyList) {
        return;
      }
      if (this.useWorldGravity) {
        this.gravity.copy(this.bodyList.body.getWorld().gravity);
      }
      for (let i: ControllerEdge = this.bodyList; i; i = i.nextBody) {
        const body = i.body;
        if (!body.isAwake()) {
          //Buoyancy force is just a function of position,
          //so unlike most forces, it is safe to ignore sleeping bodes
          continue;
        }
        const areac = new Vec2();
        const massc = new Vec2();
        let area = 0;
        let mass = 0;
        for (let fixture = body.getFixtureList(); fixture; fixture = fixture.next) {
          const sc = new Vec2();
          const sarea = fixture.getShape().computeSubmergedArea(this.normal, this.offset, body.getTransform(), sc);
          area += sarea;
          areac.x += sarea * sc.x;
          areac.y += sarea * sc.y;
          let shapeDensity = 0;
          if (this.useDensity) {
            //TODO: Expose density publicly
            shapeDensity = fixture.getDensity();
          } else {
            shapeDensity = 1;
          }
          mass += sarea * shapeDensity;
          massc.x += sarea * sc.x * shapeDensity;
          massc.y += sarea * sc.y * shapeDensity;
        }
        areac.x /= area;
        areac.y /= area;
        //    Vec2 localCentroid = MulT(body->GetXForm(),areac);
        massc.x /= mass;
        massc.y /= mass;
        if (area < epsilon) {
          continue;
        }
        //Buoyancy
        const buoyancyForce = this.gravity.clone().selfNeg();
        buoyancyForce.selfMul(this.density * area);
        body.applyForce(buoyancyForce, massc);
        //Linear drag
        const dragForce = body.getLinearVelocityFromWorldPoint(areac, new Vec2());
        dragForce.selfSub(this.velocity);
        dragForce.selfMul((-this.linearDrag * area));
        body.applyForce(dragForce, areac);
        //Angular drag
        //TODO: Something that makes more physical sense?
        body.applyTorque((-body.getInertia() / body.getMass() * area * body.getAngularVelocity() * this.angularDrag));
      }
    }

    public draw(debugDraw: Draw) {
      const r = 100;
      const p1 = new Vec2();
      const p2 = new Vec2();
      p1.x = this.normal.x * this.offset + this.normal.y * r;
      p1.y = this.normal.y * this.offset - this.normal.x * r;
      p2.x = this.normal.x * this.offset - this.normal.y * r;
      p2.y = this.normal.y * this.offset + this.normal.x * r;

      const color = new Color(0, 0, 0.8);

      debugDraw.drawSegment(p1, p2, color);
    }
  }

}

