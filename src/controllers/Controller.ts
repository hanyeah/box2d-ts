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
   * A controller edge is used to connect bodies and controllers
   * together in a bipartite graph.
   */
  export class ControllerEdge {
    public readonly controller: Controller; ///< provides quick access to other end of this edge.
    public readonly body: Body; ///< the body
    public prevBody: ControllerEdge = null; ///< the previous controller edge in the controllers's joint list
    public nextBody: ControllerEdge = null; ///< the next controller edge in the controllers's joint list
    public prevController: ControllerEdge = null; ///< the previous controller edge in the body's joint list
    public nextController: ControllerEdge = null; ///< the next controller edge in the body's joint list
    constructor(controller: Controller, body: Body) {
      this.controller = controller;
      this.body = body;
    }
  }

  /**
   * Base class for controllers. Controllers are a convience for
   * encapsulating common per-step functionality.
   */
  export abstract class Controller {
    // world: World;
    public bodyList: ControllerEdge = null;
    public bodyCount: number = 0;
    public prev: Controller = null;
    public next: Controller = null;

    /**
     * Controllers override this to implement per-step functionality.
     */
    public abstract step(step: TimeStep): void;

    /**
     * Controllers override this to provide debug drawing.
     */
    public draw(debugDraw: Draw): void{

    }

    /**
     * Adds a body to the controller list.
     */
    public addBody(body: Body): void {
      const edge = new ControllerEdge(this, body);
      //Add edge to controller list
      edge.nextBody = this.bodyList;
      edge.prevBody = null;
      if (this.bodyList) {
        this.bodyList.prevBody = edge;
      }
      this.bodyList = edge;
      ++this.bodyCount;

      //Add edge to body list
      edge.nextController = body.controllerList;
      edge.prevController = null;
      if (body.controllerList) {
        body.controllerList.prevController = edge;
      }
      body.controllerList = edge;
      ++body.controllerCount;
    }

    /**
     * Removes a body from the controller list.
     */
    public removeBody(body: Body): void {
      //Assert that the controller is not empty
      if (this.bodyCount <= 0) { throw new Error(); }

      //Find the corresponding edge
      /*ControllerEdge*/
      let edge = this.bodyList;
      while (edge && edge.body !== body) {
        edge = edge.nextBody;
      }

      //Assert that we are removing a body that is currently attached to the controller
      if (edge === null) { throw new Error(); }

      //Remove edge from controller list
      if (edge.prevBody) {
        edge.prevBody.nextBody = edge.nextBody;
      }
      if (edge.nextBody) {
        edge.nextBody.prevBody = edge.prevBody;
      }
      if (this.bodyList === edge) {
        this.bodyList = edge.nextBody;
      }
      --this.bodyCount;

      //Remove edge from body list
      if (edge.nextController) {
        edge.nextController.prevController = edge.prevController;
      }
      if (edge.prevController) {
        edge.prevController.nextController = edge.nextController;
      }
      if (body.controllerList === edge) {
        body.controllerList = edge.nextController;
      }
      --body.controllerCount;
    }

    /**
     * Removes all bodies from the controller list.
     */
    public clear(): void {
      while (this.bodyList) {
        this.removeBody(this.bodyList.body);
      }

      this.bodyCount = 0;
    }
  }

}

