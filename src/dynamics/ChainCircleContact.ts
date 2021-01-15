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
  export class ChainAndCircleContact extends Contact<ChainShape, CircleShape> {
    public static create(): Contact {
      return new ChainAndCircleContact();
    }

    public static destroy(contact: Contact): void {
    }

    private static evaluate_s_edge = new EdgeShape();
    public evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void {
      const edge: EdgeShape = ChainAndCircleContact.evaluate_s_edge;
      this.getShapeA().getChildEdge(edge, this.indexA);
      collideEdgeAndCircle(manifold, edge, xfA, this.getShapeB(), xfB);
    }
  }

}
