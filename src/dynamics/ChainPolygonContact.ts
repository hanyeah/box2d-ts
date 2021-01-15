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
  export class ChainAndPolygonContact extends Contact<ChainShape, PolygonShape> {
    public static Create(): Contact {
      return new ChainAndPolygonContact();
    }

    public static Destroy(contact: Contact): void {
    }

    private static Evaluate_s_edge = new EdgeShape();
    public Evaluate(manifold: Manifold, xfA: Transform, xfB: Transform): void {
      const edge: EdgeShape = ChainAndPolygonContact.Evaluate_s_edge;
      this.GetShapeA().GetChildEdge(edge, this.indexA);
      CollideEdgeAndPolygon(manifold, edge, xfA, this.GetShapeB(), xfB);
    }
  }
}
