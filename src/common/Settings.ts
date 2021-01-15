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

// Tunable Constants

/// You can use this to change the length scale used by your game.
/// For example for inches you could use 39.4.
// export const lengthUnitsPerMeter: number = 1.0;

/// The maximum number of vertices on a convex polygon. You cannot increase
/// this too much because BlockAllocator has a maximum object size.
// export const maxPolygonVertices: number = 8;

/// Logging function.
  export function log(message: string, ...args: any[]): void {
    console.log(message, ...args);
  }

}

