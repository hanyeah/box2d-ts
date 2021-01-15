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
  export class Pair<T> {
    constructor(public proxyA: TreeNode<T>, public proxyB: TreeNode<T>) { }
  }

/// The broad-phase is used for computing pairs and performing volume queries and ray casts.
/// This broad-phase does not persist pairs. Instead, this reports potentially new pairs.
/// It is up to the client to consume the new pairs and to track subsequent overlap.
  export class BroadPhase<T> {
    public readonly tree: DynamicTree<T> = new DynamicTree<T>();
    public proxyCount: number = 0;
    // public moveCapacity: number = 16;
    public moveCount: number = 0;
    public readonly moveBuffer: Array<TreeNode<T>> = [];
    // public pairCapacity: number = 16;
    public pairCount: number = 0;
    public readonly pairBuffer: Array<Pair<T>> = [];
    // public queryProxyId: number = 0;

    /// Create a proxy with an initial AABB. Pairs are not reported until
    /// UpdatePairs is called.
    public createProxy(aabb: AABB, userData: T): TreeNode<T> {
      const proxy: TreeNode<T> = this.tree.createProxy(aabb, userData);
      ++this.proxyCount;
      this.bufferMove(proxy);
      return proxy;
    }

    /// Destroy a proxy. It is up to the client to remove any pairs.
    public destroyProxy(proxy: TreeNode<T>): void {
      this.unBufferMove(proxy);
      --this.proxyCount;
      this.tree.destroyProxy(proxy);
    }

    /// Call MoveProxy as many times as you like, then when you are done
    /// call UpdatePairs to finalized the proxy pairs (for your time step).
    public moveProxy(proxy: TreeNode<T>, aabb: AABB, displacement: Vec2): void {
      const buffer: boolean = this.tree.moveProxy(proxy, aabb, displacement);
      if (buffer) {
        this.bufferMove(proxy);
      }
    }

    /// Call to trigger a re-processing of it's pairs on the next call to UpdatePairs.
    public touchProxy(proxy: TreeNode<T>): void {
      this.bufferMove(proxy);
    }

    /// Get the fat AABB for a proxy.
    // public GetFatAABB(proxy: TreeNode<T>): AABB {
    //   return this.tree.GetFatAABB(proxy);
    // }

    /// Get user data from a proxy. Returns NULL if the id is invalid.
    // public GetUserData(proxy: TreeNode<T>): T {
    //   return this.tree.GetUserData(proxy);
    // }

    /// Test overlap of fat AABBs.
    // public TestOverlap(proxyA: TreeNode<T>, proxyB: TreeNode<T>): boolean {
    //   const aabbA: AABB = this.tree.GetFatAABB(proxyA);
    //   const aabbB: AABB = this.tree.GetFatAABB(proxyB);
    //   return TestOverlapAABB(aabbA, aabbB);
    // }

    /// Get the number of proxies.
    public getProxyCount(): number {
      return this.proxyCount;
    }

    /// Update the pairs. This results in pair callbacks. This can only add pairs.
    public updatePairs(callback: (a: T, b: T) => void): void {
      // Reset pair buffer
      this.pairCount = 0;

      // Perform tree queries for all moving proxies.
      for (let i: number = 0; i < this.moveCount; ++i) {
        const queryProxy: TreeNode<T> = this.moveBuffer[i];
        if (queryProxy === null) {
          continue;
        }

        // This is called from .DynamicTree::Query when we are gathering pairs.
        // boolean BroadPhase::QueryCallback(int32 proxyId);

        // We have to query the tree with the fat AABB so that
        // we don't fail to create a pair that may touch later.
        const fatAABB: AABB = queryProxy.aabb; // this.tree.GetFatAABB(queryProxy);

        // Query tree, create pairs and add them pair buffer.
        this.tree.query(fatAABB, (proxy: TreeNode<T>): boolean => {
          // A proxy cannot form a pair with itself.
          if (proxy.id === queryProxy.id) {
            return true;
          }

          const moved: boolean = proxy.moved; // this.tree.WasMoved(proxy);
          if (moved && proxy.id > queryProxy.id) {
            // Both proxies are moving. Avoid duplicate pairs.
            return true;
          }

          // const proxyA = proxy < queryProxy ? proxy : queryProxy;
          // const proxyB = proxy >= queryProxy ? proxy : queryProxy;
          let proxyA: TreeNode<T>;
          let proxyB: TreeNode<T>;
          if (proxy.id < queryProxy.id) {
            proxyA = proxy;
            proxyB = queryProxy;
          } else {
            proxyA = queryProxy;
            proxyB = proxy;
          }

          // Grow the pair buffer as needed.
          if (this.pairCount === this.pairBuffer.length) {
            this.pairBuffer[this.pairCount] = new Pair(proxyA, proxyB);
          } else {
            const pair: Pair<T> = this.pairBuffer[this.pairCount];
            pair.proxyA = proxyA;
            pair.proxyB = proxyB;
          }

          ++this.pairCount;

          return true;
        });
      }

      // Send pairs to caller
      for (let i = 0; i < this.pairCount; ++i) {
        const primaryPair: Pair<T> = this.pairBuffer[i];
        const userDataA: T = primaryPair.proxyA.userData; // this.tree.GetUserData(primaryPair.proxyA);
        const userDataB: T = primaryPair.proxyB.userData; // this.tree.GetUserData(primaryPair.proxyB);

        callback(userDataA, userDataB);
      }

      // Clear move flags
      for (let i = 0; i < this.moveCount; ++i) {
        const proxy: TreeNode<T> = this.moveBuffer[i];
        if (proxy === null) {
          continue;
        }

        proxy.moved = false; // this.tree.ClearMoved(proxy);
      }

      // Reset move buffer
      this.moveCount = 0;
    }

    /// Query an AABB for overlapping proxies. The callback class
    /// is called for each proxy that overlaps the supplied AABB.
    public query(aabb: AABB, callback: (node: TreeNode<T>) => boolean): void {
      this.tree.query(aabb, callback);
    }

    public queryPoint(point: XY, callback: (node: TreeNode<T>) => boolean): void {
      this.tree.queryPoint(point, callback);
    }

    /// Ray-cast against the proxies in the tree. This relies on the callback
    /// to perform a exact ray-cast in the case were the proxy contains a shape.
    /// The callback also performs the any collision filtering. This has performance
    /// roughly equal to k * log(n), where k is the number of collisions and n is the
    /// number of proxies in the tree.
    /// @param input the ray-cast input data. The ray extends from p1 to p1 + maxFraction * (p2 - p1).
    /// @param callback a callback class that is called for each proxy that is hit by the ray.
    public rayCast(input: RayCastInput, callback: (input: RayCastInput, node: TreeNode<T>) => number): void {
      this.tree.rayCast(input, callback);
    }

    /// Get the height of the embedded tree.
    public getTreeHeight(): number {
      return this.tree.getHeight();
    }

    /// Get the balance of the embedded tree.
    public getTreeBalance(): number {
      return this.tree.getMaxBalance();
    }

    /// Get the quality metric of the embedded tree.
    public getTreeQuality(): number {
      return this.tree.getAreaRatio();
    }

    /// Shift the world origin. Useful for large worlds.
    /// The shift formula is: position -= newOrigin
    /// @param newOrigin the new origin with respect to the old origin
    public shiftOrigin(newOrigin: XY): void {
      this.tree.shiftOrigin(newOrigin);
    }

    public bufferMove(proxy: TreeNode<T>): void {
      this.moveBuffer[this.moveCount] = proxy;
      ++this.moveCount;
    }

    public unBufferMove(proxy: TreeNode<T>): void {
      const i: number = this.moveBuffer.indexOf(proxy);
      this.moveBuffer[i] = null;
    }
  }
}

