/*
* Copyright (c) 2009 Erin Catto http://www.box2d.org
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
  function verify<T>(value: T): T {
    if (value === null) { throw new Error(); }
    return value;
  }

/// A node in the dynamic tree. The client does not interact with this directly.
  export class TreeNode<T> {
    public readonly id: number = 0;
    public readonly aabb: AABB = new AABB();
    private _userData: T = null;
    public get userData(): T {
      if (this._userData === null) { throw new Error(); }
      return this._userData;
    }
    public set userData(value: T) {
      if (this._userData !== null) { throw new Error(); }
      this._userData = value;
    }
    public parent: TreeNode<T> = null; // or next
    public child1: TreeNode<T> = null;
    public child2: TreeNode<T> = null;
    public height: number = 0; // leaf = 0, free node = -1

    public moved: boolean = false;

    constructor(id: number = 0) {
      this.id = id;
    }

    public reset(): void {
      this._userData = null;
    }

    public isLeaf(): boolean {
      return this.child1 === null;
    }
  }

  export class DynamicTree<T> {
    public root: TreeNode<T> = null;

    // TreeNode* public nodes;
    // int32 public nodeCount;
    // int32 public nodeCapacity;

    public freeList: TreeNode<T> = null;

    public insertionCount: number = 0;

    public readonly stack = new GrowableStack<TreeNode<T>>(256);
    public static readonly sR = new Vec2();
    public static readonly sV = new Vec2();
    public static readonly absV = new Vec2();
    public static readonly segmentAABB = new AABB();
    public static readonly subInput = new RayCastInput();
    public static readonly combinedAABB = new AABB();
    public static readonly sAabb = new AABB();

    // public GetUserData(node: TreeNode<T>): T {
    //   // DEBUG: Assert(node !== null);
    //   return node.userData;
    // }

    // public WasMoved(node: TreeNode<T>): boolean {
    //   return node.moved;
    // }

    // public ClearMoved(node: TreeNode<T>): void {
    //   node.moved = false;
    // }

    // public GetFatAABB(node: TreeNode<T>): AABB {
    //   // DEBUG: Assert(node !== null);
    //   return node.aabb;
    // }

    public query(aabb: AABB, callback: (node: TreeNode<T>) => boolean): void {
      const stack: GrowableStack<TreeNode<T>> = this.stack.reset();
      stack.push(this.root);

      while (stack.getCount() > 0) {
        const node: TreeNode<T> = stack.pop();
        if (node === null) {
          continue;
        }

        if (node.aabb.testOverlap(aabb)) {
          if (node.isLeaf()) {
            const proceed: boolean = callback(node);
            if (!proceed) {
              return;
            }
          } else {
            stack.push(node.child1);
            stack.push(node.child2);
          }
        }
      }
    }

    public queryPoint(point: XY, callback: (node: TreeNode<T>) => boolean): void {
      const stack: GrowableStack<TreeNode<T>> = this.stack.reset();
      stack.push(this.root);

      while (stack.getCount() > 0) {
        const node: TreeNode<T> = stack.pop();
        if (node === null) {
          continue;
        }

        if (node.aabb.testContain(point)) {
          if (node.isLeaf()) {
            const proceed: boolean = callback(node);
            if (!proceed) {
              return;
            }
          } else {
            stack.push(node.child1);
            stack.push(node.child2);
          }
        }
      }
    }

    public rayCast(input: RayCastInput, callback: (input: RayCastInput, node: TreeNode<T>) => number): void {
      const p1: Vec2 = input.p1;
      const p2: Vec2 = input.p2;
      const r: Vec2 = Vec2.SubVV(p2, p1, DynamicTree.sR);
      // DEBUG: Assert(r.LengthSquared() > 0);
      r.normalize();

      // v is perpendicular to the segment.
      const v: Vec2 = Vec2.CrossOneV(r, DynamicTree.sV);
      const abs_v: Vec2 = Vec2.AbsV(v, DynamicTree.absV);

      // Separating axis for segment (Gino, p80).
      // |dot(v, p1 - c)| > dot(|v|, h)

      let maxFraction: number = input.maxFraction;

      // Build a bounding box for the segment.
      const segmentAABB: AABB = DynamicTree.segmentAABB;
      let t_x: number = p1.x + maxFraction * (p2.x - p1.x);
      let t_y: number = p1.y + maxFraction * (p2.y - p1.y);
      segmentAABB.lowerBound.x = Min(p1.x, t_x);
      segmentAABB.lowerBound.y = Min(p1.y, t_y);
      segmentAABB.upperBound.x = Max(p1.x, t_x);
      segmentAABB.upperBound.y = Max(p1.y, t_y);

      const stack: GrowableStack<TreeNode<T>> = this.stack.reset();
      stack.push(this.root);

      while (stack.getCount() > 0) {
        const node: TreeNode<T> = stack.pop();
        if (node === null) {
          continue;
        }

        if (!testOverlapAABB(node.aabb, segmentAABB)) {
          continue;
        }

        // Separating axis for segment (Gino, p80).
        // |dot(v, p1 - c)| > dot(|v|, h)
        const c: Vec2 = node.aabb.getCenter();
        const h: Vec2 = node.aabb.getExtents();
        const separation: number = Abs(Vec2.DotVV(v, Vec2.SubVV(p1, c, Vec2.s_t0))) - Vec2.DotVV(abs_v, h);
        if (separation > 0) {
          continue;
        }

        if (node.isLeaf()) {
          const subInput: RayCastInput = DynamicTree.subInput;
          subInput.p1.copy(input.p1);
          subInput.p2.copy(input.p2);
          subInput.maxFraction = maxFraction;

          const value: number = callback(subInput, node);

          if (value === 0) {
            // The client has terminated the ray cast.
            return;
          }

          if (value > 0) {
            // Update segment bounding box.
            maxFraction = value;
            t_x = p1.x + maxFraction * (p2.x - p1.x);
            t_y = p1.y + maxFraction * (p2.y - p1.y);
            segmentAABB.lowerBound.x = Min(p1.x, t_x);
            segmentAABB.lowerBound.y = Min(p1.y, t_y);
            segmentAABB.upperBound.x = Max(p1.x, t_x);
            segmentAABB.upperBound.y = Max(p1.y, t_y);
          }
        } else {
          stack.push(node.child1);
          stack.push(node.child2);
        }
      }
    }

    public static sNodeId: number = 0;

    public allocateNode(): TreeNode<T> {
      // Expand the node pool as needed.
      if (this.freeList !== null) {
        const node: TreeNode<T> = this.freeList;
        this.freeList = node.parent; // this.freeList = node.next;
        node.parent = null;
        node.child1 = null;
        node.child2 = null;
        node.height = 0;
        node.moved = false;
        return node;
      }

      return new TreeNode<T>(DynamicTree.sNodeId++);
    }

    public freeNode(node: TreeNode<T>): void {
      node.parent = this.freeList; // node.next = this.freeList;
      node.child1 = null;
      node.child2 = null;
      node.height = -1;
      node.reset();
      this.freeList = node;
    }

    public createProxy(aabb: AABB, userData: T): TreeNode<T> {
      const node: TreeNode<T> = this.allocateNode();

      // Fatten the aabb.
      const r_x: number = aabbExtension;
      const r_y: number = aabbExtension;
      node.aabb.lowerBound.x = aabb.lowerBound.x - r_x;
      node.aabb.lowerBound.y = aabb.lowerBound.y - r_y;
      node.aabb.upperBound.x = aabb.upperBound.x + r_x;
      node.aabb.upperBound.y = aabb.upperBound.y + r_y;
      node.userData = userData;
      node.height = 0;
      node.moved = true;

      this.insertLeaf(node);

      return node;
    }

    public destroyProxy(node: TreeNode<T>): void {
      // DEBUG: Assert(node.IsLeaf());

      this.removeLeaf(node);
      this.freeNode(node);
    }

    private static moveProxy_s_fatAABB = new AABB();
    private static moveProxy_s_hugeAABB = new AABB();
    public moveProxy(node: TreeNode<T>, aabb: AABB, displacement: Vec2): boolean {
      // DEBUG: Assert(node.IsLeaf());

      // Extend AABB
      const fatAABB: AABB = DynamicTree.moveProxy_s_fatAABB;
      const r_x: number = aabbExtension;
      const r_y: number = aabbExtension;
      fatAABB.lowerBound.x = aabb.lowerBound.x - r_x;
      fatAABB.lowerBound.y = aabb.lowerBound.y - r_y;
      fatAABB.upperBound.x = aabb.upperBound.x + r_x;
      fatAABB.upperBound.y = aabb.upperBound.y + r_y;

      // Predict AABB movement
      const d_x: number = aabbMultiplier * displacement.x;
      const d_y: number = aabbMultiplier * displacement.y;

      if (d_x < 0.0) {
        fatAABB.lowerBound.x += d_x;
      } else {
        fatAABB.upperBound.x += d_x;
      }

      if (d_y < 0.0) {
        fatAABB.lowerBound.y += d_y;
      } else {
        fatAABB.upperBound.y += d_y;
      }

      const treeAABB = node.aabb; // nodes[proxyId].aabb;
      if (treeAABB.contains(aabb)) {
        // The tree AABB still contains the object, but it might be too large.
        // Perhaps the object was moving fast but has since gone to sleep.
        // The huge AABB is larger than the new fat AABB.
        const hugeAABB: AABB = DynamicTree.moveProxy_s_hugeAABB;
        hugeAABB.lowerBound.x = fatAABB.lowerBound.x - 4.0 * r_x;
        hugeAABB.lowerBound.y = fatAABB.lowerBound.y - 4.0 * r_y;
        hugeAABB.upperBound.x = fatAABB.upperBound.x + 4.0 * r_x;
        hugeAABB.upperBound.y = fatAABB.upperBound.y + 4.0 * r_y;

        if (hugeAABB.contains(treeAABB)) {
          // The tree AABB contains the object AABB and the tree AABB is
          // not too large. No tree update needed.
          return false;
        }

        // Otherwise the tree AABB is huge and needs to be shrunk
      }

      this.removeLeaf(node);

      node.aabb.copy(fatAABB); // nodes[proxyId].aabb = fatAABB;

      this.insertLeaf(node);

      node.moved = true;

      return true;
    }

    public insertLeaf(leaf: TreeNode<T>): void {
      ++this.insertionCount;

      if (this.root === null) {
        this.root = leaf;
        this.root.parent = null;
        return;
      }

      // Find the best sibling for this node
      const leafAABB: AABB = leaf.aabb;
      let sibling: TreeNode<T> = this.root;
      while (!sibling.isLeaf()) {
        const child1: TreeNode<T> = verify(sibling.child1);
        const child2: TreeNode<T> = verify(sibling.child2);

        const area: number = sibling.aabb.getPerimeter();

        const combinedAABB: AABB = DynamicTree.combinedAABB;
        combinedAABB.combine2(sibling.aabb, leafAABB);
        const combinedArea: number = combinedAABB.getPerimeter();

        // Cost of creating a new parent for this node and the new leaf
        const cost: number = 2 * combinedArea;

        // Minimum cost of pushing the leaf further down the tree
        const inheritanceCost: number = 2 * (combinedArea - area);

        // Cost of descending into child1
        let cost1: number;
        const aabb: AABB = DynamicTree.sAabb;
        let oldArea: number;
        let newArea: number;
        if (child1.isLeaf()) {
          aabb.combine2(leafAABB, child1.aabb);
          cost1 = aabb.getPerimeter() + inheritanceCost;
        } else {
          aabb.combine2(leafAABB, child1.aabb);
          oldArea = child1.aabb.getPerimeter();
          newArea = aabb.getPerimeter();
          cost1 = (newArea - oldArea) + inheritanceCost;
        }

        // Cost of descending into child2
        let cost2: number;
        if (child2.isLeaf()) {
          aabb.combine2(leafAABB, child2.aabb);
          cost2 = aabb.getPerimeter() + inheritanceCost;
        } else {
          aabb.combine2(leafAABB, child2.aabb);
          oldArea = child2.aabb.getPerimeter();
          newArea = aabb.getPerimeter();
          cost2 = newArea - oldArea + inheritanceCost;
        }

        // Descend according to the minimum cost.
        if (cost < cost1 && cost < cost2) {
          break;
        }

        // Descend
        if (cost1 < cost2) {
          sibling = child1;
        } else {
          sibling = child2;
        }
      }

      // Create a parent for the siblings.
      const oldParent: TreeNode<T> = sibling.parent;
      const newParent: TreeNode<T> = this.allocateNode();
      newParent.parent = oldParent;
      newParent.aabb.combine2(leafAABB, sibling.aabb);
      newParent.height = sibling.height + 1;

      if (oldParent !== null) {
        // The sibling was not the root.
        if (oldParent.child1 === sibling) {
          oldParent.child1 = newParent;
        } else {
          oldParent.child2 = newParent;
        }

        newParent.child1 = sibling;
        newParent.child2 = leaf;
        sibling.parent = newParent;
        leaf.parent = newParent;
      } else {
        // The sibling was the root.
        newParent.child1 = sibling;
        newParent.child2 = leaf;
        sibling.parent = newParent;
        leaf.parent = newParent;
        this.root = newParent;
      }

      // Walk back up the tree fixing heights and AABBs
      let node: TreeNode<T> = leaf.parent;
      while (node !== null) {
        node = this.balance(node);

        const child1: TreeNode<T> = verify(node.child1);
        const child2: TreeNode<T> = verify(node.child2);

        node.height = 1 + Max(child1.height, child2.height);
        node.aabb.combine2(child1.aabb, child2.aabb);

        node = node.parent;
      }

      // this.Validate();
    }

    public removeLeaf(leaf: TreeNode<T>): void {
      if (leaf === this.root) {
        this.root = null;
        return;
      }

      const parent: TreeNode<T> = verify(leaf.parent);
      const grandParent: TreeNode<T> = parent && parent.parent;
      const sibling: TreeNode<T> = verify(parent.child1 === leaf ? parent.child2 : parent.child1);

      if (grandParent !== null) {
        // Destroy parent and connect sibling to grandParent.
        if (grandParent.child1 === parent) {
          grandParent.child1 = sibling;
        } else {
          grandParent.child2 = sibling;
        }
        sibling.parent = grandParent;
        this.freeNode(parent);

        // Adjust ancestor bounds.
        let index: TreeNode<T> = grandParent;
        while (index !== null) {
          index = this.balance(index);

          const child1: TreeNode<T> = verify(index.child1);
          const child2: TreeNode<T> = verify(index.child2);

          index.aabb.combine2(child1.aabb, child2.aabb);
          index.height = 1 + Max(child1.height, child2.height);

          index = index.parent;
        }
      } else {
        this.root = sibling;
        sibling.parent = null;
        this.freeNode(parent);
      }

      // this.Validate();
    }

    public balance(A: TreeNode<T>): TreeNode<T> {
      // DEBUG: Assert(A !== null);

      if (A.isLeaf() || A.height < 2) {
        return A;
      }

      const B: TreeNode<T> = verify(A.child1);
      const C: TreeNode<T> = verify(A.child2);

      const balance: number = C.height - B.height;

      // Rotate C up
      if (balance > 1) {
        const F: TreeNode<T> = verify(C.child1);
        const G: TreeNode<T> = verify(C.child2);

        // Swap A and C
        C.child1 = A;
        C.parent = A.parent;
        A.parent = C;

        // A's old parent should point to C
        if (C.parent !== null) {
          if (C.parent.child1 === A) {
            C.parent.child1 = C;
          } else {
            // DEBUG: Assert(C.parent.child2 === A);
            C.parent.child2 = C;
          }
        } else {
          this.root = C;
        }

        // Rotate
        if (F.height > G.height) {
          C.child2 = F;
          A.child2 = G;
          G.parent = A;
          A.aabb.combine2(B.aabb, G.aabb);
          C.aabb.combine2(A.aabb, F.aabb);

          A.height = 1 + Max(B.height, G.height);
          C.height = 1 + Max(A.height, F.height);
        } else {
          C.child2 = G;
          A.child2 = F;
          F.parent = A;
          A.aabb.combine2(B.aabb, F.aabb);
          C.aabb.combine2(A.aabb, G.aabb);

          A.height = 1 + Max(B.height, F.height);
          C.height = 1 + Max(A.height, G.height);
        }

        return C;
      }

      // Rotate B up
      if (balance < -1) {
        const D: TreeNode<T> = verify(B.child1);
        const E: TreeNode<T> = verify(B.child2);

        // Swap A and B
        B.child1 = A;
        B.parent = A.parent;
        A.parent = B;

        // A's old parent should point to B
        if (B.parent !== null) {
          if (B.parent.child1 === A) {
            B.parent.child1 = B;
          } else {
            // DEBUG: Assert(B.parent.child2 === A);
            B.parent.child2 = B;
          }
        } else {
          this.root = B;
        }

        // Rotate
        if (D.height > E.height) {
          B.child2 = D;
          A.child1 = E;
          E.parent = A;
          A.aabb.combine2(C.aabb, E.aabb);
          B.aabb.combine2(A.aabb, D.aabb);

          A.height = 1 + Max(C.height, E.height);
          B.height = 1 + Max(A.height, D.height);
        } else {
          B.child2 = E;
          A.child1 = D;
          D.parent = A;
          A.aabb.combine2(C.aabb, D.aabb);
          B.aabb.combine2(A.aabb, E.aabb);

          A.height = 1 + Max(C.height, D.height);
          B.height = 1 + Max(A.height, E.height);
        }

        return B;
      }

      return A;
    }

    public getHeight(): number {
      if (this.root === null) {
        return 0;
      }

      return this.root.height;
    }

    private static getAreaNode<T>(node: TreeNode<T>): number {
      if (node === null) {
        return 0;
      }

      if (node.isLeaf()) {
        return 0;
      }

      let area: number = node.aabb.getPerimeter();
      area += DynamicTree.getAreaNode(node.child1);
      area += DynamicTree.getAreaNode(node.child2);
      return area;
    }

    public getAreaRatio(): number {
      if (this.root === null) {
        return 0;
      }

      const root: TreeNode<T> = this.root;
      const rootArea: number = root.aabb.getPerimeter();

      const totalArea: number = DynamicTree.getAreaNode(this.root);

      /*
      float32 totalArea = 0.0;
      for (int32 i = 0; i < nodeCapacity; ++i) {
        const TreeNode<T>* node = nodes + i;
        if (node.height < 0) {
          // Free node in pool
          continue;
        }

        totalArea += node.aabb.GetPerimeter();
      }
      */

      return totalArea / rootArea;
    }

    public static computeHeightNode<T>(node: TreeNode<T>): number {
      if (node === null) {
        return 0;
      }

      if (node.isLeaf()) {
        return 0;
      }

      const height1: number = DynamicTree.computeHeightNode(node.child1);
      const height2: number = DynamicTree.computeHeightNode(node.child2);
      return 1 + Max(height1, height2);
    }

    public computeHeight(): number {
      const height: number = DynamicTree.computeHeightNode(this.root);
      return height;
    }

    public validateStructure(node: TreeNode<T>): void {
      if (node === null) {
        return;
      }

      if (node === this.root) {
        // DEBUG: Assert(node.parent === null);
      }

      if (node.isLeaf()) {
        // DEBUG: Assert(node.child1 === null);
        // DEBUG: Assert(node.child2 === null);
        // DEBUG: Assert(node.height === 0);
        return;
      }

      const child1: TreeNode<T> = verify(node.child1);
      const child2: TreeNode<T> = verify(node.child2);

      // DEBUG: Assert(child1.parent === index);
      // DEBUG: Assert(child2.parent === index);

      this.validateStructure(child1);
      this.validateStructure(child2);
    }

    public validateMetrics(node: TreeNode<T>): void {
      if (node === null) {
        return;
      }

      if (node.isLeaf()) {
        // DEBUG: Assert(node.child1 === null);
        // DEBUG: Assert(node.child2 === null);
        // DEBUG: Assert(node.height === 0);
        return;
      }

      const child1: TreeNode<T> = verify(node.child1);
      const child2: TreeNode<T> = verify(node.child2);

      // DEBUG: const height1: number = child1.height;
      // DEBUG: const height2: number = child2.height;
      // DEBUG: const height: number = 1 + Max(height1, height2);
      // DEBUG: Assert(node.height === height);

      const aabb: AABB = DynamicTree.sAabb;
      aabb.combine2(child1.aabb, child2.aabb);

      // DEBUG: Assert(aabb.lowerBound === node.aabb.lowerBound);
      // DEBUG: Assert(aabb.upperBound === node.aabb.upperBound);

      this.validateMetrics(child1);
      this.validateMetrics(child2);
    }

    public validate(): void {
      // DEBUG: this.ValidateStructure(this.root);
      // DEBUG: this.ValidateMetrics(this.root);

      // let freeCount: number = 0;
      // let freeIndex: TreeNode<T> = this.freeList;
      // while (freeIndex !== null) {
      //   freeIndex = freeIndex.parent; // freeIndex = freeIndex.next;
      //   ++freeCount;
      // }

      // DEBUG: Assert(this.GetHeight() === this.ComputeHeight());

      // Assert(this.nodeCount + freeCount === this.nodeCapacity);
    }

    private static getMaxBalanceNode<T>(node: TreeNode<T>, maxBalance: number): number {
      if (node === null) {
        return maxBalance;
      }

      if (node.height <= 1) {
        return maxBalance;
      }

      // DEBUG: Assert(!node.IsLeaf());

      const child1: TreeNode<T> = verify(node.child1);
      const child2: TreeNode<T> = verify(node.child2);
      const balance: number = Abs(child2.height - child1.height);
      return Max(maxBalance, balance);
    }

    public getMaxBalance(): number {
      const maxBalance: number = DynamicTree.getMaxBalanceNode(this.root, 0);

      /*
      int32 maxBalance = 0;
      for (int32 i = 0; i < nodeCapacity; ++i) {
        const TreeNode<T>* node = nodes + i;
        if (node.height <= 1) {
          continue;
        }

        Assert(!node.IsLeaf());

        int32 child1 = node.child1;
        int32 child2 = node.child2;
        int32 balance = Abs(nodes[child2].height - nodes[child1].height);
        maxBalance = Max(maxBalance, balance);
      }
      */

      return maxBalance;
    }

    public rebuildBottomUp(): void {
      /*
      int32* nodes = (int32*)Alloc(nodeCount * sizeof(int32));
      int32 count = 0;

      // Build array of leaves. Free the rest.
      for (int32 i = 0; i < nodeCapacity; ++i) {
        if (nodes[i].height < 0) {
          // free node in pool
          continue;
        }

        if (nodes[i].IsLeaf()) {
          nodes[i].parent = nullNode;
          nodes[count] = i;
          ++count;
        } else {
          FreeNode(i);
        }
      }

      while (count > 1) {
        float32 minCost = maxFloat;
        int32 iMin = -1, jMin = -1;
        for (int32 i = 0; i < count; ++i) {
          AABB aabbi = nodes[nodes[i]].aabb;

          for (int32 j = i + 1; j < count; ++j) {
            AABB aabbj = nodes[nodes[j]].aabb;
            AABB b;
            b.Combine(aabbi, aabbj);
            float32 cost = b.GetPerimeter();
            if (cost < minCost) {
              iMin = i;
              jMin = j;
              minCost = cost;
            }
          }
        }

        int32 index1 = nodes[iMin];
        int32 index2 = nodes[jMin];
        TreeNode<T>* child1 = nodes + index1;
        TreeNode<T>* child2 = nodes + index2;

        int32 parentIndex = AllocateNode();
        TreeNode<T>* parent = nodes + parentIndex;
        parent.child1 = index1;
        parent.child2 = index2;
        parent.height = 1 + Max(child1.height, child2.height);
        parent.aabb.Combine(child1.aabb, child2.aabb);
        parent.parent = nullNode;

        child1.parent = parentIndex;
        child2.parent = parentIndex;

        nodes[jMin] = nodes[count-1];
        nodes[iMin] = parentIndex;
        --count;
      }

      root = nodes[0];
      Free(nodes);
      */

      this.validate();
    }

    private static shiftOriginNode<T>(node: TreeNode<T>, newOrigin: XY): void {
      if (node === null) {
        return;
      }

      if (node.height <= 1) {
        return;
      }

      // DEBUG: Assert(!node.IsLeaf());

      const child1: TreeNode<T> = node.child1;
      const child2: TreeNode<T> = node.child2;
      DynamicTree.shiftOriginNode(child1, newOrigin);
      DynamicTree.shiftOriginNode(child2, newOrigin);

      node.aabb.lowerBound.selfSub(newOrigin);
      node.aabb.upperBound.selfSub(newOrigin);
    }

    public shiftOrigin(newOrigin: XY): void {

      DynamicTree.shiftOriginNode(this.root, newOrigin);

      /*
      // Build array of leaves. Free the rest.
      for (int32 i = 0; i < nodeCapacity; ++i) {
        nodes[i].aabb.lowerBound -= newOrigin;
        nodes[i].aabb.upperBound -= newOrigin;
      }
      */
    }
  }
}
