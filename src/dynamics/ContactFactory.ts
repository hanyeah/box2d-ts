namespace b2 {

  export class ContactRegister {
    public pool: Contact[] = [];
    public createFcn: (() => Contact) = null;
    public destroyFcn: ((contact: Contact) => void) = null;
    public primary: boolean = false;
  }

  export class ContactFactory {
    public readonly registers: ContactRegister[][] = [];

    constructor() {
      this.initializeRegisters();
    }

    private addType(createFcn: () => Contact, destroyFcn: (contact: Contact) => void, typeA: ShapeType, typeB: ShapeType): void {
      const pool: Contact[] = [];

      function poolCreateFcn(): Contact {
        return pool.pop() || createFcn();
      }

      function poolDestroyFcn(contact: Contact): void {
        pool.push(contact);
      }

      this.registers[typeA][typeB].pool = pool;
      this.registers[typeA][typeB].createFcn = poolCreateFcn; // createFcn;
      this.registers[typeA][typeB].destroyFcn = poolDestroyFcn; // destroyFcn;
      this.registers[typeA][typeB].primary = true;

      if (typeA !== typeB) {
        this.registers[typeB][typeA].pool = pool;
        this.registers[typeB][typeA].createFcn = poolCreateFcn; // createFcn;
        this.registers[typeB][typeA].destroyFcn = poolDestroyFcn; // destroyFcn;
        this.registers[typeB][typeA].primary = false;
      }
    }

    private initializeRegisters(): void {
      for (let i: number = 0; i < ShapeType.ShapeTypeCount; i++) {
        this.registers[i] = [];
        for (let j: number = 0; j < ShapeType.ShapeTypeCount; j++) {
          this.registers[i][j] = new ContactRegister();
        }
      }

      this.addType(          CircleContact.create,           CircleContact.destroy, ShapeType.CircleShape,  ShapeType.CircleShape);
      this.addType(PolygonAndCircleContact.create, PolygonAndCircleContact.destroy, ShapeType.PolygonShape, ShapeType.CircleShape);
      this.addType(         PolygonContact.create,          PolygonContact.destroy, ShapeType.PolygonShape, ShapeType.PolygonShape);
      this.addType(   EdgeAndCircleContact.create,    EdgeAndCircleContact.destroy, ShapeType.EdgeShape,    ShapeType.CircleShape);
      this.addType(  EdgeAndPolygonContact.create,   EdgeAndPolygonContact.destroy, ShapeType.EdgeShape,    ShapeType.PolygonShape);
      this.addType(  ChainAndCircleContact.create,   ChainAndCircleContact.destroy, ShapeType.ChainShape,   ShapeType.CircleShape);
      this.addType( ChainAndPolygonContact.create,  ChainAndPolygonContact.destroy, ShapeType.ChainShape,   ShapeType.PolygonShape);
    }

    public create(fixtureA: Fixture, indexA: number, fixtureB: Fixture, indexB: number): Contact {
      const typeA: ShapeType = fixtureA.getType();
      const typeB: ShapeType = fixtureB.getType();

      // DEBUG: Assert(0 <= typeA && typeA < ShapeType.e_shapeTypeCount);
      // DEBUG: Assert(0 <= typeB && typeB < ShapeType.e_shapeTypeCount);

      const reg: ContactRegister = this.registers[typeA][typeB];
      if (reg.createFcn) {
        const c: Contact = reg.createFcn();
        if (reg.primary) {
          c.reset(fixtureA, indexA, fixtureB, indexB);
        } else {
          c.reset(fixtureB, indexB, fixtureA, indexA);
        }
        return c;
      } else {
        return null;
      }
    }

    public destroy(contact: Contact): void {
      const typeA: ShapeType = contact.fixtureA.getType();
      const typeB: ShapeType = contact.fixtureB.getType();

      // DEBUG: Assert(0 <= typeA && typeB < ShapeType.e_shapeTypeCount);
      // DEBUG: Assert(0 <= typeA && typeB < ShapeType.e_shapeTypeCount);

      const reg: ContactRegister = this.registers[typeA][typeB];
      if (reg.destroyFcn) {
        reg.destroyFcn(contact);
      }
    }
  }

}
