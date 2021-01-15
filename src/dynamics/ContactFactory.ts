namespace b2 {

  export class ContactRegister {
    public pool: Contact[] = [];
    public createFcn: (() => Contact) | null = null;
    public destroyFcn: ((contact: Contact) => void) | null = null;
    public primary: boolean = false;
  }

  export class ContactFactory {
    public readonly registers: ContactRegister[][] = [];

    constructor() {
      this.InitializeRegisters();
    }

    private AddType(createFcn: () => Contact, destroyFcn: (contact: Contact) => void, typeA: ShapeType, typeB: ShapeType): void {
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

    private InitializeRegisters(): void {
      for (let i: number = 0; i < ShapeType.e_shapeTypeCount; i++) {
        this.registers[i] = [];
        for (let j: number = 0; j < ShapeType.e_shapeTypeCount; j++) {
          this.registers[i][j] = new ContactRegister();
        }
      }

      this.AddType(          CircleContact.Create,           CircleContact.Destroy, ShapeType.e_circleShape,  ShapeType.e_circleShape);
      this.AddType(PolygonAndCircleContact.Create, PolygonAndCircleContact.Destroy, ShapeType.e_polygonShape, ShapeType.e_circleShape);
      this.AddType(         PolygonContact.Create,          PolygonContact.Destroy, ShapeType.e_polygonShape, ShapeType.e_polygonShape);
      this.AddType(   EdgeAndCircleContact.Create,    EdgeAndCircleContact.Destroy, ShapeType.e_edgeShape,    ShapeType.e_circleShape);
      this.AddType(  EdgeAndPolygonContact.Create,   EdgeAndPolygonContact.Destroy, ShapeType.e_edgeShape,    ShapeType.e_polygonShape);
      this.AddType(  ChainAndCircleContact.Create,   ChainAndCircleContact.Destroy, ShapeType.e_chainShape,   ShapeType.e_circleShape);
      this.AddType( ChainAndPolygonContact.Create,  ChainAndPolygonContact.Destroy, ShapeType.e_chainShape,   ShapeType.e_polygonShape);
    }

    public Create(fixtureA: Fixture, indexA: number, fixtureB: Fixture, indexB: number): Contact | null {
      const typeA: ShapeType = fixtureA.GetType();
      const typeB: ShapeType = fixtureB.GetType();

      // DEBUG: Assert(0 <= typeA && typeA < ShapeType.e_shapeTypeCount);
      // DEBUG: Assert(0 <= typeB && typeB < ShapeType.e_shapeTypeCount);

      const reg: ContactRegister = this.registers[typeA][typeB];
      if (reg.createFcn) {
        const c: Contact = reg.createFcn();
        if (reg.primary) {
          c.Reset(fixtureA, indexA, fixtureB, indexB);
        } else {
          c.Reset(fixtureB, indexB, fixtureA, indexA);
        }
        return c;
      } else {
        return null;
      }
    }

    public Destroy(contact: Contact): void {
      const typeA: ShapeType = contact.fixtureA.GetType();
      const typeB: ShapeType = contact.fixtureB.GetType();

      // DEBUG: Assert(0 <= typeA && typeB < ShapeType.e_shapeTypeCount);
      // DEBUG: Assert(0 <= typeA && typeB < ShapeType.e_shapeTypeCount);

      const reg: ContactRegister = this.registers[typeA][typeB];
      if (reg.destroyFcn) {
        reg.destroyFcn(contact);
      }
    }
  }

}
