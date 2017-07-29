"""Base module providing common routines and data for both
the single and dual alchemical approaches."""
import vpython as v

atom_colours = {'h': v.color.white,
                'n': v.color.blue,
                'c': v.color.cyan,
                'o': v.color.red,
                's': v.color.yellow,
                'f': v.color.green,
                'd': v.color.magenta}

atom_sizes = {'h': 0.2,
              'n': 0.3,
              'c': 0.3,
              'o': 0.3,
              's': 0.3,
              'f': 0.3,
              'd': 0.3}


class Atom(object):
    """Atom object. Stores the atom name, type (based on name),
    position and molecule name. Takes a pdb line as input."""
    def __init__(self, line):
        """line: String, pdb file line containing an ATOM record"""

        self.name = line[13:16].lower()
        self.type = self.name[0]

        self.xyz = v.vector(
            *map(float, [line[31:38], line[39:46], line[47:54]]))
        self.mol_name = line[17:20].lower()


class Molecule(object):
    """Molecule object. Stores constituent atoms."""
    def __init__(self):
        self.atoms = []

    def add_atom(self, line):
        """Add a new atom to the molecule, based on a pdb line."""
        self.atoms.append(Atom(line))
        self.name = self.atoms[-1].mol_name

    def get_atom_named(self, name):
        """Find and return the atom within this molecule with
        the provided name."""
        for i in self.atoms:
            if i.name == name:
                return i

    def test_display(self, lim):
        """Check whether or not to display the molecule based its,
        position and set self.display as True or False correspondingly.
        True if all atoms have x, y and z coordinates -lim < X < lim,
        othersiwe False.
        """
        if False not in [in_range(i.xyz, lim) for i in self.atoms]:
            self.display = True
        else:
            self.display = False


class System(object):
    """Container for a molecular system. Reads a pdb file at
    file_loc and constructs Molecule objects for each residue
    as determined by TER cards. Determine whether loaded
    molecules should be displayed based on their position and
    the passed value of lim. See Molecule.test_display for
    details."""
    def __init__(self, file_loc, lim):
        """file_loc : String providing pdb file location
        lim: All molecules with atoms positioned outside of the range of
             [-lim,lim] in any dimension are not included in the System.
        """
        self.mols = []

        mol = Molecule()
        with open(file_loc) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATOM'):
                    mol.add_atom(line)

                if line.startswith('TER'):
                    self.mols.append(mol)
                    mol = Molecule()

        if len(mol.atoms) != 0:
            self.mols.append(mol)

        # recenter system coordinates on the origin
        center = v.vector(0., 0., 0,)
        count = 0
        for mol in self.mols:
            for at in mol.atoms:
                center += at.xyz
                count += 1

        for mol in self.mols:
            for at in mol.atoms:
                at.xyz -= center / count

        # set display values for all molecules
        for mol in self.mols:
            mol.test_display(lim)


def in_range(vals, lim):
    """Test the provided vector vals to ensure that all constituent
    values lie in the range [-lim,lim]."""
    for i in (vals.x, vals.y, vals.z):
        if i < -lim or i > lim:
            return False
    return True


def draw_bond_between_atoms(mol, name1, name2):
    """Draw a bond between the two atoms with names name1 and name2
    belonging to mol."""

    at1 = mol.get_atom_named(name1)
    at2 = mol.get_atom_named(name2)

    pos, axis = get_values_for_bond(at1, at2)
    return v.cylinder(pos=pos,
                      axis=axis,
                      radius=0.04,
                      color=v.color.white)


def get_values_for_bond(at1, at2):
    """Based on two input atoms (or spheres) provide the position
    and axis vectors for the bond between them."""

    # the below provides plasticity for both Atom objects and
    # spheres
    try:
        pos1, pos2 = at1.xyz, at2.xyz
        r1, r2 = atom_sizes[at1.type], atom_sizes[at2.type]
    except AttributeError:
        pos1, pos2 = at1.pos, at2.pos
        r1, r2 = at1.radius, at2.radius

    diff = pos1 - pos2

    # we want to draw the bond between the outer points of the
    # atom spheres, this stops things looking bad with transparency
    pos = pos1 - diff/diff.mag*r1
    axis = -diff/diff.mag*(diff.mag-r1-r2)
    return pos, axis
