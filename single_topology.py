"""Using visual Python provides an illustratative animation of a single
topology alchemical perturbation."""

from numpy import linspace
from base import atom_colours, atom_sizes, draw_bond_between_atoms
from base import get_values_for_bond, System
import vpython as v

def draw_scene(sys, show_waters=True):
    """Routine to draw the provided system. Contains hardcoded directions
    for which atoms to bond together. Returns a list of spheres,
    corresponding to atoms, and cylinders corresponding to bonds, for those
    parts of the structure that are transformed."""
    key_atoms = []
    key_bonds = []
    for mol in sys.mols:
        if mol.display:
            if mol.name == 't4p':
                if not show_waters:
                    continue
                for at in mol.atoms:
                    if at.type != 'm':
                        a = v.sphere(pos=at.xyz, radius=atom_sizes[at.type],
                                     color=atom_colours[at.type], opacity=0.3)
            else:
                for at in mol.atoms:
                    a = v.sphere(pos=at.xyz, radius=atom_sizes[at.type],
                                 color=atom_colours[at.type])
                    if at.name in ["c14", "c34", "h38", "h39", "h40"]:
                        key_atoms.append(a)

            if mol.name == 't4p':
                draw_bond_between_atoms(mol, 'o00', 'h01')
                draw_bond_between_atoms(mol, 'o00', 'h02')

            if mol.name in ['ce8', 'ce1']:
                draw_bond_between_atoms(mol, 'n3 ', 'hn1')
                draw_bond_between_atoms(mol, 'n3 ', 'hn2')
                draw_bond_between_atoms(mol, 'n3 ', 's1 ')
                draw_bond_between_atoms(mol, 's1 ', 'o1 ')
                draw_bond_between_atoms(mol, 's1 ', 'o2 ')
                draw_bond_between_atoms(mol, 's1 ', 'c8 ')
                draw_bond_between_atoms(mol, 'c8 ', 'c9 ')
                draw_bond_between_atoms(mol, 'c9 ', 'c10')
                draw_bond_between_atoms(mol, 'c8 ', 'c7 ')
                draw_bond_between_atoms(mol, 'c7 ', 'c6 ')
                draw_bond_between_atoms(mol, 'c6 ', 'c5 ')
                draw_bond_between_atoms(mol, 'c10', 'c5 ')
                draw_bond_between_atoms(mol, 'c9 ', 'h9 ')
                draw_bond_between_atoms(mol, 'c10', 'h10')
                draw_bond_between_atoms(mol, 'c7 ', 'h7 ')
                draw_bond_between_atoms(mol, 'c6 ', 'h6 ')
                draw_bond_between_atoms(mol, 'c8 ', 'c9 ')
                draw_bond_between_atoms(mol, 'c5 ', 'n1 ')
                draw_bond_between_atoms(mol, 'n1 ', 'n2 ')
                draw_bond_between_atoms(mol, 'n2 ', 'c3 ')
                draw_bond_between_atoms(mol, 'c3 ', 'c1 ')
                draw_bond_between_atoms(mol, 'c1 ', 'c2 ')
                draw_bond_between_atoms(mol, 'c2 ', 'n1 ')
                draw_bond_between_atoms(mol, 'c1 ', 'h1 ')
                draw_bond_between_atoms(mol, 'c3 ', 'c4 ')
                draw_bond_between_atoms(mol, 'c1 ', 'c2 ')
                draw_bond_between_atoms(mol, 'c4 ', 'f1 ')
                draw_bond_between_atoms(mol, 'c4 ', 'f2 ')
                draw_bond_between_atoms(mol, 'c4 ', 'f3 ')
                draw_bond_between_atoms(mol, 'c1 ', 'c2 ')
                draw_bond_between_atoms(mol, 'c2 ', 'c11')
                draw_bond_between_atoms(mol, 'c11', 'c12')
                draw_bond_between_atoms(mol, 'c12', 'c13')
                draw_bond_between_atoms(mol, 'c13', 'c14')
                draw_bond_between_atoms(mol, 'c14', 'c15')
                draw_bond_between_atoms(mol, 'c15', 'c16')
                draw_bond_between_atoms(mol, 'c16', 'c11')
                draw_bond_between_atoms(mol, 'c12', 'h12')
                draw_bond_between_atoms(mol, 'c13', 'h13')
                draw_bond_between_atoms(mol, 'c15', 'h15')
                draw_bond_between_atoms(mol, 'c16', 'h16')

            if mol.name == 'ce1':
                key_bonds.append(draw_bond_between_atoms(mol, 'c14', 'c34'))
                key_bonds.append(draw_bond_between_atoms(mol, 'c34', 'h38'))
                key_bonds.append(draw_bond_between_atoms(mol, 'c34', 'h39'))
                key_bonds.append(draw_bond_between_atoms(mol, 'c34', 'h40'))
    return key_atoms, key_bonds


v.scene.background = v.color.white
v.scene.width = 800
v.scene.height = 800

limit = 10.0
sys1 = System('single_topology.pdb', limit)
atoms, bonds = draw_scene(sys1, show_waters=False)

c14, c34, h38, h39, h40 = atoms
c_2_h39 = h39.pos - c34.pos
c_2_h40 = h40.pos - c34.pos
c_2_h38 = h38.pos - c34.pos
c14_2_c34 = c34.pos - c14.pos

c_2_h40_dist = c_2_h40.mag
c14_2_c34_dist = c14_2_c34.mag
unit_c_2_h40 = c_2_h40 / c_2_h40_dist
unit_c14_2_c34 = c14_2_c34 / c14_2_c34_dist

c14_2_c34, c34_2_h38, c34_2_h39, c34_2_h40 = bonds
c2o = 1.40
o2h = 0.95

label = v.label(pos=v.vector(3.00, 4.0, 1.0),
                xoffset=0.0,
                # yoffset=300.0,
                text="lambda = 0.00",
                line=False,
                height=50,
                color=v.color.red)

v.scene.autoscale = False
v.scene.center = v.vector(4.0, 2.5, 1.0)
v.scene.range = 3

import povexport

while True:
    for lam in linspace(0., 1., 101):
        v.rate(10)

        # get correct bond lengths for this lambda value
        c2o_len = (1 - lam) * c14_2_c34_dist + lam * c2o
        o2h_len = (1 - lam) * c_2_h40_dist + lam * o2h

        # update atom positions
        pos_c34 = c14.pos + c2o_len * unit_c14_2_c34
        c34.pos = pos_c34
        h38.pos = pos_c34 + ((1 - lam) * 0.5 + 0.5) * c_2_h38
        h39.pos = pos_c34 + ((1 - lam) * 0.5 + 0.5) * c_2_h39
        h40.pos = pos_c34 + o2h_len * unit_c_2_h40

        # update radius and opacity of dummy atoms
        h38.radius = (1 - lam) * 0.1 + 0.1
        h39.radius = (1 - lam) * 0.1 + 0.1
        h38.opacity = ((1-lam) * 0.5) + 0.5
        h39.opacity = ((1-lam) * 0.5) + 0.5
        c34_2_h38.opacity = ((1-lam) * 0.5) + 0.5
        c34_2_h39.opacity = ((1-lam) * 0.5) + 0.5

        # update bonds
        pos, axis = get_values_for_bond(c14, c34)
        c14_2_c34.axis = axis  # pos_c34 - c14.pos

        pos, axis = get_values_for_bond(c34, h40)
        c34_2_h40.pos = pos
        c34_2_h40.axis = axis

        pos, axis = get_values_for_bond(c34, h39)
        c34_2_h39.pos = pos
        c34_2_h39.axis = axis

        pos, axis = get_values_for_bond(c34, h38)
        c34_2_h38.pos = pos
        c34_2_h38.axis = axis

        # change colour to morph from cyan to red
        c34.color = (1 - lam) * v.color.cyan + lam * v.color.red

        label.text = "lambda = %.2f" % lam
        povexport.export(filename='images/st%03d.pov' % round(lam*100))
    break
