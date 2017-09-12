"""Using visual Python provides an illustratative animation of a dual
topology alchemical perturbation."""

from numpy import linspace
import vpython as v
from base import atom_colours, atom_sizes, draw_bond_between_atoms, System

ce8_morph = ['o34', 'h38']
ce1_morph = ['c34', 'h38', 'h39', 'h40']


def morph(sys, show_waters=True):
    atoms = {'ce1': [],
             'ce1_morph': [],
             'ce8': [],
             'ce8_morph': []}
    bonds = {'ce1': [],
             'ce1_morph': [],
             'ce8': [],
             'ce8_morph': []}
    for mol in sys.mols:
        if mol.display:
            if mol.name == 't4p':
                if not show_waters:
                    continue
                for at in mol.atoms:
                    if at.type != 'm':
                        a = v.sphere(pos=at.xyz, radius=0.1,
                                     color=atom_colours[at.type], opacity=0.3)
            else:
                for at in mol.atoms:
                    if mol.name == 'ce8' and at.name not in ce8_morph:
                        continue
                    a = v.sphere(pos=at.xyz, radius=atom_sizes[at.type],
                                 color=atom_colours[at.type])
                    if at.name in ce8_morph + ce1_morph:
                        atoms['%s_morph' % mol.name].append(a)
                    else:
                        atoms[mol.name].append(a)

            if mol.name == 't4p':
                o = mol.atoms[0].xyz
                h1 = mol.atoms[1].xyz
                h2 = mol.atoms[2].xyz
                v.cylinder(pos=o, axis=h1 - o, radius=0.02,
                           color=v.color.white, opacity=0.3)
                v.cylinder(pos=o, axis=h2 - o, radius=0.02,
                           color=v.color.white, opacity=0.3)

            if mol.name in ['ce1']:
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'n3 ', 'hn1'))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'n3 ', 'hn2'))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'n3 ', 's1 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 's1 ', 'o1 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 's1 ', 'o2 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 's1 ', 'c8 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c8 ', 'c9 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c9 ', 'c10'))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c8 ', 'c7 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c7 ', 'c6 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c6 ', 'c5 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c10', 'c5 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c9 ', 'h9 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c10', 'h10'))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c7 ', 'h7 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c6 ', 'h6 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c8 ', 'c9 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c5 ', 'n1 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'n1 ', 'n2 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'n2 ', 'c3 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c3 ', 'c1 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c1 ', 'c2 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c2 ', 'n1 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c1 ', 'h1 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c3 ', 'c4 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c1 ', 'c2 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c4 ', 'f1 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c4 ', 'f2 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c4 ', 'f3 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c1 ', 'c2 '))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c2 ', 'c11'))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c11', 'c12'))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c12', 'c13'))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c13', 'c14'))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c14', 'c15'))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c15', 'c16'))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c16', 'c11'))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c12', 'h12'))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c13', 'h13'))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c15', 'h15'))
                bonds[mol.name].append(
                    draw_bond_between_atoms(mol, 'c16', 'h16'))
            if mol.name == 'ce8':
                bonds['ce8_morph'].append(
                    draw_bond_between_atoms(mol, 'c14', 'o34'))
                bonds['ce8_morph'].append(
                    draw_bond_between_atoms(mol, 'o34', 'h38'))
            if mol.name == 'ce1':
                bonds['ce1_morph'].append(
                    draw_bond_between_atoms(mol, 'c14', 'c34'))
                bonds['ce1_morph'].append(
                    draw_bond_between_atoms(mol, 'c34', 'h38'))
                bonds['ce1_morph'].append(
                    draw_bond_between_atoms(mol, 'c34', 'h39'))
                bonds['ce1_morph'].append(
                    draw_bond_between_atoms(mol, 'c34', 'h40'))
    return atoms, bonds


v.scene.background = v.vector(1.0, 1.0, 1.0)
v.scene.width = 800
v.scene.height = 800

limit = 10.0
sys1 = System('dual_topology.pdb', limit)
atoms, bonds = morph(sys1, show_waters=False)

label = v.label(pos=v.vector(3.0, 4.0, 1.0),
                xoffset=0.0,
                yoffset=300.0,
                text="lambda = 0.00",
                line=False,
                height=50,
                color=v.color.red)

v.scene.autoscale = False
v.scene.center = v.vector(4.0, 2.5, 1.0)
v.scene.range = 3

import povexport

while True:
    for lam in linspace(0.0, 1.0, 101):
        v.rate(10)
        # doing lambda**2 looks better for the opacity
        for at in atoms['ce8_morph']:
            at.opacity = lam**2 * 0.7 + 0.3
        for b in bonds['ce8_morph']:
            b.opacity = lam**2 * 0.7 + 0.3
        for at in atoms['ce1_morph']:
            at.opacity = (1 - lam)**2 * 0.7 + 0.3
        for b in bonds['ce1_morph']:
            b.opacity = (1 - lam)**2 * 0.7 + 0.3

        label.text = "lambda = %.2f" % lam
        povexport.export(filename='images/dt%03d.pov' % round(lam*100))
    break
