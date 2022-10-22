from ase import Atoms
from element import *
import math

A_site = ['Li', 'Na', 'K', 'Rb', 'Cs']
X_site = ['F', 'Cl', 'Br', 'I']
M_site_mono = ['Pd', 'Li', 'Na', 'K', 'Rb', 'Cs', 'Cu', 'Ag', 'Au', 'Hg', 'In', 'Tl']
M_site_tri = ['Pd', 'Ir', 'Pr', 'Rh', 'Ru', 'La', 'Mo', 'Nd', 'Ni', 'Nb', 'Lu', 'Ce', 'Mn', 'Co', 'Cr', 'Dy', 'Er',
              'Sc', 'Ta', "Tb", 'Eu', 'Y', 'Al', 'Gd', 'Ga', 'In', 'As', 'Sb', 'Bi', 'Fe', "Sb", "Sc", "Sm", "Ti",
              "Tl", "Tm", "V", "Y", 'Au']

M_site_mono_exclusive = [x for x in M_site_mono if x not in M_site_tri]
M_site_tri_exclusive = [x for x in M_site_tri if x not in M_site_mono]
M_site_variable = [x for x in M_site_tri if x in M_site_mono]


def chemical_classifier(atoms: Atoms) -> dict:
    symbols = atoms.get_chemical_symbols()
    atom_dict = {s:symbols.count(s) for s in list(set(symbols))}
    all_elements = list(atom_dict.keys())
    stochiometry = list(sorted([atom_dict[k] for k in all_elements]))
    chemical_dict = {'A_cation': None, 'M_cation_mono': None, 'M_cation_tri': None, 'X_anion': None}

    for e in all_elements:
        if e in X_site:
            chemical_dict['X_anion'] = e
            all_elements.remove(e)

    if (stochiometry == [1, 1, 3]) or (stochiometry == [2, 2, 6]):
        for e in all_elements:
            if (e in M_site_mono) and (e in M_site_tri):
                chemical_dict['M_cation_mono'] = e
                chemical_dict['M_cation_tri'] = e
                all_elements.remove(e)
        chemical_dict['A_cation'] = all_elements[-1]
        assert (None not in [chemical_dict[k] for k in chemical_dict.keys()])
    elif stochiometry == [1, 3, 6]:
        for e in all_elements:
            if (e in A_site) and (e in M_site_mono):
                chemical_dict['A_cation'] = e
                chemical_dict['M_cation_mono'] = e
                all_elements.remove(e)
        chemical_dict['M_cation_tri'] = all_elements[-1]
        assert (None not in [chemical_dict[k] for k in chemical_dict.keys()])
    elif stochiometry == [1, 1, 2, 6]:
        for e in all_elements:
            if (e in A_site) and (atom_dict[e] == 2):
                chemical_dict['A_cation'] = e
                all_elements.remove(e)

        M_site_elements = all_elements

        for e in M_site_elements:
            if e in M_site_mono_exclusive:
                chemical_dict['M_cation_mono'] = e
                M_site_elements.remove(e)
                if len(M_site_elements) == 1:
                    chemical_dict['M_cation_tri'] = M_site_elements[-1]
            elif e in M_site_tri_exclusive:
                chemical_dict['M_cation_tri'] = e
                M_site_elements.remove(e)
                if len(M_site_elements) == 1:
                    chemical_dict['M_cation_mono'] = M_site_elements[-1]
        if len(M_site_elements) == 2:
            if all([m in M_site_variable for m in M_site_elements]):
                # cannot really tell which one in which charge state, randomly assign one
                print('variable valence, randomly assigned')
                if 'Pd' not in M_site_elements:
                    chemical_dict['M_cation_mono'] = M_site_elements[0]
                    chemical_dict['M_cation_tri'] = M_site_elements[1]
                else:
                    chemical_dict['M_cation_tri'] = 'Pd'
                    M_site_elements.remove('Pd')
                    chemical_dict['M_cation_mono'] = M_site_elements[-1]
                M_site_elements = []
        assert (None not in [chemical_dict[k] for k in chemical_dict.keys()])

    return chemical_dict

def geometric_fingerprint(atoms: Atoms):
    chemistry = chemical_classifier(atoms)
    r_a = shannon_radii[chemistry['A_cation']]["1"]["VI"]['r_ionic']
    r_m = shannon_radii[chemistry['M_cation_mono']]["1"]["VI"]['r_ionic']
    r_mp = shannon_radii[chemistry['M_cation_tri']]["3"]["VI"]['r_ionic']
    r_x = shannon_radii[chemistry['X_anion']]["-1"]["VI"]['r_ionic']
    return chemistry, octahedral_factor(r_m, r_mp, r_x), octahedral_mismatch(r_m, r_mp,
                                                                             r_x), generalised_tolerance_factor(r_a,
                                                                                                                r_m,
                                                                                                                r_mp,
                                                                                                                r_x)

def octahedral_factor(r_m, r_mp, r_x):
    return (r_m + r_mp) / (2.0 * r_x)


def octahedral_mismatch(r_m, r_mp, r_x):
    return abs(r_m - r_mp) / (2.0 * r_x)


def generalised_tolerance_factor(r_a, r_m, r_mp, r_x):
    nominator = r_a + r_x
    denominator = (r_m + r_mp) / 2.0 + r_x
    denominator = denominator ** 2 + (r_m - r_mp) ** 2 / 4
    denominator = math.sqrt(denominator) * math.sqrt(2)
    return nominator / denominator

