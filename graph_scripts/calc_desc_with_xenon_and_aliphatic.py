"""Calculate modified synthon descriptors that include xenon"""
import csv
import gzip
import os
import tempfile
from argparse import ArgumentParser
from typing import List, Tuple

from joblib import Parallel, delayed
from oddt import shape, toolkit
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures, rdMolDescriptors, Mol, rdmolops
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from tqdm import tqdm


def num_rings(mol):
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    return num_rings


def num_atoms(mol):
    return mol.GetNumAtoms()


def is_in_same_ring(idx1: int, idx2: int, bond_rings) -> bool:
    """
    Checks if two atoms (using atom indices) are in the same ring
    """
    for bond_ring in bond_rings:
        if idx1 in bond_ring and idx2 in bond_ring:
            return True
    return False


def get_fragments(mol: Mol, bonds: list) -> Tuple[list, list]:
    """
    Cleaves the molecule into fragments at specified bonds.
    """
    _fragments = Chem.FragmentOnBonds(mol, bonds)
    fragments = Chem.GetMolFrags(_fragments, asMols=True)
    nonring_fragments = []

    for fragment in fragments:
        # check the fragment does not contain rings
        numRings = rdMolDescriptors.CalcNumRings(fragment)
        if numRings == 0:
            nonring_fragments.append(fragment)

    return nonring_fragments


def max_path(mol, synthon=True) -> int:
        """
        Get the maximum path length within the substructure.
        """
        if synthon:
            try:
                xe_atom = mol.GetSubstructMatch(Chem.MolFromSmiles('[Xe]'))[0]
            except:
                xe_atom = None
        maxPath = 0
        numAtoms = mol.GetNumHeavyAtoms()
        for i in range(numAtoms):  # look for paths up the length of number of atoms
            # get max consecutive path length in linker/sidechain
            if synthon and xe_atom:
                paths = rdmolops.FindAllPathsOfLengthN(mol, i, useBonds=True, rootedAtAtom=xe_atom)
            else:
                paths = rdmolops.FindAllPathsOfLengthN(mol, i, useBonds=True)
            if len(paths) > 0:
                maxPath = i

        return maxPath


def get_bonds_to_cleave(mol: Mol):
    """
    Get the bonds that are attached to rings (cleaving at these bonds will give the linker and sidechain
    substructures.
    """
    cleave_bonds = []

    # record original indices of atoms and bonds as properties
    for atom in mol.GetAtoms():
        atom.SetIntProp("orig_idx", atom.GetIdx())
    for bond in mol.GetBonds():
        bond.SetIntProp("orig_idx", bond.GetIdx())

    # get set of bonds that are in rings
    ring_info = mol.GetRingInfo()
    bond_rings = (
        ring_info.BondRings()
    )  # gives tuples with bonds that are in the rings
    ring_bonds = set()
    for ring_bond_idxs in bond_rings:  # get set
        for idx in ring_bond_idxs:
            ring_bonds.add(idx)

    # get the non-ring bonds
    all_bonds_idx = [bond.GetIdx() for bond in mol.GetBonds()]
    non_ring_bonds = set(all_bonds_idx) - ring_bonds

    for bond_idx in non_ring_bonds:
        # get the indices of the atoms on each end of the bond
        bgn_idx = mol.GetBondWithIdx(bond_idx).GetBeginAtomIdx()
        end_idx = mol.GetBondWithIdx(bond_idx).GetEndAtomIdx()
        if (
            mol.GetBondWithIdx(bond_idx).GetBondTypeAsDouble() == 1.0
        ):  # if a single bond
            # GetBondTypeAsDouble() - 1.0 for single, 1.5 for aromatic, 2.0 for double
            # if one of the atoms on the side of the bond is in the ring, append to cleave_bonds
            if (
                mol.GetAtomWithIdx(bgn_idx).IsInRing()
                + mol.GetAtomWithIdx(end_idx).IsInRing()
                == 1
            ):
                bond = mol.GetBondWithIdx(bond_idx)
                orig_idx = bond.GetIntProp("orig_idx")
                cleave_bonds.append(orig_idx)
            # if the bond atoms are in two different rings, append to cleave_bonds
            elif (
                not is_in_same_ring(bgn_idx, end_idx, bond_rings)
                and mol.GetAtomWithIdx(bgn_idx).IsInRing()
                + mol.GetAtomWithIdx(end_idx).IsInRing()
                == 2
            ):
                bond = mol.GetBondWithIdx(bond_idx)
                orig_idx = bond.GetIntProp("orig_idx")
                cleave_bonds.append(orig_idx)
    return cleave_bonds


def max_path_nonring(mol):
    bonds = self.get_bonds_to_cleave(mol)

    maxpath = 0
    if len(bonds) == 0:
        return 0
    else:
        nonring_fragments = self.get_fragments(mol, bonds)
        if len(nonring_fragments) > 0:
            max_paths = [max_path(f) for f in nonring_fragments]
            maxpath = max(max_paths)

    return maxpath

def calculate_usrcat(mol, idx, tmp_dir, xenon):
    tmp_file = os.path.join(tmp_dir, f"usrcatxe_{idx}.mol")
    try:
        Chem.MolToMolFile(mol, tmp_file)
        oddt_mol = next(toolkit.readfile('mol', tmp_file))
        oddt_mol.addh()
        usrcat = list(shape.usr_cat(oddt_mol, xenon))
        os.remove(tmp_file)
        return ";".join(map(str, usrcat))
    except Exception as e:
        print(e)
        return None


def calc_pharmfp(mol, fdefName):
    featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    sigFactory = SigFactory(featFactory, maxPointCount=2)
    sigFactory.SetBins([(0, 2), (2, 5), (5, 8)])
    sigFactory.Init()
    sigFactory.GetSigSize()
    fp = Generate.Gen2DFingerprint(mol, sigFactory)
    fp = list(fp)
    return ";".join(map(str, fp))


def calculate_properties(smiles, idx, tmp_dir, fdefName):
    mol = Chem.MolFromSmiles(smiles)
    rings = num_rings(mol)
    atoms = num_atoms(mol)
    maxpathlength = max_path(mol, True)
    usrcat = calculate_usrcat(mol, idx, tmp_dir, xenon=True)
    pharm_fp = calc_pharmfp(mol, fdefName)
    return smiles, str(rings), str(atoms), str(maxpathlength), usrcat, pharm_fp


def main():
    parser = ArgumentParser()
    parser.add_argument('--synthon_gzip_file')
    parser.add_argument('--fdef_file')
    parser.add_argument('--n_cpus', type=int)
    parser.add_argument('--output_dir')
    args = parser.parse_args()
    synthons = []

    print('reading synthons')

    with gzip.GzipFile(filename=args.synthon_gzip_file) as file:
        for line in tqdm(file):
            l = line.decode().strip()
            synthons.append(l)

    print('synthons read')
    print('calculating descriptors')

    idxs = [i for i in range(len(synthons))]
    tmpdir = tempfile.TemporaryDirectory(prefix='/tmp/usrcat_')
    properties = Parallel(n_jobs=args.n_cpus, backend="multiprocessing")(
        delayed(calculate_properties)(smiles, idx, tmpdir.name, args.fdef_file)
        for smiles, idx in tqdm(zip(synthons, idxs))
    )

    print('properties calculated')
    print('writing descriptor file')

    desc_fname = os.path.join(args.output_dir, f"synthons_xenon_aliphatic_ring_descriptors.csv.gz")
    desc_f = gzip.open(desc_fname, 'wt')
    desc_writer = csv.writer(desc_f, delimiter=',')

    for prop in tqdm(properties):
        desc_writer.writerow(prop)

    desc_f.close()
    print('descriptor file written')


if __name__ == "__main__":
    main()
