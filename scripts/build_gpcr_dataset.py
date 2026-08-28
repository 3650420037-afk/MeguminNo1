import argparse
import csv
import os
import pickle
import random
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from utils.protein_ligand import PDBProtein


def read_smiles(csv_path):
    with open(csv_path, encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle)
        field = "canonical_smiles"
        if field not in (rows.fieldnames or []):
            raise ValueError("CSV must contain canonical_smiles")
        return [row[field].strip() for row in rows if row.get(field, "").strip()]


def build_conformer(smiles, center, seed):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("invalid SMILES")
    supported_elements = {6, 7, 8, 9, 15, 16, 17}
    unsupported = sorted({atom.GetAtomicNum() for atom in molecule.GetAtoms()} - supported_elements)
    if unsupported:
        raise ValueError("unsupported atomic numbers: %s" % unsupported)
    molecule = Chem.AddHs(molecule)
    status = AllChem.EmbedMolecule(molecule, randomSeed=seed, useRandomCoords=True)
    if status != 0:
        raise ValueError("3D embedding failed")
    try:
        AllChem.MMFFOptimizeMolecule(molecule, maxIters=200)
    except Exception:
        AllChem.UFFOptimizeMolecule(molecule, maxIters=200)
    conformer = molecule.GetConformer()
    coordinates = np.asarray([list(conformer.GetAtomPosition(i)) for i in range(molecule.GetNumAtoms())])
    heavy_mask = np.asarray([atom.GetAtomicNum() != 1 for atom in molecule.GetAtoms()])
    translation = np.asarray(center, dtype=np.float32) - coordinates[heavy_mask].mean(axis=0)
    for index in range(molecule.GetNumAtoms()):
        position = conformer.GetAtomPosition(index)
        shifted = np.asarray([position.x, position.y, position.z]) + translation
        conformer.SetAtomPosition(index, shifted.tolist())
    molecule = Chem.RemoveHs(molecule)
    if molecule.GetNumConformers() != 1:
        raise ValueError("missing final conformer")
    if Descriptors.HeavyAtomCount(molecule) < 2:
        raise ValueError("too few heavy atoms")
    return molecule


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--center", required=True, type=lambda value: [float(x) for x in value.split(",")])
    parser.add_argument("--output", required=True)
    parser.add_argument("--pocket-radius", type=float, default=12.0)
    parser.add_argument("--val-ratio", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=2021)
    args = parser.parse_args()

    if len(args.center) != 3:
        raise ValueError("center must contain three coordinates")
    random.seed(args.seed)
    output = os.path.abspath(args.output)
    os.makedirs(os.path.join(output, "ligands"), exist_ok=True)

    protein = PDBProtein(args.pdb)
    residues = protein.query_residues_radius(args.center, args.pocket_radius, criterion="center_of_mass")
    if not residues:
        raise ValueError("no protein residues found near the pocket center")
    pocket_name = "pocket.pdb"
    with open(os.path.join(output, pocket_name), "w", encoding="ascii") as handle:
        handle.write(protein.residues_to_pdb_block(residues, name="GPCR_POCKET"))

    smiles_values = read_smiles(args.csv)
    writer = Chem.SDWriter(os.path.join(output, "ligands", "placeholder.sdf"))
    writer.close()
    os.remove(os.path.join(output, "ligands", "placeholder.sdf"))
    entries = []
    skipped = []
    for index, smiles in enumerate(smiles_values):
        try:
            molecule = build_conformer(smiles, args.center, args.seed + index)
            ligand_name = os.path.join("ligands", "ligand_%04d.sdf" % index)
            writer = Chem.SDWriter(os.path.join(output, ligand_name))
            writer.write(molecule)
            writer.close()
            entries.append((pocket_name, ligand_name, None, "0.0"))
        except Exception as error:
            skipped.append((index, str(error)))

    if len(entries) < 5:
        raise RuntimeError("fewer than five valid ligand pairs were generated")
    random.shuffle(entries)
    split_at = max(1, int(len(entries) * (1.0 - args.val_ratio)))
    train_entries = entries[:split_at]
    val_entries = entries[split_at:]
    if not val_entries:
        val_entries = train_entries[-1:]
        train_entries = train_entries[:-1]
    index_path = os.path.join(output, "index.pkl")
    with open(index_path, "wb") as handle:
        pickle.dump(entries, handle)
    split = {
        "train": [(item[0], item[1]) for item in train_entries],
        "test": [(item[0], item[1]) for item in val_entries],
    }
    torch.save(split, os.path.join(output, "split_by_name.pt"))
    with open(os.path.join(output, "build_report.txt"), "w", encoding="utf-8") as handle:
        handle.write("source_csv=%s\nsource_pdb=%s\n" % (args.csv, args.pdb))
        handle.write("total_smiles=%d\nvalid_pairs=%d\nskipped=%d\n" % (len(smiles_values), len(entries), len(skipped)))
        handle.write("train=%d\nval=%d\n" % (len(train_entries), len(val_entries)))
        for index, reason in skipped:
            handle.write("skipped_%d=%s\n" % (index, reason))
    print("valid_pairs=%d train=%d val=%d skipped=%d output=%s" % (
        len(entries), len(train_entries), len(val_entries), len(skipped), output
    ))


if __name__ == "__main__":
    main()
