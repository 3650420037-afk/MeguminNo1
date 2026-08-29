# -*- coding: utf-8 -*-
"""Compute pocket centers from co-crystal ligands in RCSB mmCIF files.

Parses the _atom_site loop of a PDBx/mmCIF file, groups HETATM records by
(auth chain, comp id), computes the unweighted geometric centroid of the
heavy atoms of the target ligand copy, and validates the center by counting
receptor residues within 12 A using the same PDBProtein helper as
build_gpcr_dataset.py.

Usage:
    python calc_pocket_centers.py
"""
import os
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from utils.protein_ligand import PDBProtein

# target: (cif_path, receptor_pdb_path, ligand_comp_id)
TARGETS = {
    # 3PBL co-crystal ligand: eticlopride, CCD code ETQ
    "D3": (
        r"D:\MMModel\靶点结构\3PBL.cif",
        r"D:\MMModel\靶点结构\3PBL_D3多巴胺受体.pdb",
        "ETQ",
    ),
    # 4IB4 co-crystal ligand: ergotamine (C33 H35 N5 O5), CCD code ERM
    "5HT2B": (
        r"D:\MMModel\靶点结构\4IB4.cif",
        r"D:\MMModel\靶点结构\4IB4_5HT2B受体.pdb",
        "ERM",
    ),
}

SKIP_ALTLOC = {"", ".", "?", "A"}


def parse_atom_site_loop(cif_path):
    """Yield dicts for every atom_site record of an mmCIF file."""
    records = []
    with open(cif_path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    index = 0
    while index < len(lines):
        if lines[index].strip().startswith("loop_"):
            block_start = index + 1
            names = []
            cursor = block_start
            while cursor < len(lines):
                token = lines[cursor].strip()
                if token.startswith("_atom_site."):
                    names.append(token)
                    cursor += 1
                elif token == "" or token.startswith("#"):
                    cursor += 1
                else:
                    break
            if names:
                columns = [name.split(".", 1)[1] for name in names]
                while cursor < len(lines):
                    line = lines[cursor].strip()
                    if line == "" or line == "#" or line == "loop_" or line.startswith("_"):
                        break
                    parts = line.split()
                    if len(parts) != len(columns):
                        break
                    records.append(dict(zip(columns, parts)))
                    cursor += 1
                index = cursor
                continue
        index += 1
    return records


def pdb_chain_ids(pdb_path):
    chains = set()
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                chains.add(line[21:22].strip())
    return chains


def heavy_atoms(records, comp_id):
    """Group heavy atoms of comp_id by auth chain. Returns {chain: [(x,y,z)]}."""
    groups = {}
    for row in records:
        if row.get("group_PDB") != "HETATM":
            continue
        if row.get("comp_id", row.get("label_comp_id")) != comp_id:
            continue
        if row.get("pdbx_PDB_model_num", "1") != "1":
            continue
        if row.get("alt_id", row.get("label_alt_id", ".")) not in SKIP_ALTLOC:
            continue
        element = row.get("type_symbol", row.get("label_atom_id", "?")).upper()
        if element in ("H", "D"):
            continue
        chain = row.get("auth_asym_id", row.get("label_asym_id", "?"))
        groups.setdefault(chain, []).append(
            (float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"]))
        )
    return groups


def centroid(points):
    count = len(points)
    return (
        sum(p[0] for p in points) / count,
        sum(p[1] for p in points) / count,
        sum(p[2] for p in points) / count,
    )


def main():
    for target, (cif_path, pdb_path, comp_id) in TARGETS.items():
        records = parse_atom_site_loop(cif_path)
        het_comps = {}
        for row in records:
            if row.get("group_PDB") == "HETATM" and row.get("pdbx_PDB_model_num", "1") == "1":
                name = row.get("comp_id", row.get("label_comp_id"))
                het_comps[name] = het_comps.get(name, 0) + 1
        print("=== %s ===" % target)
        print("cif=%s" % cif_path)
        print("hetatm_comp_ids=%s" % sorted(het_comps.items()))
        chains = pdb_chain_ids(pdb_path)
        print("receptor_pdb_chains=%s" % sorted(chains))
        groups = heavy_atoms(records, comp_id)
        if not groups:
            raise SystemExit("ligand %s not found in %s" % (comp_id, cif_path))
        selected = None
        for chain in sorted(groups):
            points = groups[chain]
            center = centroid(points)
            print(
                "copy chain=%s heavy_atoms=%d centroid=(%.6f, %.6f, %.6f)"
                % (chain, len(points), center[0], center[1], center[2])
            )
            if selected is None and chain in chains:
                selected = (chain, center)
        if selected is None:
            selected = (sorted(groups)[0], centroid(groups[sorted(groups)[0]]))
        chain, center = selected
        protein = PDBProtein(pdb_path)
        with open(os.devnull, "w") as devnull:
            original_stdout = sys.stdout
            sys.stdout = devnull
            try:
                residues = protein.query_residues_radius(center, 12.0, criterion="center_of_mass")
            finally:
                sys.stdout = original_stdout
        print(
            "SELECTED %s chain=%s center=\"%.3f,%.3f,%.3f\" residues_within_12A=%d"
            % (comp_id, chain, center[0], center[1], center[2], len(residues))
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
