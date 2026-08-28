import argparse
import os
import pickle
import shutil

import torch


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--val-ratio", type=float, default=0.2)
    args = parser.parse_args()

    output = os.path.abspath(args.output)
    os.makedirs(os.path.join(output, "ligands"), exist_ok=True)
    entries = []
    source_reports = []
    for source in args.inputs:
        source = os.path.abspath(source)
        with open(os.path.join(source, "index.pkl"), "rb") as handle:
            source_entries = pickle.load(handle)
        prefix = os.path.basename(source)
        source_reports.append((prefix, len(source_entries)))
        source_pocket = os.path.join(source, "pocket.pdb")
        target_pocket = "%s_pocket.pdb" % prefix
        shutil.copyfile(source_pocket, os.path.join(output, target_pocket))
        for index, (_, ligand_name, _, rmsd) in enumerate(source_entries):
            target_ligand = os.path.join("ligands", "%s_ligand_%04d.sdf" % (prefix, index))
            shutil.copyfile(os.path.join(source, ligand_name), os.path.join(output, target_ligand))
            entries.append((target_pocket, target_ligand, None, rmsd))

    if len(entries) < 10:
        raise RuntimeError("merged dataset is too small")
    split_at = max(1, int(len(entries) * (1.0 - args.val_ratio)))
    if split_at >= len(entries):
        split_at = len(entries) - 1
    train_entries = entries[:split_at]
    val_entries = entries[split_at:]
    with open(os.path.join(output, "index.pkl"), "wb") as handle:
        pickle.dump(entries, handle)
    torch.save({
        "train": [(item[0], item[1]) for item in train_entries],
        "test": [(item[0], item[1]) for item in val_entries],
    }, os.path.join(output, "split_by_name.pt"))
    with open(os.path.join(output, "build_report.txt"), "w", encoding="utf-8") as handle:
        handle.write("sources=%s\n" % source_reports)
        handle.write("total=%d\ntrain=%d\nval=%d\n" % (len(entries), len(train_entries), len(val_entries)))
    print("merged=%d train=%d val=%d output=%s" % (len(entries), len(train_entries), len(val_entries), output))


if __name__ == "__main__":
    main()
