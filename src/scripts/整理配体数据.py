# -*- coding: utf-8 -*-
"""整理下载的 ChEMBL 配体数据为 CSV,并生成训练文件清单"""
import json, csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import KNOWN_DRUGS, LIGAND_DATA, LIGANDS, src_on_path
src_on_path()

def chembl_activities_to_smiles(json_path, target_name, out_csv):
    """从 ChEMBL activity JSON 中提取分子 SMILES, 去重后保存 CSV"""
    with open(json_path, encoding='utf-8') as f:
        data = json.load(f)
    acts = data.get('activities', [])
    seen = set()
    rows = []
    for a in acts:
        mol_id = a.get('molecule_chembl_id')
        smiles = a.get('canonical_smiles')
        if not mol_id or not smiles:
            continue
        if smiles in seen:
            continue
        seen.add(smiles)
        rows.append({
            'target': target_name,
            'molecule_chembl_id': mol_id,
            'canonical_smiles': smiles,
            'pchembl_value': a.get('pchembl_value', ''),
        })
    with open(out_csv, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=['target', 'molecule_chembl_id', 'canonical_smiles', 'pchembl_value'])
        w.writeheader()
        w.writerows(rows)
    return len(rows)

# 处理各靶点活性数据
pairs = [
    (os.path.join(LIGAND_DATA, 'A2A_活性分子_1uM.json'), 'A2A', os.path.join(LIGANDS, 'A2A_活性配体.csv')),
    (os.path.join(LIGAND_DATA, 'A2A_活性分子_1uM_p2.json'), 'A2A', None),
    (os.path.join(LIGAND_DATA, 'A2A_活性分子_1uM_p3.json'), 'A2A', None),
    (os.path.join(LIGAND_DATA, 'B2AR_活性分子_1uM.json'), 'B2AR', os.path.join(LIGANDS, 'B2AR_活性配体.csv')),
    (os.path.join(LIGAND_DATA, 'D3_活性分子_1uM.json'), 'D3', os.path.join(LIGANDS, 'D3_活性配体.csv')),
    (os.path.join(LIGAND_DATA, '5HT2B_活性分子_1uM.json'), '5HT2B', os.path.join(LIGANDS, '5HT2B_活性配体.csv')),
]

# 先把 A2A 的 3 页合并
import glob
all_a2a = []
for p in [os.path.join(LIGAND_DATA, 'A2A_活性分子_1uM.json'),
          os.path.join(LIGAND_DATA, 'A2A_活性分子_1uM_p2.json'),
          os.path.join(LIGAND_DATA, 'A2A_活性分子_1uM_p3.json')]:
    with open(p, encoding='utf-8') as f:
        all_a2a += json.load(f).get('activities', [])
a2a_merged = {'activities': all_a2a}
with open(os.path.join(LIGAND_DATA, 'A2A_活性分子_合并.json'), 'w', encoding='utf-8') as f:
    json.dump(a2a_merged, f, ensure_ascii=False)
print('A2A 合并后记录数:', len(all_a2a))

counts = {}
for jp, tn, csvp in pairs:
    if csvp:
        n = chembl_activities_to_smiles(jp, tn, csvp)
        counts[tn] = n
        print(f'{tn}: {n} 个唯一配体 → {csvp}')

# 处理药物库
with open(os.path.join(KNOWN_DRUGS, 'ChEMBL_已知药物样本.json'), encoding='utf-8') as f:
    d = json.load(f)
mols = d.get('molecules', [])
drug_rows = []
for m in mols:
    smiles = (m.get('molecule_structures') or {}).get('canonical_smiles')
    if smiles:
        drug_rows.append({
            'molecule_chembl_id': m.get('molecule_chembl_id'),
            'pref_name': m.get('pref_name', ''),
            'canonical_smiles': smiles,
        })
with open(os.path.join(KNOWN_DRUGS, 'ChEMBL_已知药物_100.csv'), 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.DictWriter(f, fieldnames=['molecule_chembl_id', 'pref_name', 'canonical_smiles'])
    w.writeheader()
    w.writerows(drug_rows)
print(f'已知药物库: {len(drug_rows)} 个已上市药物 → CSV')

print('\n=== 所有配体数据汇总 ===')
for k, v in counts.items():
    print(f'  {k}: {v} 个唯一配体')
