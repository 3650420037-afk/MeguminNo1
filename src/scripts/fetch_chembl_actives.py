# -*- coding: utf-8 -*-
"""大批量抓取 ChEMBL 活性分子, 用于把"每口袋样本密度"从百级提到千级。

背景: 已实测证明**每口袋样本密度是决定微调成败的关键**(1.04/口袋时分子永不终止、
完成数 0-30; 40-127/口袋时完成数 61-64)。当前 `data/ligands/*_活性配体.csv` 每个靶点
只有 99-191 个分子, 本脚本把效价阈值放宽到 10 uM(pchembl>=5)并完整分页, 目标千级。

用法:
    python src/scripts/fetch_chembl_actives.py --pchembl-min 5 --limit-targets 4
产物:
    data/ligands/<靶点>_活性配体_full.csv   (target, molecule_chembl_id, canonical_smiles, pchembl_value)
    data/ligand_data/chembl_<靶点>_raw.json (原始分页记录, 供溯源)
"""
import argparse
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import DATA, LIGANDS, LIGAND_DATA, ensure_dir  # noqa: E402

API = "https://www.ebi.ac.uk/chembl/api/data"
# 靶点 -> 人类(Homo sapiens) ChEMBL 靶点 ID。
# 用固定 ID 而不是按名称猜: 名称检索会先命中小鼠/大鼠等同源靶点(实测), 且各受体
# 在 ChEMBL 里的 pref_name 写法不一(如 "Adenosine receptor A2a")。脚本仍会调 API
# 回读 pref_name/organism 做校验, 对不上就报警。
TARGETS = {
    "A2A": "CHEMBL251",     # Adenosine receptor A2a (human)
    "B2AR": "CHEMBL210",    # Beta-2 adrenergic receptor (human)
    "D3": "CHEMBL234",      # Dopamine D3 receptor (human)
    "5HT2B": "CHEMBL1833",  # 5-hydroxytryptamine receptor 2B (human)
}


def get_json(url, tries=6, sleep=15.0):
    """带长退避的 GET。EBI ChEMBL 在被连续大量拉取时会返回 HTTP 500 或直接超时
    (实测: 连续拉 9 页 x1000 条后, 其余靶点全部 500), 因此用分钟级退避而不是秒级。"""
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json",
                                                       "User-Agent": "7-eonmol/1.0"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:
            if i == tries - 1:
                print("    请求失败(%d/%d): %s" % (i + 1, tries, e), flush=True)
                return None
            wait = sleep * (i + 1)
            print("    第 %d 次失败(%s), %.0f 秒后重试..." % (i + 1, str(e)[:50], wait),
                  flush=True)
            time.sleep(wait)
    return None


def resolve_target(tid):
    """回读该 ChEMBL 靶点的真实名称与物种, 用于校验硬编码 ID 是否正确。"""
    d = get_json("%s/target/%s.json" % (API, tid))
    if not d:
        return None, None
    return d.get("pref_name"), d.get("organism")


def fetch_activities(tid, pchembl_min, max_pages=40):
    """分页抓取该靶点的活性记录, 返回去重后的 [(mol_id, smiles, pchembl)] 与原始页。"""
    best, raw_pages = {}, []
    offset = 0
    limit = 1000
    for page in range(max_pages):
        url = ("%s/activity.json?target_chembl_id=%s&pchembl_value__gte=%s"
               "&limit=%d&offset=%d" % (API, tid, pchembl_min, limit, offset))
        d = get_json(url)
        if not d:
            break
        acts = d.get("activities", [])
        raw_pages.append(acts)
        for a in acts:
            smi = (a.get("canonical_smiles") or "").strip()
            if not smi:
                continue
            mid = a.get("molecule_chembl_id")
            try:
                pv = float(a.get("pchembl_value"))
            except (TypeError, ValueError):
                continue
            if mid not in best or pv > best[mid][1]:
                best[mid] = (smi, pv)
        nxt = d.get("page_meta", {}).get("next")
        print("    第 %d 页: %d 条 | 累计唯一分子 %d" % (page + 1, len(acts), len(best)),
              flush=True)
        if not nxt:
            break
        offset += limit
        time.sleep(2.0)
    return [(m, s, p) for m, (s, p) in best.items()], raw_pages


def main():
    ap = argparse.ArgumentParser(description="批量抓取 ChEMBL 活性分子")
    ap.add_argument("--pchembl-min", type=float, default=5.0,
                    help="效价阈值(pchembl); 6 约等于 1 uM, 5 约等于 10 uM。默认 5 (更宽)")
    ap.add_argument("--max-pages", type=int, default=40, help="每个靶点最多抓多少页")
    ap.add_argument("--targets", nargs="+", default=list(TARGETS))
    args = ap.parse_args()

    ensure_dir(LIGAND_DATA)
    summary = []
    for name in args.targets:
        tid = TARGETS.get(name)
        if not tid:
            print("[%s] 无对应靶点 ID, 跳过" % name)
            continue
        print("=" * 70)
        print("[%s] ChEMBL 靶点 %s, 校验中..." % (name, tid), flush=True)
        real, org = resolve_target(tid)
        if not real:
            print("  校验失败(API 无响应), 跳过")
            continue
        print("  -> %s | 物种 %s" % (real, org), flush=True)
        if not (org or "").lower().startswith("homo"):
            print("  !! 警告: 该靶点不是人类来源, 仍继续但请核对")
        rows, raws = fetch_activities(tid, args.pchembl_min, args.max_pages)
        if not rows:
            print("  未取到记录")
            continue
        rawp = os.path.join(LIGAND_DATA, "chembl_%s_raw.json" % name)
        with open(rawp, "w", encoding="utf-8") as f:
            json.dump({"target_chembl_id": tid, "pref_name": real,
                       "pchembl_min": args.pchembl_min, "pages": raws}, f, ensure_ascii=False)
        out = os.path.join(LIGANDS, "%s_活性配体_full.csv" % name)
        with open(out, "w", encoding="utf-8-sig", newline="") as f:
            w = csv.writer(f)
            w.writerow(["target", "molecule_chembl_id", "canonical_smiles", "pchembl_value"])
            for mid, smi, pv in sorted(rows, key=lambda x: -x[2]):
                w.writerow([name, mid, smi, pv])
        pv = sorted(r[2] for r in rows)
        print("  [%s] 唯一活性分子 %d 个 (pchembl 中位 %.2f, 范围 %.1f-%.1f) -> %s"
              % (name, len(rows), pv[len(pv) // 2], pv[0], pv[-1],
                 os.path.relpath(out, os.path.dirname(LIGANDS))), flush=True)
        summary.append((name, tid, len(rows)))

    print("=" * 70)
    print("汇总:", ", ".join("%s=%d" % (a, c) for a, _, c in summary))
    print("合计 %d 个活性分子" % sum(c for _, _, c in summary))


if __name__ == "__main__":
    main()
