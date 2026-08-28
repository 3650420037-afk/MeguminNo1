# -*- coding: utf-8 -*-
"""方向1 引导束搜索: 对未完成候选计算化学引导分数(QED/SA), 用于束排序.
零侵入: 不修改任何原有函数行为, 仅新增; 默认 config.sample.guided 不存在时走原路径."""
import sys, os, types, pickle, importlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED

# --- rdkit.six 兼容 shim: 新版 RDKit 移除了 rdkit.six, 而 evaluation/sascorer.py 依赖它 ---
try:
    importlib.import_module("rdkit.six")
except ImportError:
    _six = types.ModuleType("rdkit.six")
    _moves = types.ModuleType("rdkit.six.moves")
    _moves.cPickle = pickle
    _six.moves = _moves
    _six.iteritems = lambda d, **kw: iter(d.items(**kw))
    sys.modules.setdefault("rdkit.six", _six)
    sys.modules.setdefault("rdkit.six.moves", _moves)

from evaluation.sascorer import calculateScore

RDLogger.DisableLog("rdApp.*")
_SCORE_CACHE = {}

def build_partial_mol(data):
    """由未完成候选(部分分子)构建宽松 RDKit 分子; 失败返回 None."""
    try:
        pos = data.ligand_context_pos.numpy()
        els = data.ligand_context_element.numpy().tolist()
        bi = data.ligand_context_bond_index.numpy().tolist()
        bt = data.ligand_context_bond_type.numpy().tolist()
    except Exception:
        return None
    m = Chem.RWMol()
    for z in els:
        m.AddAtom(Chem.Atom(int(z)))
    conf = Chem.Conformer(len(els))
    for i, p in enumerate(pos):
        conf.SetAtomPosition(i, (float(p[0]), float(p[1]), float(p[2])))
    m.AddConformer(conf)
    btmap = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
    for k, t in enumerate(bt):
        i, j = bi[0][k], bi[1][k]
        if i < j and int(t) in btmap:
            try:
                m.AddBond(i, j, btmap[int(t)])
            except Exception:
                return None
    mol = m.GetMol()
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE
                         ^ Chem.SANITIZE_SETAROMATICITY ^ Chem.SANITIZE_PROPERTIES,
                         catchErrors=True)
    except Exception:
        return None
    return mol

def chem_score(data, qed_w=1.0, sa_w=1.0, sa_max=10.0, fallback=0.5):
    """化学引导分数, 区间约 [0, qed_w+sa_w], 越大越类药/越易合成."""
    mol = build_partial_mol(data)
    if mol is None:
        return fallback * (qed_w + sa_w)
    try:
        smi = Chem.MolToSmiles(mol)
    except Exception:
        return fallback * (qed_w + sa_w)
    if smi in _SCORE_CACHE:
        return _SCORE_CACHE[smi]
    try:
        q = QED.qed(Chem.MolFromSmiles(smi) or mol)
    except Exception:
        q = fallback
    try:
        s = calculateScore(mol)  # 1(易合成)~10(难)
    except Exception:
        s = sa_max * (1.0 - fallback)
    score = q * qed_w + (1.0 - min(s, sa_max) / sa_max) * sa_w
    if len(_SCORE_CACHE) < 200000:
        _SCORE_CACHE[smi] = score
    return score

def logp_to_rank_prob_guided(logp, chem_scores, weight=1.0, lam=1.0):
    """引导版束排序: prob ∝ (exp(sum logp)+1) * weight * exp(lam * chem_score).
    chem_scores 与 logp 等长; lam=0 时退化为原始排序."""
    logp_sum = np.array([np.sum(l) for l in logp])
    prob = np.exp(logp_sum) + 1
    prob = prob * np.array(weight)
    prob = prob * np.exp(lam * np.asarray(chem_scores, dtype=float))
    total = prob.sum()
    if not np.isfinite(total) or total <= 0:
        return np.ones(len(prob)) / len(prob)
    return prob / total
