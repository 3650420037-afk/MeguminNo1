# -*- coding: utf-8 -*-
"""get_tri_edges 向量化改造的一致性测试 + 性能对比
验证: 新实现(向量化/回退) 与原版逐 query 循环实现在多组随机输入下输出逐位一致。
"""
import os, sys, time
sys.path.insert(0, r"D:\MMModel\Pocket2Mol")
import torch
from models.maskfill import MaskFillModelVN

def original_tri_edges(edge_index_query, pos_query, idx_ligand, ligand_bond_index, ligand_bond_type, device="cpu"):
    """原版实现逐字复刻 (device 参数化, 逻辑不变)"""
    row, col = edge_index_query
    acc_num_edges = 0
    li, lj = [], []
    for node in torch.arange(pos_query.size(0)):
        num_edges = (row == node).sum()
        index_edge_i = torch.arange(num_edges, dtype=torch.long, device=device) + acc_num_edges
        index_edge_i, index_edge_j = torch.meshgrid(index_edge_i, index_edge_i, indexing=None)
        li.append(index_edge_i.flatten()); lj.append(index_edge_j.flatten())
        acc_num_edges += num_edges
    index_real_cps_edge_i = torch.cat(li, dim=0)
    index_real_cps_edge_j = torch.cat(lj, dim=0)
    node_a = col[index_real_cps_edge_i]; node_b = col[index_real_cps_edge_j]
    n_context = len(idx_ligand)
    adj = (torch.zeros([n_context, n_context], dtype=torch.long) - torch.eye(n_context, dtype=torch.long)).to(device)
    adj[ligand_bond_index[0], ligand_bond_index[1]] = ligand_bond_type
    t = adj[node_a, node_b]
    feat = (t.view([-1, 1]) == torch.tensor([[-1, 0, 1, 2, 3]]).to(device)).long()
    return torch.stack([index_real_cps_edge_i, index_real_cps_edge_j], 0), torch.stack([node_a, node_b], 0), feat

def make_case(n_query, n_context, seed):
    g = torch.Generator().manual_seed(seed)
    # query x context 全连接 (与 query_position 的 meshgrid(arange(n_query), arange(n_context)) 一致)
    mesh = torch.stack(torch.meshgrid(torch.arange(n_query), torch.arange(n_context), indexing="ij"), dim=0)
    edge_index_query = mesh.reshape(2, -1)
    pos_query = torch.randn(n_query, 3, generator=g)
    idx_ligand = torch.arange(n_context)
    # 随机键: 30% 概率加一条键
    nb = int(0.3 * n_context)
    if nb > 0:
        b0 = torch.randint(0, n_context, (nb,), generator=g)
        b1 = torch.randint(0, n_context, (nb,), generator=g)
        bt = torch.randint(1, 4, (nb,), generator=g)
        bond_index = torch.stack([b0, b1]); bond_type = bt
    else:
        bond_index = torch.empty(2, 0, dtype=torch.long); bond_type = torch.empty(0, dtype=torch.long)
    return edge_index_query, pos_query, idx_ligand, bond_index, bond_type

def check(a, b, name):
    if a.shape != b.shape:
        return "%s SHAPE MISMATCH %s vs %s" % (name, tuple(a.shape), tuple(b.shape))
    if not torch.equal(a, b):
        diff = (a != b).sum().item()
        return "%s VALUE MISMATCH (%d diffs)" % (name, diff)
    return None

fails = []
cases = []
for n_context in (1, 2, 5, 12, 30):
    for n_query in (1, 3, 8):
        cases.append((n_query, n_context))
print("=== 一致性测试 (规则全连接分组) ===")
for i, (nq, nc) in enumerate(cases):
    eiq, pq, il, bi, bt = make_case(nq, nc, seed=100 + i)
    o = original_tri_edges(eiq, pq, il, bi, bt)
    n = MaskFillModelVN.get_tri_edges(MaskFillModelVN, eiq, pq, il, bi, bt)
    for k, nm in enumerate(("edge_index_pair", "tri_edge_index", "tri_edge_feat")):
        err = check(o[k], n[k], "nq=%d nc=%d %s" % (nq, nc, nm))
        if err: fails.append(err)
    print("  nq=%2d nc=%2d sizes=%s/%s featshape=%s" % (nq, nc, tuple(o[0].shape), tuple(n[0].shape), tuple(o[2].shape)))

# 非常规分组 -> 回退路径
print("=== 回退路径测试 (非规则分组) ===")
g = torch.Generator().manual_seed(7)
n_context = 6
rows = torch.cat([torch.zeros(4, dtype=torch.long), torch.ones(2, dtype=torch.long), torch.full((3,), 2)])
cols = torch.randint(0, n_context, rows.shape, generator=g)
eiq = torch.stack([rows, cols])
pq = torch.randn(3, 3); il = torch.arange(n_context)
bi = torch.stack([torch.tensor([0, 1]), torch.tensor([1, 2])]); bt = torch.tensor([1, 2])
o = original_tri_edges(eiq, pq, il, bi, bt)
n = MaskFillModelVN.get_tri_edges(MaskFillModelVN, eiq, pq, il, bi, bt)
for k, nm in enumerate(("edge_index_pair", "tri_edge_index", "tri_edge_feat")):
    err = check(o[k], n[k], "fallback %s" % nm)
    if err: fails.append(err)
print("  回退路径 sizes=%s" % (tuple(n[0].shape),))

# 性能对比 (beam=100 场景的典型规模)
print("=== 性能对比 (n_query=100, n_context=30) ===")
eiq, pq, il, bi, bt = make_case(100, 30, seed=1)
for name, fn in (("original", original_tri_edges), ("vectorized", lambda *a: MaskFillModelVN.get_tri_edges(MaskFillModelVN, *a))):
    t0 = time.time()
    for _ in range(5):
        fn(eiq, pq, il, bi, bt)
    dt = (time.time() - t0) / 5
    print("  %-11s %.4f s/call" % (name, dt))

print()
print("结果:", "ALL PASS" if not fails else "FAILURES:\n  " + "\n  ".join(fails))
