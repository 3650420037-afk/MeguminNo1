# -*- coding: utf-8 -*-
"""等价性证明: 改造后的 FrontierLayerVN 在默认容量下必须与官方实现**逐位相同**。

为什么必须证明
--------------
`frontier.py` / `maskfill.py` 被**训练和采样共用**。为了扩容 frontier_pred, 我把它的
构造改成了可配置(hidden_sca/hidden_vec/layers)。既然采样也走同一份代码, 那么
"默认配置(128/32/1 层)下行为与官方一字不差"就必须是**被证明的事实**, 而不是我的断言 ——
否则正在跑的药库生成与 A/B 对照都可能被悄悄改变。

证明三件事
----------
1. 参数名与形状集合完全一致(否则加载官方权重会失败或错位);
2. 把同一组权重分别灌进两个实现, 前向输出 `torch.equal`(**逐位相同**, 不是 allclose);
3. 用官方真实权重加载官方实现 vs 新实现, 在同一输入上输出逐位相同。

用法
----
    python src/scripts/verify_frontier_equivalence.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import ROOT, PRETRAINED, src_on_path  # noqa: E402

src_on_path()
import torch  # noqa: E402
from torch.nn import Module, Sequential  # noqa: E402
from models.invariant import GVPerceptronVN, GVLinear  # noqa: E402
from models.frontier import FrontierLayerVN  # noqa: E402

FAILS = []


def chk(name, ok, detail=""):
    print("  %s %s%s" % ("[OK]  " if ok else "[FAIL]", name, ("  -> " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)
    return ok


class OfficialFrontierLayerVN(Module):
    """官方实现的**原样复制**(改容量的改造之前), 用作对照基准。"""

    def __init__(self, in_sca, in_vec, hidden_dim_sca, hidden_dim_vec):
        super().__init__()
        self.net = Sequential(
            GVPerceptronVN(in_sca, in_vec, hidden_dim_sca, hidden_dim_vec),
            GVLinear(hidden_dim_sca, hidden_dim_vec, 1, 1)
        )

    def forward(self, h_att, idx_ligans):
        h_att_ligand = [h_att[0][idx_ligans], h_att[1][idx_ligans]]
        pred = self.net(h_att_ligand)
        pred = pred[0]
        return pred


def main():
    torch.manual_seed(0)
    IN_SCA, IN_VEC = 256, 64

    print("=" * 70)
    print("1/3 参数名与形状")
    print("=" * 70)
    old = OfficialFrontierLayerVN(IN_SCA, IN_VEC, 128, 32)
    new = FrontierLayerVN(IN_SCA, IN_VEC, 128, 32, n_layers=1)
    so = {k: tuple(v.shape) for k, v in old.state_dict().items()}
    sn = {k: tuple(v.shape) for k, v in new.state_dict().items()}
    chk("键集合一致", set(so) == set(sn),
        "仅旧有=%s 仅新有=%s" % (sorted(set(so) - set(sn))[:4], sorted(set(sn) - set(so))[:4]))
    diff = [k for k in so if k in sn and so[k] != sn[k]]
    chk("同名张量形状一致", not diff, str(diff[:4]))
    print("      张量数 %d, 参数量 %.6f M" % (len(so), sum(v.numel() for v in old.parameters()) / 1e6))

    print("=" * 70)
    print("2/3 同一权重下前向输出逐位相同")
    print("=" * 70)
    new.load_state_dict(old.state_dict())
    n_lig = 37
    h_att = [torch.randn(120, IN_SCA), torch.randn(120, IN_VEC, 3)]
    idx = torch.arange(n_lig)
    with torch.no_grad():
        yo, yn = old(h_att, idx), new(h_att, idx)
    chk("逐位相同 (torch.equal)", torch.equal(yo, yn),
        "最大差 %.3e" % float((yo - yn).abs().max()))
    chk("输出形状一致", tuple(yo.shape) == tuple(yn.shape), "%s vs %s" % (tuple(yo.shape), tuple(yn.shape)))

    print("=" * 70)
    print("3/3 用官方真实权重加载两种实现, 输出仍逐位相同")
    print("=" * 70)
    ck = PRETRAINED if os.path.exists(PRETRAINED) else os.path.join(ROOT, "models", "pretrained_Pocket2Mol.pt")
    if not chk("官方权重存在", os.path.exists(ck), ck):
        return 1
    sd = torch.load(ck, map_location="cpu", weights_only=False)["model"]
    fo = {k: v for k, v in sd.items() if k.startswith("frontier_pred.")}
    if not chk("官方权重含 frontier_pred", len(fo) > 0, "%d 个张量" % len(fo)):
        return 1
    o2 = OfficialFrontierLayerVN(IN_SCA, IN_VEC, 128, 32)
    n2 = FrontierLayerVN(IN_SCA, IN_VEC, 128, 32, n_layers=1)
    remap = lambda d: {k[len("frontier_pred."):]: v for k, v in d.items()}
    o2.load_state_dict(remap(fo))
    n2.load_state_dict(remap(fo))
    with torch.no_grad():
        yo2, yn2 = o2(h_att, idx), n2(h_att, idx)
    chk("官方权重下逐位相同", torch.equal(yo2, yn2),
        "最大差 %.3e" % float((yo2 - yn2).abs().max()))
    print("      官方 frontier_pred 参数量 %.6f M" % (sum(v.numel() for v in o2.parameters()) / 1e6))

    print("=" * 70)
    if FAILS:
        print("结果: FAIL -> %s" % FAILS)
        return 1
    print("结果: PASS —— 默认容量(128/32/1层)下与官方实现逐位等价, 采样行为未被改变")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
