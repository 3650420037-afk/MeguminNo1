import os
import sys
import time
import random
import logging
import torch
import numpy as np
import yaml
from easydict import EasyDict
from logging import Logger
from tqdm.auto import tqdm

# ============================================================
# PyTorch 2.6+ 兼容性补丁
# 新版 PyTorch 的 torch.load 默认 weights_only=True,而 2022 年的
# 旧权重文件中包含 easydict.EasyDict、numpy 等对象,会导致
# "Weights only load failed" 错误。此处显式放行这些类型。
# ============================================================
try:
    torch.serialization.add_safe_globals([EasyDict])
    import numpy as _np
    torch.serialization.add_safe_globals([_np.ndarray, _np.dtype, _np.float32,
                                          _np.float64, _np.int64, _np.int32,
                                          _np.bool_])
except Exception:
    pass

# ============================================================
# PyG 2.8 + Windows 兼容性补丁
# PyG 2.8 的 knn_graph 依赖 pyg-lib 的 knn 算子,但 Windows 版
# pyg-lib 编译时未包含该算子,导致 WITH_KNN=False 而报错。
# 此处强制启用 KNN,并让 torch_geometric 的回退实现使用
# torch_cluster 的 knn_graph(功能完全一致)。
# ============================================================
try:
    import torch_geometric.typing as _tg_typing
    if not _tg_typing.WITH_KNN:
        _tg_typing.WITH_KNN = True
except Exception:
    pass


class BlackHole(object):
    def __setattr__(self, name, value):
        pass
    def __call__(self, *args, **kwargs):
        return self
    def __getattr__(self, name):
        return self


def load_config(path):
    # utf-8-sig: 同时兼容带 BOM 与不带 BOM 的 YAML (PowerShell/记事本写出的配置可能带
    # BOM, 旧写法会把 BOM 并进首个键名, 使 config.model 变成 config['\ufeffmodel'],
    # 进而报 AttributeError: 'EasyDict' object has no attribute 'model')
    with open(path, 'r', encoding='utf-8-sig') as f:
        return EasyDict(yaml.safe_load(f))


def _utf8_stream():
    """给 StreamHandler 一个 UTF-8 包装的 stdout。

    为什么不直接用默认流: 中文 Windows 上默认流按 GBK 编码, 而文件日志已改为 UTF-8,
    于是"同一份日志两种编码"; 且控制台重定向到文件时中文会变成乱码,
    与文件日志内容不一致, 给排障制造假信息。
    取不到 fileno(如某些 IDE/嵌入式环境)则返回 None, 由调用方回退到默认流。
    """
    try:
        return open(sys.stdout.fileno(), mode='w', encoding='utf-8',
                    errors='replace', buffering=1, closefd=False)
    except Exception:
        return None


def get_logger(name, log_dir=None):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('[%(asctime)s::%(name)s::%(levelname)s] %(message)s')

    _s = _utf8_stream()
    stream_handler = logging.StreamHandler(_s) if _s is not None else logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    if log_dir is not None:
        # encoding='utf-8' 是**必须**的: 不指定时 FileHandler 用 locale 默认编码
        # (中文 Windows = GBK), 日志里的中文行就不是合法 UTF-8 ——
        # grep / read / CI 会看到乱码甚至直接报 "line is not valid UTF-8",
        # 中英混编的日志无法被工具稳定解析(实测差点因此漏看关键的诊断行)。
        file_handler = logging.FileHandler(os.path.join(log_dir, 'log.txt'),
                                           encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def get_new_log_dir(root='./logs', prefix='', tag=''):
    fn = time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())
    if prefix != '':
        fn = prefix + '_' + fn
    if tag != '':
        fn = fn + '_' + tag
    log_dir = os.path.join(root, fn)
    os.makedirs(log_dir)
    return log_dir


def seed_all(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)


def log_hyperparams(writer, args):
    from torch.utils.tensorboard.summary import hparams
    vars_args = {k:v if isinstance(v, str) else repr(v) for k, v in vars(args).items()}
    exp, ssi, sei = hparams(vars_args, {})
    writer.file_writer.add_summary(exp)
    writer.file_writer.add_summary(ssi)
    writer.file_writer.add_summary(sei)


def int_tuple(argstr):
    return tuple(map(int, argstr.split(',')))


def str_tuple(argstr):
    return tuple(argstr.split(','))

def unique(x, dim=None):
    """Unique elements of x and indices of those unique elements
    https://github.com/pytorch/pytorch/issues/36748#issuecomment-619514810

    e.g.

    unique(tensor([
        [1, 2, 3],
        [1, 2, 4],
        [1, 2, 3],
        [1, 2, 5]
    ]), dim=0)
    => (tensor([[1, 2, 3],
                [1, 2, 4],
                [1, 2, 5]]),
        tensor([0, 1, 3]))
    """
    unique, inverse = torch.unique(
        x, sorted=True, return_inverse=True, dim=dim)
    perm = torch.arange(inverse.size(0), dtype=inverse.dtype,
                        device=inverse.device)
    inverse, perm = inverse.flip([0]), perm.flip([0])
    return unique, inverse.new_empty(unique.size(dim)).scatter_(0, inverse, perm)