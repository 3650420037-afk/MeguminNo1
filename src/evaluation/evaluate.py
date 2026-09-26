import os
import argparse
import sys
from copy import deepcopy

# 本文件在 src/evaluation/: 先把上层的 src/ 放进 sys.path, 再取统一路径。
# 这样无论从仓库根、src/ 还是 src/evaluation/ 启动, `utils.* / evaluation.*` 都能导入。
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paths import OUTPUTS, TARGETS, src_on_path  # noqa: E402

src_on_path()

import torch
from tqdm.auto import tqdm
from rdkit.Chem.QED import qed

from utils.reconstruct import reconstruct_from_generated_with_edges
# 裸导入只在 cwd=evaluation 时成立; 补包内导入兼容 (见 scoring_func.py 同处说明)
try:
    from .sascorer import compute_sa_score
except ImportError:
    from sascorer import compute_sa_score
from evaluation.docking import *
from utils.misc import *
from evaluation.scoring_func import *



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('exp_name', type=str)
    parser.add_argument('--result_root', type=str, default=OUTPUTS)
    parser.add_argument('--protein_root', type=str, default=TARGETS)
    args = parser.parse_args()

    exp_dir = os.path.join(args.result_root, args.exp_name)
    save_path = os.path.join(exp_dir, 'samples_all.pt') # get_saved_filename(exp_dir))

    logger = get_logger('evaluate', exp_dir)
    logger.info(args)
    logger.info(save_path)

    samples = torch.load(save_path, map_location='cpu', weights_only=False)  # torch2.6 默认 weights_only=True 会拒载自产 pickle

    sim_with_train = SimilarityWithTrain()
    results = []
    for i, data in enumerate(tqdm(samples.finished, desc='All')):
        try:
            mol = reconstruct_from_generated_with_edges(data)
            vina_task = QVinaDockingTask.from_generated_data(data, protein_root=args.protein_root)
            results.append({
                'mol': mol,
                'vina': vina_task.run_sync(),
                'qed': qed(mol),
                'sa': compute_sa_score(mol),
                'lipinski': obey_lipinski(mol),
                'logp': get_logp(mol),
            })
        except Exception as e:
            logger.warning('Failed %d' % i)
            logger.warning(e)

    logger.info('Number of results: %d' % len(results))

    result_path = os.path.join(exp_dir, 'results.pt')
    torch.save(results, result_path)
