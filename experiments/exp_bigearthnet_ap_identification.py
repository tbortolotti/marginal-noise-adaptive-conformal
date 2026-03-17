import numpy as np
from sklearn.model_selection import train_test_split
import torch
import json
import os
import hydra
import pytorch_lightning as pl
from hydra.utils import instantiate
from hydra.core.global_hydra import GlobalHydra

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
import sys

sys.path.append("..")
sys.path.append("../third_party")

from cln import contamination
from cln.AP_identification import AnchorPointsIdentification
from cln.T_estimation import evaluate_estimate, TMatrixEstimation
from third_party import arc
from third_party.bigearthnet.datamodules.bigearthnet_datamodule import BigEarthNetDataModule
from third_party.bigearthnet.models.bigearthnet_module import BigEarthNetModule, BigEarthNetFeatureExtractor

# Define default parameters
exp_num = 1101
epsilon = 0.1
nu = 0
contamination_model = "uniform"
n_train1 = 1000
n_train2 = 1000
seed = 1

# Parse input parameters
if True:
    print('Number of arguments:', len(sys.argv), 'arguments.')
    print('Argument List:', str(sys.argv))
    if len(sys.argv) != 8:
        print("Error: incorrect number of parameters.")
        quit()
    sys.stdout.flush()
    exp_num = int(sys.argv[1])
    epsilon = float(sys.argv[2])
    nu = float(sys.argv[3])
    contamination_model = sys.argv[4]
    n_train1 = int(sys.argv[5])
    n_train2 = int(sys.argv[6])
    seed = int(sys.argv[7])

# Define other constant parameters
batch_size = n_train1 + n_train2
data_name = "bigearthnet"
K = 6
num_exp = 5

# Initialize noise contamination process
if contamination_model == "uniform":
    T = contamination.construct_T_matrix_simple(K, epsilon)
elif contamination_model == "block":
    T = contamination.construct_T_matrix_block(K, epsilon)
elif contamination_model == "RRB":
    T = contamination.construct_T_matrix_block_RR(K, epsilon, nu)
elif contamination_model == "random":
    T = contamination.construct_T_matrix_random(K, epsilon, random_state=seed)
else:
    print("Unknown contamination model!")
    sys.stdout.flush()
    exit(-1)

# BigEarthNet setup

# Clear any previous Hydra instances
if GlobalHydra.instance().is_initialized():
    GlobalHydra.instance().clear()
hydra.initialize(config_path="../third_party/bigearthnet/configs", version_base="1.2")
cfg = hydra.compose(config_name="config")
transforms = instantiate(cfg.transforms.obj)

# Load the clean label mapping (v1 -> v2 grouped labels)
with open('../third_party/bigearthnet/data/label_mapping.json', 'r') as f:
    label_mapping = json.load(f)

dataset_name = cfg.datamodule.dataset_name
file_path_train = os.path.join(
    '../third_party/bigearthnet/data',
    f'train_{dataset_name}.csv'
)
v1v2_corresp_train = pd.read_csv(file_path_train, header=0)

# Load the pre-trained BigEarthNet model and use it as a feature extractor.
# We call it in eval mode and extract the penultimate-layer embeddings,
# mirroring the role of ImageNetResNet18Features in the CIFAR-10 script.
black_box_model = BigEarthNetModule(cfg)
mod_path = os.path.join(cfg.out_directory.dir, "trained_model.pth")
black_box_model.load_state_dict(torch.load(mod_path))
black_box_model.eval()

feature_extractor = BigEarthNetFeatureExtractor(black_box_model)

# SVC black box, same as CIFAR-10 script
black_box_SVC = arc.black_boxes.SVC(clip_proba_factor=1e-5)

# Define dataloader
datamodule = BigEarthNetDataModule(cfg.datamodule.dataset_dir,
                                   cfg.datamodule.dataset_name,
                                   batch_size,
                                   cfg.datamodule.num_workers,
                                   transforms,
                                   label_mapping)
datamodule.setup()

# Add important parameters to table of results
header = pd.DataFrame({'data': [data_name], 'K': [K],
                       'n_train1': [n_train1], 'n_train2': [n_train2],
                       'epsilon': [epsilon], 'contamination': [contamination_model],
                       'seed': [seed]})

# Output file
outfile_prefix = "exp" + str(exp_num) + "/" + data_name + "_n" + str(batch_size)
outfile_prefix += "_eps" + str(epsilon) + "_" + contamination_model
outfile_prefix += "_nt1_" + str(n_train1) + "_nt2_" + str(n_train2) + "_seed" + str(seed)
print("Output file: {:s}.".format("results/" + outfile_prefix), end="\n")
sys.stdout.flush()


def run_experiment(random_state):
    print("\nRunning experiment in batch {:d}...".format(random_state))
    sys.stdout.flush()

    print("\nGenerating data...", end=' ')
    sys.stdout.flush()
    dataloader = datamodule.train_dataloader(seed=random_state)
    batch = next(iter(dataloader))
    X_batch = batch['data']

    # Retrieve the corresponding clean (v2-grouped) labels from the CSV,
    # using the same generator the datamodule used to shuffle its dataset.
    shuffled_csv = v1v2_corresp_train.iloc[datamodule.last_train_indices].reset_index(drop=True)
    batch_csv = shuffled_csv.iloc[:batch_size]
    Y_clean = batch_csv['v2-labels-grouped'].to_numpy()

    # OLD VERSION
    #generator = torch.Generator().manual_seed(datamodule.random_seed)
    #indices = torch.randperm(len(datamodule.train_dataset), generator=generator).tolist()
    #shuffled_csv = v1v2_corresp_train.iloc[indices].reset_index(drop=True)

    # Drop rows where clean label is NaN
    valid = ~np.isnan(Y_clean)
    X_batch = X_batch[valid]
    Y_clean = Y_clean[valid].astype(int)
    print(f"Done. Batch size after NaN removal: {len(Y_clean)}")
    sys.stdout.flush()

    # Feature extraction
    X_features = feature_extractor.transform(X_batch).numpy()
    del X_batch
    torch.cuda.empty_cache()

    # Artificial label contamination
    print("Generating contaminated labels...", end=' ')
    sys.stdout.flush()
    contamination_process = contamination.LinearContaminationModel(T, random_state=random_state + 1)
    Yt = contamination_process.sample_labels(Y_clean)
    print("Done.")
    sys.stdout.flush()

    # Train/test split into stage-1 and stage-2 sets
    X_features_train1, X_features_train2, _, Y_train2, Yt_train1, Yt_train2 = train_test_split(X_features, Y_clean, Yt, test_size=n_train2, random_state=random_state+2)
    del X_features

    # Methods
    methods = {
        "SVC": lambda: AnchorPointsIdentification(
            X_features_train1, Yt_train1,
            X_features_train2, Yt_train2, K,
            use_classifier=True, black_box=black_box_SVC,
            calibrate_gamma=True
        ),
        "IF": lambda: AnchorPointsIdentification(
            X_features_train1, Yt_train1,
            X_features_train2, Yt_train2, K,
            outlier_detection=True,
            outlier_detection_method="isolation_forest",
            selection="accuracy"
        ),
        "optimal": lambda: AnchorPointsIdentification(
            X_features_train1, Yt_train1,
            X_features_train2, Yt_train2, K,
            black_box=black_box_SVC,
            optimal_method=True,
            random_state=random_state + 3
        )
    }

    res_list = []

    # Baseline: clean-sample estimate using true Y/Yt correspondence
    T_method = TMatrixEstimation(Y_train2, Yt_train2, K, estimation_method="empirical")
    T_hat = T_method.get_estimate()
    performances = evaluate_estimate(T, T_hat, Y_train2, Y_train2, Yt_train2, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(
        Method='Clean sample',
        n_train1=n_train1, n_train2=n_train2,
        size=np.sum(Y_train2 != -1),
        **performances
    )
    res_list.append(res_update)

    # Anchor-point methods
    for method_name, method_func in methods.items():
        print(f"Applying {method_name} method...", end=' ')
        sys.stdout.flush()

        method = method_func()
        Ya_train2, _, _, _ = method.get_anchor_points()

        T_method = TMatrixEstimation(Ya_train2, Yt_train2, K, estimation_method="empirical_parametricRR")
        T_hat = T_method.get_estimate()

        print("Done.")
        sys.stdout.flush()

        performances = evaluate_estimate(T, T_hat, Y_train2, Ya_train2, Yt_train2, K)
        res_update = header.copy()
        res_update = res_update.assign(
            Method=method_name,
            n_train1=n_train1, n_train2=n_train2,
            size=np.sum(Ya_train2 != -1),
            **performances
        )
        res_list.append(res_update)

    res = pd.concat(res_list, ignore_index=True)
    return res


# Run all experiments
results_list = []
for batch in np.arange(1, num_exp + 1):
    # Seed formula: keeps random_state > 0 for all batch/seed combinations,
    # fixing the zero-seed bug present in the original bigearthnet script.
    res = run_experiment(1000 * seed + batch)
    results_list.append(res)

    # Save incrementally
    results = pd.concat(results_list, ignore_index=True)
    outfile = "results/" + outfile_prefix + ".txt"
    results.to_csv(outfile, index=False, float_format="%.5f")

print("\nFinished.\nResults written to {:s}\n".format(outfile))
sys.stdout.flush()