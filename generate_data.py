import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader

from tqdm import tqdm

from sage.combinat.set_partition_ordered import OrderedSetPartition
from sage.combinat.set_partition import SetPartitions
from sage.combinat.composition import Compositions
from sage.combinat.subset import Subsets
from sage.sets.set import Set
from bullock_early_jiang import hDOSP, sJ_XSr, XJ_XSr

def generate_dataset(n, k):
    dataset = []
    for Sr in hDOSP(n, k):
        for J in Subsets(n, k):
            label = XJ_XSr(Sr, J)
            if label not in [-1, 0, 1]:
                continue  # skip invalid labels
            S_python = [[int(e) for e in block] for block in Sr[0]]  # list of lists
            r_python = [int(x) for x in Sr[1]]
            J_python = [int(j) for j in J]
            dataset.append(((S_python, r_python, J_python), int(label)))
    return dataset

def generate_all_datasets(n_min, n_max):
    all_datasets = {}
    for n in tqdm(range(n_min, n_max + 1), desc="n-loop"):
        for k in range(2, n):
            dataset = generate_dataset(n, k)
            all_datasets[(n, k)] = dataset
    return all_datasets

def save_all_datasets(all_datasets, filename="multi_dataset.pt"):
    torch.save(all_datasets, filename)
    print(f"Datasets saved to {filename}")

def load_all_datasets(filename="multi_dataset.pt"):
    return torch.load(filename)

save_all_datasets(generate_all_datasets(3,7), '37.pt')