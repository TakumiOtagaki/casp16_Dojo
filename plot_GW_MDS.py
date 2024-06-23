import os
import numpy as np
from Bio.PDB import PDBParser
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
import ot

import multiprocessing
# CPU の数
num_cores = 30

print("Hello, World!")
# PDBファイルのディレクトリ
pdb_dir = '/large/otgk/casp/casp16/R1224s2/R1224s_2_pdb100'

# PDBファイルのリストを取得
pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]

# RNAの主鎖の原子名
main_chain_atoms = ['C5\'', 'C4\'', 'C3\'']

# 残基間距離を計算する関数
def calculate_distances(coords):
    n = len(coords)
    distances = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            distances[i, j] = np.linalg.norm(coords[i] - coords[j])
            distances[j, i] = distances[i, j]
    return distances

# PDBファイルから主鎖の座標を抽出する関数
def extract_main_chain_coords(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('', pdb_file)
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() in main_chain_atoms:
                        coords.append(atom.get_coord())
    return np.array(coords)




# for pdb_file in pdb_files:
#     print(f'Processing {pdb_file}..., {count}/{len(pdb_files)}')
#     count += 1
#     coords = extract_main_chain_coords(os.path.join(pdb_dir, pdb_file))
#     distance_matrix = calculate_distances(coords)
#     distance_matrices.append(distance_matrix)

# Gromov-Wasserstein距離を計算する関数
def calculate_gw_distance(D1, D2):
    print("GW")
    p = ot.unif(D1.shape[0])
    q = ot.unif(D2.shape[0])
    gw_distance = ot.gromov.gromov_wasserstein2(D1, D2, p, q, 'square_loss', max_iter=100000)
    return gw_distance

def rmsd(coords1, coords2):
    n = len(coords1)
    assert n == len(coords2)
    return np.sqrt(np.sum(np.linalg.norm(coords1 - coords2, axis=1)**2) / n)


if __name__ == "__main__":
    # 主鎖座標を抽出し、距離行列を計算
    distance_matrices = []
    count = 0
    with multiprocessing.Pool(num_cores) as pool:
        distance_matrices = pool.map(extract_main_chain_coords, [os.path.join(pdb_dir, pdb_file) for pdb_file in pdb_files])
        distance_matrices = pool.map(calculate_distances, distance_matrices)

    # 全ての距離行列間のGromov-Wasserstein距離を計算
    n = len(distance_matrices)
    gw_distances = np.zeros((n, n))
    print("Calculating Gromov-Wasserstein distances...")

    with multiprocessing.Pool(num_cores) as pool:
        results = pool.starmap(calculate_gw_distance, [(distance_matrices[i], distance_matrices[j]) for i in range(n) for j in range(i+1, n)])
        for i in range(n):
            for j in range(i+1, n):
                gw_distances[i, j] = results.pop(0)
                gw_distances[j, i] = gw_distances[i, j]


    # MDSを用いて次元削減と可視化
    mds = MDS(n_components=2, dissimilarity='precomputed')
    coords_2d = mds.fit_transform(gw_distances)

    # 可視化
    plt.figure(figsize=(10, 7))
    plt.scatter(coords_2d[:, 0], coords_2d[:, 1])
    for i, pdb_file in enumerate(pdb_files):
        plt.annotate(pdb_file, (coords_2d[i, 0], coords_2d[i, 1]))
    plt.xlabel('MDS Dimension 1')
    plt.ylabel('MDS Dimension 2')
    plt.title('MDS of RNA Structures based on Gromov-Wasserstein Distances')
    plt.savefig("/large/otgk/casp/casp16/R1224s2/figures/GW_MDS.png")

    # RMSD でもやる
    rmsd_distances = np.zeros((n, n))
    print("Calculating RMSD distances...")
    with multiprocessing.Pool(num_cores) as pool:
        results = pool.starmap(rmsd, [(distance_matrices[i], distance_matrices[j]) for i in range(n) for j in range(i+1, n)])
        for i in range(n):
            for j in range(i+1, n):
                rmsd_distances[i, j] = results.pop(0)
                rmsd_distances[j, i] = rmsd_distances[i, j]

    # MDSを用いて次元削減と可視化
    mds = MDS(n_components=2, dissimilarity='precomputed')
    coords_2d = mds.fit_transform(rmsd_distances)

    # 可視化
    plt.figure(figsize=(10, 7))
    plt.scatter(coords_2d[:, 0], coords_2d[:, 1])
    for i, pdb_file in enumerate(pdb_files):
        plt.annotate(pdb_file, (
            coords_2d[i, 0], coords_2d[i, 1]))
    plt.xlabel('MDS Dimension 1')
    plt.ylabel('MDS Dimension 2')
    plt.title('MDS of RNA Structures based on RMSD Distances')
    plt.savefig("/large/otgk/casp/casp16/R1224s2/figures/RMSD_MDS.png")
    print("Done!")
    
