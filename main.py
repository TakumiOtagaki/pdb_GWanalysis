import os
import numpy as np
from Bio.PDB import PDBParser
import pandas as pd
# import multiprocessing
from multiprocessing import Pool
import sys
import re
import glob
# import ot
from scipy.special import logsumexp

# NUM_CPU = multiprocessing.cpu_count()   
NUM_CPU = 1
SKIP_DOWNLOADING = True
print(f"Number of CPUs: {NUM_CPU}")

def parse_and_compute_distance_matrix(args):
    pdb_file, output_directory = args
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('RNA', pdb_file)
    pdb_file_wo_suffix = re.sub(r'^pdb', '', os.path.basename(pdb_file))
    # もし RNA じゃなくて DNA なら、やめる
    # output_file = os.path.join(output_directory, f"{os.path.basename(pdb_file).split('.')[0]}_{*}_{*}.csv")
    # もし存在して、空じゃないのならやめる。
    output_file_regex = os.path.join(output_directory, f"{os.path.basename(pdb_file_wo_suffix).split('.')[0]}_*.csv")
    # return None
    # if os.path.exists(output_file_regex):
    if glob.glob(output_file_regex):
        # print(f"File {output_file_regex.split('/')[-1].split('.csv')[0]} already exists and is not empty. Skipping...")
        return None
    print(f"Processing {pdb_file}...")
    
    if not any(residue.get_resname() in ['A', 'C', 'G', 'U'] for model in structure for chain in model for residue in chain):
        return None
    for model in structure:
        for chain in model:
            atoms = [atom for residue in chain if residue.get_resname() in ['A', 'C', 'G', 'U']
                     for atom in residue if atom.get_name() in ["P",  "C3'",  "O3'"]]
            if atoms:
                num_atoms = len(atoms)
                dist_matrix = np.zeros((num_atoms, num_atoms))
                for i in range(num_atoms):
                    for j in range(i + 1, num_atoms):
                        distance = atoms[i].coord - atoms[j].coord
                        dist_matrix[i, j] = dist_matrix[j, i] = np.sqrt(np.sum(distance**2))
                
                # Save the distance matrix to a CSV file
                # pdb_file がたまに pdb から始まっていることがあるのでその場合は suffix から"pdb"という文字列を取り除く。
                # base = f"{os.path.basename(pdb_file).split('.')[0]}_{model.id}_{chain.id}"
                base = f"{os.path.basename(pdb_file_wo_suffix).split('.')[0].replace('pdb', '')}_{model.id}_{chain.id}"
                output_file = os.path.join(output_directory, f"{base}.csv")
                pd.DataFrame(dist_matrix).to_csv(output_file, index=False)
        break # Only process the first model
    return None

def my_sinkhorn(a, b, C, epsilon, max_iter=100000):
    n, m = C.shape
    f = np.ones(n)
    g = np.ones(m)
    for t in range(max_iter):
        if t % 1000 == 0:
            print(f"Iteration {t} / {max_iter}")
        # f は g から
        f = epsilon * np.log(a) - epsilon * logsumexp((-C + g.reshape(1, m)) / epsilon, axis=1) 
        # g は f から
        #         g = -epsilon * torch.logsumexp((-C + f.reshape(n, 1)) / epsilon, dim=0) + epsilon * torch.log(q)
        g = epsilon * np.log(b) - epsilon * logsumexp((-C + f.reshape(n, 1)) / epsilon, axis=0)
    P = np.exp((-C + f.reshape(n, 1) + g.reshape(1, m)) / epsilon)
    return P

def check_convergence(P, C, small_value = 1e-3):
    return np.linalg.norm(P.sum(axis=1) - 1) < small_value and np.linalg.norm(P.sum(axis=0) - 1) < small_value

def my_gromov_wasserstein_distance2(matrix1, matrix2, max_iter=100000, epsilon=1e-4, lambda_=1e1):
    n, m = matrix1.shape[0], matrix2.shape[0]
    a = np.ones(n) / n
    b = np.ones(m) / m

    P = np.ones((n, m)) / (n * m)
    for t in range(max_iter):
        C = -4 * matrix1 @ P @ matrix2.T + (epsilon - lambda_) * np.log(P)
        P = my_sinkhorn(a, b, C, lambda_)
        if check_convergence(P, C, epsilon):
            break
    return np.sum(P * C)
    
        



        
def compute_gw_distance(pair):
    matrix1_path, matrix2_path = pair
    print(f"Computing GW distance between {matrix1_path} and {matrix2_path}...")
    matrix1 = pd.read_csv(matrix1_path, header=0).values
    matrix2 = pd.read_csv(matrix2_path, header=0).values
    # numpy に
    print(matrix1.shape, matrix2.shape)
    
    
    # gw_distance = ot.gromov.gromov_wasserstein2(
    #     a=ot.unif(n), b=ot.unif(m), M=cost_matrix, p=2, q=2, log=True, max_iter=1000)
    gw_distance = my_gromov_wasserstein_distance2(matrix1, matrix2)
    return matrix1_path, matrix2_path, gw_distance

if __name__ == "__main__":
    pdb_directory = '/work/gs58/s58007/rna_GWclustering/data/pdb'
    output_directory = '/work/gs58/s58007/rna_GWclustering/data/internal_distmat'
    os.makedirs(output_directory, exist_ok=True)

    pdb_files = [os.path.join(pdb_directory, f) for f in os.listdir(pdb_directory) if f.endswith('.ent')]
    pdb_ids = [os.path.basename(f).split('.')[0] for f in pdb_files]
    if not SKIP_DOWNLOADING:
        with Pool() as pool:
            pool.map(parse_and_compute_distance_matrix, [(pdb_file, output_directory) for pdb_file in pdb_files])

    # List all distance matrix files
    matrix_files = [os.path.join(output_directory, f) for f in os.listdir(output_directory) if f.endswith('.csv')]
    pairs = [(matrix_files[i], matrix_files[j]) for i in range(len(matrix_files)) for j in range(i+1, len(matrix_files))]

    # Compute GW distances in parallel
    with Pool() as pool:
        results = pool.map(compute_gw_distance, pairs[0:1])
    # Cx_, Cy_ = matrix_files[0], matrix_files[1]
    # Cx = pd.read_csv(Cx_, header=0).values
    # Cy = pd.read_csv(Cy_, header=0).values
    # print(Cx.shape, Cy.shape)

    # You can further process the GW distance results as needed, e.g., save them to a file or analyze them.
    with open("/work/gs58/s58007/rna_GWclustering/data/gwdistmat.csv") as f:
        for result in results:
            m1, m2, d = result
            m1_base = os.path.basename(m1).split('.')[0]
            m2_base = os.path.basename(m2).split('.')[0]
            print(f"{m1_base},{m2_base},{d}")
            f.write(f"{m1_base},{m2_base},{d}\n")
