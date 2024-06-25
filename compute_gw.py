from modules.gw import my_gromov_wasserstein_distance2, my_sinkhorn, check_convergence

import os
import numpy as np
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='Compute Gromov-Wasserstein distance between RNA structures')
    parser.add_argument('matrix1 path', type=str, help='Path to the first distance matrix')
    parser.add_argument('matrix2 path', type=str, help='Path to the second distance matrix')
    parser.add_argument('output_directory', type=str, help='Directory to save distance matrices')
    parser.add_argument('max_iter', type=int, default=10000, help='Maximum number of iterations')
    return parser.parse_args()



def main():
    args = parse_args()

    m1 = pd.read_csv(args.matrix1_path, header=0).values
    m2 = pd.read_csv(args.matrix2_path, header=0).values
    gw_distance = my_gromov_wasserstein_distance2(m1, m2, max_iter=args.max_iter)
    
    with open(os.path.join(args.output_directory, f"{os.path.basename(args.matrix1_path)}_{os.path.basename(args.matrix2_path)}.txt"), 'w') as f:
        f.write(f"{m1},{m2},{gw_distance}")
    

if __name__ == "__main__":
    main()