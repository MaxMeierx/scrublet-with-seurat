import sys
import os
import glob
import pandas as pd
import scipy.io
import scipy.sparse
import matplotlib.pyplot as plt
import numpy as np
import scrublet as scr

def process_files(input_dir):
    csv_files = glob.glob(f'{input_dir}*.csv')
    for i in csv_files:
        print(f'Processing {i}')
        df = pd.read_csv(i, index_col=0)
        df = df.apply(pd.to_numeric, errors='coerce')
        # 将 DataFrame 转换为稀疏矩阵
        counts_matrix = scipy.sparse.csr_matrix(df.values).T

        scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)

        # 准备输出数据
        out_df = pd.DataFrame({
            'doublet_scores': doublet_scores,
            'predicted_doublets': predicted_doublets
        }, index=df.T.index)

        # 构建输出目录和文件路径
        output_dir = os.path.join(input_dir, 'result')
        os.makedirs(output_dir, exist_ok=True)
        output_file_name = os.path.basename(i).replace('.csv', '_doublet.txt')
        output_file_path = os.path.join(output_dir, output_file_name)

        # 保存输出数据
        out_df.to_csv(output_file_path, index=True, header=True)
        print(f'Results saved to {output_file_path}')

if __name__ == '__main__':
    input_dir = './'  # 默认当前目录
    if len(sys.argv) > 1:
        input_dir = sys.argv[1]  # 允许从命令行参数指定目录
    process_files(input_dir)
  
