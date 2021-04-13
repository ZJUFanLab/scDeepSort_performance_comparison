from sklearn.svm import SVC
import os

import argparse
import pandas as pd
from sklearn.decomposition import PCA
from time import time
from scipy.sparse import csr_matrix, vstack
from pathlib import Path
import numpy as np


def get_map_dict(species_data_path: Path, tissue):
    map_df = pd.read_excel(species_data_path / 'map.xlsx')
    # {num: {test_cell1: {train_cell1, train_cell2}, {test_cell2:....}}, num_2:{}...}
    map_dic = dict()
    for idx, row in enumerate(map_df.itertuples()):
        if getattr(row, 'Tissue') == tissue:
            num = getattr(row, 'num')
            test_celltype = getattr(row, 'Celltype')
            train_celltype = getattr(row, '_5')
            if map_dic.get(getattr(row, 'num')) is None:
                map_dic[num] = dict()
                map_dic[num][test_celltype] = set()
            elif map_dic[num].get(test_celltype) is None:
                map_dic[num][test_celltype] = set()
            map_dic[num][test_celltype].add(train_celltype)
    return map_dic


def get_id_2_gene(gene_statistics_path, species_data_path, tissue, train_dir: str):
    if not gene_statistics_path.exists():
        data_path = species_data_path / train_dir
        data_files = data_path.glob(f'*{tissue}*_data.csv')
        genes = None
        for file in data_files:
            data = pd.read_csv(file, dtype=np.str, header=0).values[:, 0]
            if genes is None:
                genes = set(data)
            else:
                genes = genes | set(data)
        id2gene = list(genes)
        id2gene.sort()
        with open(gene_statistics_path, 'w', encoding='utf-8') as f:
            for gene in id2gene:
                f.write(gene + '\r\n')
    else:
        id2gene = []
        with open(gene_statistics_path, 'r', encoding='utf-8') as f:
            for line in f:
                id2gene.append(line.strip())
    return id2gene


def get_id_2_label(cell_statistics_path, species_data_path, tissue, train_dir: str):
    if not cell_statistics_path.exists():
        data_path = species_data_path / train_dir
        cell_files = data_path.glob(f'*{tissue}*_celltype.csv')
        cell_types = set()
        for file in cell_files:
            df = pd.read_csv(file, dtype=np.str, header=0)
            df['Cell_type'] = df['Cell_type'].map(str.strip)
            cell_types = set(df.values[:, 2]) | cell_types
            # cell_types = set(pd.read_csv(file, dtype=np.str, header=0).values[:, 2]) | cell_types
        id2label = list(cell_types)
        with open(cell_statistics_path, 'w', encoding='utf-8') as f:
            for cell_type in id2label:
                f.write(cell_type + '\r\n')
    else:
        id2label = []
        with open(cell_statistics_path, 'r', encoding='utf-8') as f:
            for line in f:
                id2label.append(line.strip())
    return id2label


def load_data(params):
    random_seed = params.random_seed
    dense_dim = params.dense_dim
    train = params.train_dataset
    test = params.test_dataset
    tissue = params.tissue

    proj_path = Path(__file__).parent.resolve().parent.resolve()
    species_data_path = proj_path / 'data' / params.species
    statistics_path = species_data_path / 'statistics'
    map_dict = get_map_dict(species_data_path, tissue)

    gene_statistics_path = statistics_path / (tissue + '_genes.txt')  # train+test gene
    cell_statistics_path = statistics_path / (tissue + '_cell_type.txt')  # train labels

    # generate gene statistics file
    id2gene = get_id_2_gene(gene_statistics_path, species_data_path, tissue, params.train_dir)
    # generate cell label statistics file
    id2label = get_id_2_label(cell_statistics_path, species_data_path, tissue, params.train_dir)

    train_num, test_num = 0, 0
    # prepare unified genes
    gene2id = {gene: idx for idx, gene in enumerate(id2gene)}
    num_genes = len(id2gene)
    # prepare unified labels
    num_labels = len(id2label)
    label2id = {label: idx for idx, label in enumerate(id2label)}
    print(f"totally {num_genes} genes, {num_labels} labels.")

    train_labels = []
    test_label_dict = dict()  # test label dict
    test_index_dict = dict()  # test-num: [begin-index, end-index]
    test_cell_id_dict = dict()  # test-num: ['c1', 'c2'...]
    # TODO
    matrices = []

    for num in train + test:
        start = time()
        if num in train:
            data_path = species_data_path / (params.train_dir + f'/{params.species}_{tissue}{num}_data.csv')
            type_path = species_data_path / (params.train_dir + f'/{params.species}_{tissue}{num}_celltype.csv')
        else:
            data_path = species_data_path / (params.test_dir + f'/{params.species}_{tissue}{num}_data.csv')
            type_path = species_data_path / (params.test_dir + f'/{params.species}_{tissue}{num}_celltype.csv')

        # load celltype file then update labels accordingly
        cell2type = pd.read_csv(type_path, index_col=0)
        cell2type.columns = ['cell', 'type']
        cell2type['type'] = cell2type['type'].map(str.strip)
        if num in train:
            cell2type['id'] = cell2type['type'].map(label2id)
            assert not cell2type['id'].isnull().any(), 'something wrong in celltype file.'
            train_labels += cell2type['id'].tolist()
        else:
            # test_labels += cell2type['type'].tolist()
            test_label_dict[num] = cell2type['type'].tolist()

        # load data file then update graph
        df = pd.read_csv(data_path, index_col=0)  # (gene, cell)
        if num in test:
            test_cell_id_dict[num] = list(df.columns)
        df = df.transpose(copy=True)  # (cell, gene)

        assert cell2type['cell'].tolist() == df.index.tolist()
        df = df.rename(columns=gene2id)
        # filter out useless columns if exists (when using gene intersection)
        col = [c for c in df.columns if c in gene2id.values()]
        df = df[col]
        print(f'Nonzero Ratio: {df.fillna(0).astype(bool).sum().sum() / df.size * 100:.2f}%')
        # maintain inter-datasets index for graph and RNA-seq values
        arr = df.to_numpy()
        row_idx, col_idx = np.nonzero(arr > params.threshold)  # intra-dataset index
        non_zeros = arr[(row_idx, col_idx)]  # non-zero values

        gene_idx = df.columns[col_idx].astype(int).tolist()  # gene_index
        info_shape = (len(df), num_genes)
        info = csr_matrix((non_zeros, (row_idx, gene_idx)), shape=info_shape)
        matrices.append(info)

        if num in train:
            train_num += len(df)
        else:
            test_index_dict[num] = list(range(train_num + test_num, train_num + test_num + len(df)))
            test_num += len(df)
        print(f'Costs {time() - start:.3f} s in total.')
    train_labels = np.array(list(map(int, train_labels)))

    # 2. create features
    sparse_feat = vstack(matrices).toarray()  # cell-wise  (cell, gene)
    test_feat_dict = dict()
    # transpose to gene-wise
    gene_pca = PCA(dense_dim, random_state=random_seed).fit(sparse_feat[:train_num].T)
    gene_feat = gene_pca.transform(sparse_feat[:train_num].T)
    gene_evr = sum(gene_pca.explained_variance_ratio_) * 100
    print(f'[PCA] Gene EVR: {gene_evr:.2f} %.')

    # do normalization
    sparse_feat = sparse_feat / (np.sum(sparse_feat, axis=1, keepdims=True) + 1e-6)
    # use weighted gene_feat as cell_feat
    cell_feat = sparse_feat.dot(gene_feat)  # [total_cell_num, d]
    train_cell_feat = cell_feat[:train_num]

    for num in test_label_dict.keys():
        test_feat_dict[num] = cell_feat[test_index_dict[num]]

    return num_labels, train_labels, train_cell_feat, map_dict, np.array(id2label, dtype=np.str), \
           test_label_dict, test_feat_dict, test_cell_id_dict


class Runner:
    def __init__(self, args):
        self.args = args
        self.prj_path = Path(__file__).parent.resolve().parent.resolve()
        self.num_labels, self.train_labels, self.train_cell_feat, self.map_dict, self.id2label, \
        self.test_label_dict, self.test_feat_dict, self.test_cell_id_dict = load_data(args)
        self.model = self.fit()

    def fit(self):
        model = SVC(random_state=self.args.random_seed, probability=True). \
            fit(self.train_cell_feat, self.train_labels)
        return model

    def evaluate(self):
        for num in self.args.test_dataset:
            score = self.model.predict_proba(self.test_feat_dict[num])  # [cell, class-num]
            pred_labels = []
            unsure_num, correct = 0, 0
            for pred, t_label in zip(score, self.test_label_dict[num]):
                pred_label = self.id2label[pred.argmax().item()]
                if pred_label in self.map_dict[num][t_label]:
                    correct += 1
                pred_labels.append(pred_label)

            acc = correct / score.shape[0]
            print(f'SVM-{self.args.species}-{self.args.tissue}-{num}-ACC: {acc:.5f}')
            self.save(num, pred_labels)

    def save(self, num, pred):
        label_map = pd.read_excel(self.prj_path / 'data' / 'celltype2subtype.xlsx',
                                  sheet_name=self.args.species, header=0,
                                  names=['species', 'old_type', 'new_type', 'new_subtype'])

        save_path = self.prj_path / self.args.save_dir
        if not save_path.exists():
            save_path.mkdir()

        label_map = label_map.fillna('N/A', inplace=False)
        oldtype2newtype = {}
        oldtype2newsubtype = {}
        for _, old_type, new_type, new_subtype in label_map.itertuples(index=False):
            oldtype2newtype[old_type] = new_type
            oldtype2newsubtype[old_type] = new_subtype
        if not os.path.exists(self.args.save_dir):
            os.mkdir(self.args.save_dir)

        df = pd.DataFrame({
            'index': self.test_cell_id_dict[num],
            'original label': self.test_label_dict[num],
            'cell type': [oldtype2newtype.get(p, p) for p in pred],
            'cell subtype': [oldtype2newsubtype.get(p, p) for p in pred]
        })
        df.to_csv(
            save_path / ('SVM_' + self.args.species + f"_{self.args.tissue}_{num}.csv"),
            index=False)
        print(f"output has been stored in {self.args.species}_{self.args.tissue}_{num}.csv")


def arg_parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--random_seed", type=int, default=10086)
    parser.add_argument("--train_dataset", nargs="+", required=True, type=int,
                        help="list of dataset id")
    parser.add_argument("--test_dataset", nargs="+", required=True, type=int,
                        help="list of dataset id")
    parser.add_argument("--species", default='mouse', type=str)
    parser.add_argument("--tissue", required=True, type=str)
    parser.add_argument("--train_dir", type=str, default='train')
    parser.add_argument("--test_dir", type=str, default='test')

    params = parser.parse_args()
    params.save_dir = 'result'
    params.dense_dim = 400
    params.threshold = 0
    return params


if __name__ == '__main__':
    params = arg_parse()

    runner = Runner(params)
    runner.evaluate()
