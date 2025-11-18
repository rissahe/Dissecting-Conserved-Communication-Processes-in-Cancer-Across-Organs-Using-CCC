import scanpy as sc
import liana as li
import pandas as pd
import os

out_path = "/data/sm979441/RCC_results/liana"
data = sc.read_h5ad("/home/larissa/Documents/Masterarbeit/objects/kidney_cancer_processed.h5ad")
data = data[(data.obs["Predicted_subgroups_PILOT"] == "Tumor_early_tarjectory")  |  (data.obs["Predicted_subgroups_PILOT"] == "Normal")]
data.raw = data

for i in set(data.obs['condition']):
    print(i)
    lr = li.method.cellphonedb(data[data.obs['sample'] == i],
                               groupby='cell_type',
                               expr_prop=0.1,
                               verbose=True,
                               resource_name='consensus', inplace=False)
    # lr.to_csv(f"{i}_lr_liana_consensus.csv")
    lr.to_csv(f"{out_path}{i}_lr_liana_consensus.csv")

data = {}
for i in os.listdir(out_path):
    if i.endswith('lr_liana_consensus.csv'):
        evfull = pd.read_csv(out_path+i)
        evfull = evfull.loc[:, ['ligand', 'receptor', 'source', 'target', 'lr_means', 'cellphone_pvals']]
        evfull['type_gene_A'] = 'Ligand'
        evfull['type_gene_B'] = 'Receptor'
        evfull['gene_A'] = evfull['ligand']
        evfull['gene_B'] = evfull['receptor']
        evfull['MeanLR'] = evfull['lr_means']
        k = i[0:i.find('_lr_')]
        evfull.loc[list(evfull.cellphone_pvals.to_numpy() <= 0.01), :].to_csv(f'{out_path}{k}_lr_ready.csv')
        data[k] = os.path.abspath(f'{out_path}{k}_lr_ready.csv')
