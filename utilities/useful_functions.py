#!python
# This file is subject to the terms and conditions defined in
# file 'LICENCE', which is part of this source code package.
# Author: Leo Guignard (leo.guignard...@AT@...gmail.com)


import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import os
import scanpy as sc
import anndata as ad

def plot_gene_distribution(adata, genes, c_name='Cluster {:s}',
                           bins=None, clustering='leiden',
                           fig_arr=None, figsize=None,
                           logscale=True, density=False):
    if bins is None:
        bins = 50
    if not isinstance(genes, list):
        genes = [genes]
    categories = adata.obs[clustering].cat.categories
    nb = len(categories)
    if fig_arr is None:
        rows = np.round(np.sqrt(nb)).astype(np.int)
        columns = np.ceil(nb/rows).astype(np.int)
    else:
        rows, columns = np.array(fig_arr).astype(np.int)
    if figsize is None:
        figsize = (4*columns, 2*rows+(rows//2))
    fig, axes = plt.subplots(rows, columns, figsize=figsize)
    axes = axes.flatten()
    ranges = []
    for g in genes:
        ranges += [(0, np.max(adata.raw[:, g].X))]
    ranges = np.max(ranges, axis=0)
    if len(genes)==1:
        histtype='bar'
    else:
        histtype='step'
    for i, index in enumerate(categories):
        tmp = adata[adata.obs[clustering]==index]
        values = []
        for g in genes:
            values += [tmp.raw[:, g].X.flatten()]
        ax = axes[i]
        ax.hist(values, bins=bins, density=density, range=ranges,
                log=logscale, histtype=histtype, label=genes)
        ax.set_title(c_name.format(index))
        if i==0:
            ax.legend()
    for ax in axes[i+1:]:
        ax.remove()
    fig.tight_layout()
    return fig

def cluster_and_plot(adata, val, c_name=None):
    if c_name is None:
        c_name = 'cluster.v{:.3f}'.format(val)
    sc.tl.leiden(adata, val)
    fig = sc.pl.umap(adata, color=['leiden', 'Bmp2', 'Sox9'], size=80, show=False, return_fig=True)
    fig.set_figwidth(20)
    fig.set_figheight(5)
    fig.savefig('figures/Umap.{:s}.pdf'.format(c_name))
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False,
                            return_fig=True, save='.t-test.{:s}.pdf'.format(c_name))
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False,
                            return_fig=True, save='.wilcoxon.{:s}.pdf'.format(c_name))

def is_number(s):
    """ Returns True is string is a number. """
    try:
        float(s)
        return True
    except ValueError:
        return False

def get_threshold(th, tab):
    while not is_number(th):
        th = input('Please enter a numeric value: ')
    th = float(th)
    filt_tab = tab < th
    rem = np.sum(filt_tab==False)
    len_tab = len(filt_tab)
    print('You are removing {:d} cells over a total of {:d} ({:.2f}%)'.format(rem, len_tab, 100*rem/len_tab))
    ans = input('Are you satisfied? (y/n) ')
    while ans!='y':
        th = input('Please enter a new threshold value: ')
        while not is_number(th):
            th = input('Please enter a numeric value: ')
        th = float(th)
        filt_tab = tab < th
        rem = np.sum(filt_tab==False)
        print('You are now removing {:d} cells over a total of {:d} ({:.2f}%)'.format(rem, len_tab, 100*rem/len_tab))
        ans = input('Are you satisfied? (y/n) ')
    return th

def get_clusters_et_al(path, size=5, filter_ncounts=False, filter_mito=False, reload_file=False):
    f = path
    ext = os.path.splitext(f)[-1]
    if not reload_file and (ext == '.h5ad' or os.path.exists(path.replace(ext, '.h5ad'))):
        results_file = path.replace(ext, '.h5ad')
        adata = ad.read_h5ad(results_file)
    elif ext == '.csv' or os.path.exists(path.replace(ext, '.csv')):
        results_file = path.replace(ext, '.h5ad')
        path = path.replace(ext, '.csv')
        adata = ad.read_csv(path).transpose()
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.filter_cells(adata, min_genes=200)
        if filter_ncounts:
            if isinstance(filter_ncounts, bool):
                adata.obs['n_counts'] = adata.X.sum(axis=1)
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
                ax1.hist(adata.obs.n_genes, bins=100,
                        range=(0, np.percentile(adata.obs.n_genes, 99)))
                ax1.set_xlabel('Number of counts')
                ax1.set_ylabel('Number of cells')
                ax2.hist(adata.obs.n_genes, bins=100, cumulative=True, density=True,
                        range=(0, np.percentile(adata.obs.n_genes, 99)))
                ax2.set_xlabel('Number of counts')
                ax2.set_ylabel('Ratio of cells')
                ax2.grid(True, axis='both')
                fig.tight_layout()
                plt.show()
                th_ncount = input('Please enter the threshold value for the maximum number of counts: ')
                th_ncount = get_threshold(th_ncount, adata.obs.n_genes)
                while not is_number(th_ncount):
                    th_ncount = input('Please enter a numeric value: ')
                th_ncount = float(th_ncount)
            else:
                th_ncount = filter_ncounts
            fig, ax = plt.subplots(1, 1)
            filter_tab_ncounts = adata.obs.n_genes<th_ncount
            ax.hist([adata.obs.n_genes[filter_tab_ncounts], adata.obs.n_genes[filter_tab_ncounts==False]],
                    color=['k', 'r'], label=['kept', 'removed'], bins=100, histtype='barstacked',
                    range=(0, np.percentile(adata.obs.n_genes, 99)))
            ax.set_xlabel('Number of counts')
            ax.set_ylabel('Number of cells')
            ax.legend()
            plt.show()
        else:
            filter_tab_ncounts = np.ones(adata.shape[0], dtype=bool)
        if filter_mito:
            if isinstance(filter_mito, bool):
                mito_genes = (adata.var_names.str.startswith('mt-') |
                              adata.var_names.str.startswith('Mt-') |
                              adata.var_names.str.startswith('MT-'))
                adata.obs['percent_mito'] = np.sum(
                    adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
                ax1.hist(adata.obs.percent_mito, bins=100,
                        range=(0, np.percentile(adata.obs.percent_mito, 99)))
                ax1.set_xlabel('Percent of mito expression')
                ax1.set_ylabel('Number of cells')
                ax2.hist(adata.obs.percent_mito, bins=100, cumulative=True, density=True,
                        range=(0, np.percentile(adata.obs.percent_mito, 99)))
                ax2.set_xlabel('Percent of mito expression')
                ax2.set_ylabel('Ratio of cells')
                ax2.grid(True, axis='both')
                fig.tight_layout()
                plt.show()
                th_mito = input('Please enter the threshold value for the maximum percent of mito expression: ')
                th_mito = get_threshold(th_mito, adata.obs.percent_mito)
                while not is_number(th_mito):
                    th_mito = input('Please enter a numeric value: ')
                th_mito = float(th_mito)
            else:
                th_mito = filter_mito
            plt.close(fig)
            filter_tab_mito = adata.obs.percent_mito<th_mito
            fig, ax = plt.subplots(1, 1)
            ax.hist([adata.obs.percent_mito[filter_tab_mito], adata.obs.percent_mito[filter_tab_mito==False]],
                    color=['k', 'r'], label=['kept', 'removed'], bins=100, histtype='barstacked',
                    range=(0, np.percentile(adata.obs.percent_mito, 99)))
            ax.set_xlabel('Number of counts')
            ax.set_ylabel('Number of cells')
            ax.legend()
            plt.show()
        else:
            filter_tab_mito = np.ones(adata.shape[0], dtype=bool)
        final_filt = np.ones(adata.shape[0], dtype=bool)
        if filter_ncounts:
            final_filt[filter_tab_ncounts==False] = False
        if filter_mito:
            final_filt[filter_tab_mito==False] = False
        both = np.sum((filter_tab_ncounts.astype(int)+filter_tab_mito)==0)
        diff = filter_tab_ncounts.astype(int)-filter_tab_mito
        nc = np.sum(diff==-1)
        mito = np.sum(diff==1)
        pie_values = [np.sum(final_filt), nc, both, mito]
        fig, ax = plt.subplots()
        ax.pie(pie_values, labels=['Kept', 'ncounts', 'both', 'mito'],
                shadow=False, startangle=90)
        ax.axis('equal')

        adata = adata[final_filt, :]
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata.raw = adata
        adata = adata[:, adata.var.highly_variable]
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, .25)
        adata.write(results_file)
    else:
        print('Can only work with .csv or .h5ad files (you gave {})'.format(path))
        return
    # sc.pl.highly_variable_genes(adata.raw)
    sc.pl.pca(adata, color=['Bmp2', 'Sox9', 'Sox17'], size=size)
    sc.pl.pca_variance_ratio(adata, log=True)

    fig = sc.pl.umap(adata, color=['Bmp2', 'Sox9', 'Wnt3'], size=size, show=False, return_fig=True)
    fig.set_figwidth(20)
    fig.set_figheight(6)

    c_name = os.path.splitext(os.path.split('data/'+results_file)[-1])[0]+'.0.25'
    fig = sc.pl.umap(adata, color=['leiden', 'Bmp2', 'Sox9'], size=size, show=False, return_fig=True)
    fig.set_figwidth(20)
    fig.set_figheight(5)
    fig.savefig('figures/Umap.{:s}.pdf'.format(c_name))

    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=True,
                            return_fig=True, save='.t-test.{:s}.pdf'.format(c_name))
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=True,
                            return_fig=True, save='.wilcoxon.{:s}.pdf'.format(c_name))
    return adata

def find_gene(adata, to_find):
    to_find_L = to_find.lower()
    found = []
    for gene_name in adata.raw.var_names:
        if to_find_L in gene_name.lower() or gene_name.lower() in to_find_L:
            found += [gene_name]
    print('{n:d} genes found corresponding to {tf:s}:'.format(n=len(found), tf=to_find))
    for f in found:
        print('\t{f:s}'.format(f=f))

def enseml2gene(fi, f_ensembl='Ensembl2Gene.txt'):
    if fi[-4:] != '.csv':
        print('Can only work with csv files (you gave: {})'.fromat(fi))
        return
    with open(f_ensembl) as f:
        lines = f.readlines()[1:]
    g2Ens = {}
    Ens2g = {}
    count = {}
    done = set()
    for l in lines:
        ens, g = l.split()
        if not ens in done:
            if g in count:
                count[g] += 1
            else:
                count[g] = 1

            gene = g if count[g]==1 else g+('='*(count[g]-1))
            g2Ens[gene]= ens
        done.add(ens)
    Ens2g = {}
    for k, v in g2Ens.items():
        Ens2g.setdefault(v, []).append(k)
    Ens2g = {k:sorted(v, key=len)[0] for k, v in Ens2g.items()}
    not_found = []
    with open(fi) as f:
        lines = f.readlines()
    new_fi = fi.replace('.csv', '.Genes.csv')
    new_lines = []
    with open(new_fi, 'w') as f:
        new_lines += [lines[0]]
        for l in lines[1:]:
            Ens = l.split(',')[0]
            new_l = l.replace(Ens, Ens2g.get(Ens, Ens))
            if not Ens in Ens2g:
                not_found += [Ens]
            new_lines += [new_l]
        f.writelines(new_lines)
        f.close()