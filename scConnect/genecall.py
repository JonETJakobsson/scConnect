import scConnect as cn
version = cn.database.version
# Module for different gene calling strategies
# Local functions


def get_adata_df(adata, transformation, layer, use_raw):
    """return a pandas DataFrame with gene expression counts and index and column lables.
    
    transformation: str or function, specify which transformation the data has been 
    subjected to so that count or UMI values can be calculated. By default, assumes that the data has been normalized using logp1.
    log1p back-transforms with expm1.
    Use a lambda call or provide a function that transforms the data back to the scale of gene counts. 
    It is fine to leave size factor normalization and linear scalings, but the data should not be logaritmised.
    layer: provide layer name to get specific layers.
    use_raw: weather to provide the stored raw data, or current adata.X matrix

    Returns: pandas DataFrame with values from specified layer
    """
    import pandas as pd
    from numpy import expm1

    # Get requested data
    if use_raw: # use raw if this is disired
        data = adata.raw.X.T
    elif layer != None: # else get the layer specified
        try:
            data = adata.layers[layer].T
        except: "layer not found"
    else: # if nothing is stated, provide current adata.X
        data = adata.X.T

    # transform data to counts
    if transformation == "log1p":
        data = expm1(data)

    if callable(transformation):
        print("using provided function")
        data = transformation(data)

    idx = adata.var_names
    col = adata.obs_names

    df = pd.DataFrame(data, index=idx, columns=col)
    return df


# Basic gene calling algorithms
def meanThreshold(adata, groupby, threshold, return_df=False, layer=None, use_raw=False, transformation="log1p"):
    """Binarize gene expression for groups aggregated by mean.

    Returns: adata object with updated uns.gene_call
    """
    from sklearn.preprocessing import binarize
    import pandas as pd
    df = get_adata_df(adata, layer=layer, use_raw=use_raw, transformation=transformation)
    result = df.groupby(by=adata.obs[groupby], axis=1).mean()
    binarize(result, threshold=threshold, copy=False)

    if return_df is True:
        return result
    else:
        adata.uns.update({"gene_call": result})
        return adata


def medianThreshold(adata, groupby, threshold, return_df=False, layer=None):
    """Binarize gene expression for groups aggregated by median.

    Returns: adata object with updated uns.gene_call
    """
    
    from sklearn.preprocessing import binarize
    import pandas as pd
    df = get_adata_df(adata, layer)
    result = df.groupby(by=adata.obs[groupby], axis=1).median()
    binarize(result, threshold=threshold, copy=False)

    if return_df is True:
        return result
    else:
        adata.uns.update({"gene_call": result})
        return adata


def meanExpression(adata, groupby, return_df=False, layer=None, use_raw=False, normalization=False, transformation="log1p"):
    """return gene expression for groups aggregated by mean.

    .. warning::
        Not compatible with certain molecular ligands as values are not between 0 and 1.
        Exclusion genes are weighted as:
        
        .. math:: Ligand or receptor call = Geom(Synthesis, Transporter, Reuptake) * (1 - Exclusion)
        
        Hence gene expression levels close to 1 contribute to 
        exlude the ligand and expression levels close to 0 leaves
        the ligand call untouched. 

    use_raw: Boolian, True if using the raw data of the adata structure. else False

    Normalization: Boolian, Normalize the gene expression over the different groups for each gene using L2 normalization.
    Note that this will lead to highly inflated values for lowly expressed genes, and under estimation of highly expressed genes.
    For experimental use only.

    Returns: adata object with updated uns.gene_call
    """

    import pandas as pd
    import numpy as np
    from sklearn.preprocessing import normalize

    df = get_adata_df(adata=adata, layer=layer, use_raw=use_raw, transformation=transformation)

    result = df.groupby(by=adata.obs[groupby], axis=1).mean()

    if normalization:
        df_norm = normalize(result, copy=True, axis=1, norm="l2")
        result = pd.DataFrame(df_norm, columns=result.columns, index=result.index)

    if return_df is True:
        return result
    else:
        adata.uns["gene_call"] = result.to_dict()
        return adata


def percentExpression(adata, groupby, return_df=False, layer=None, use_raw=False,):
    """Calculate percent of expressing cells in a cluster
    
    Returns: adata object with updated uns.gene_call
    """
    import pandas as pd
    from sklearn.preprocessing import binarize
    df = get_adata_df(adata, transformation=None, layer=layer, use_raw=use_raw)
    binarize(df, threshold=0.0, copy=False)
    result = df.groupby(by=adata.obs[groupby], axis=1).mean()

    if return_df is True:
        return result
    else:
        adata.uns.update({"gene_call": result})
        return adata


# Implementing betabinomial trinarization as seen in Zeisel et al. 2018, Cell

def betabinomialTrinarization(adata, groupby, pep=0.05, f=0.2, return_df=False, verbouse=True, layer=None):
    """Trinarize gene expression for groups.
    
    Function is adapted from `Zeisel 2018`__.

    __ https://www.sciencedirect.com/science/article/pii/S009286741830789X


    pep: posterior error probability, float (0-1)
    
    f: fraction, float (0-1)

    Returns: adata object with updated uns.gene_call
    """
    import pandas as pd
    from collections import Counter
    from scConnect.tools import printProgressBar

    df = get_adata_df(adata, layer)
    labels = adata.obs[groupby]
    clusters = dict(Counter(labels)).keys()
    gene_call = pd.DataFrame(columns=clusters)
    tot_genes = df.shape[0]
    called_gene = 0
    for gene, expression in df.iterrows():
        expression_by_label = betabinomial_trinarize_array(
            expression, labels, pep, f)
        gene_call.loc[gene] = expression_by_label
        called_gene += 1
        if verbouse:
            printProgressBar(called_gene, tot_genes,
                             prefix=f"Calling gene expression using betabinomial trinarization [Zeisel et al. 2018] gene: {called_gene} of {tot_genes}")
    if return_df is True:
        return gene_call
    else:
        adata.uns.update({"gene_call": gene_call})
        return adata


def betabinomial_trinarize_array(array, labels, pep, f, n_labels=None):
    import numpy as np
    from collections import Counter

    if n_labels is None:
        #n_labels = np.max(labels) + 1
        n_labels = len(set(labels))

    n_by_label = dict(Counter(labels))
    k_by_label = dict()

    for lbl in n_by_label.keys():
        try:
            k_by_label[lbl] = np.count_nonzero(
                array[np.where(labels == lbl)[0]])
        except KeyError:
            print(f"no cells from cluster {lbl}")
            k_by_label[lbl] = 0

    ns = [n for n in n_by_label.values()]
    ks = [k for k in k_by_label.values()]

    vfunc = np.vectorize(p_half)
    ps = vfunc(ks, ns, f)

    expr_by_label = np.zeros(n_labels) + 0.5
    expr_by_label[np.where(ps > (1 - pep))[0]] = 1
    expr_by_label[np.where(ps < pep)[0]] = 0

    return (expr_by_label)


def p_half(k, n, f):
    from math import exp, lgamma, log
    from scipy.special import beta, betainc, betaln
    import numpy as np

    # These are the prior hyperparameters beta(a,b)
    a = 1.5
    b = 2

    incb = betainc(a + k, b - k + n, f)

    if incb == 0:
        p = 1.0
    else:
        p = 1.0 - exp(log(incb) + betaln(a + k, b - k + n) +
                      lgamma(a + b + n) - lgamma(a + k) - lgamma(b - k + n))

    return p


#Calculate new mean values based on Zero inflation model
def ziMean(adata, groupby, organism="mmusculus"):
    """Calculates mean expression based on estimated dropouts given the mean expression of non zero cells.
    
    Using the decay coeficient calculated by the ZIFA algorithm, estimated dropput rates can 
    be calculated based on the mean expression of non zero counts. New mean expression values is then calculated
    from imputed counts.
    """

    import ZIFA.block_ZIFA as zf
    import scanpy as sc
    adata_copy = adata.copy()
    adata_copy = filter_genes(adata_copy, organism=organism)
    sc.pp.filter_genes(adata_copy, min_counts=1)
    countmatrix = adata_copy.X
    model, params = zf.fitModel(countmatrix, 2, singleSigma=True)
    dc = params["decay_coef"]
    df = get_adata_df(adata_copy)

    def ZImean(countmatrix):
        """Calculated new mean values given set decay coeficiant. Can be used with groupby.aggregate"""

        import numpy as np
        gene_call = dict()
        for gene, expression in countmatrix.iterrows():
            total_cells = len(expression)
            expressing_cells = len(expression[expression>0])
            if expressing_cells>0:
                non_zero_mean = np.mean(expression[expression>0])
                prob = np.exp(-dc*non_zero_mean**2)
                with_dropout = int(np.round(expressing_cells/(1-prob), decimals=0))
                if with_dropout > total_cells:
                    with_dropout = total_cells

                new_counts = [non_zero_mean]*with_dropout+[0]*(total_cells-with_dropout)
                new_mean = np.mean(new_counts)
            else:
                new_mean = 0
            gene_call[gene] = new_mean
    
        return gene_call

    
    gene_call = df.groupby(adata_copy.obs[groupby], axis=1).aggregate(ZImean)

    adata.uns.update({"gene_call": gene_call})
    return adata
# Plotting tools to investigate gene calling performance.


def geneParticipation(adata):
    """Plot histogram showing % of cluster genes are prescent in."""
    import numpy as np
    df = adata.uns["gene_call"]
    cluster_participation = df.aggregate(np.mean, axis=1)
    cluster_participation.plot(kind="hist", bins=20)


def filter_genes(adata, organism="mmusculus"):
    """Filter adata object to only have ligand and receptor genes.
    
    .. note ::
        Usefull to create ligand and receptor centric clustering and to speed up receptor and ligand calling.
    """
    import numpy as np
    import pandas as pd
    import pkg_resources

    ligands = pd.read_csv(pkg_resources.resource_filename(
        __name__, (f"data/Gene_annotation/{version}/{organism}/ligands.csv")))
    receptors = pd.read_csv(pkg_resources.resource_filename(
        __name__, (f"data/Gene_annotation/{version}/{organism}/receptors.csv")))

    def add_gene(gene, gene_list):
        try:
            gene = eval(str(gene))
        except (NameError, TypeError):
            pass
        if type(gene) == str:
            gene_list.append(gene)
            return gene_list
        elif type(gene) == list:
            for g in gene:
                gene_list.append(g)
            return gene_list
        else:
            return gene_list

    gene_list = list()
    ligands = ligands[["preprogene", "synthesis",
                       "transport", "reuptake", "excluded"]]
    for column in ligands.columns:
        for gene in ligands[column]:
            gene_list = add_gene(gene, gene_list)

    for gene in receptors["gene"]:
        gene_list = add_gene(gene, gene_list)

    gene_list = list(set(gene_list))
    filter_genes = [
        True if gene in gene_list else False for gene in adata.var_names]
    adata = adata[:, filter_genes]
    return adata

# Finding optimal settings using Glutamate and GABA segregation


def segregation_score(adata):
    """Calculate a score of segregation of gabaergic and glutamatergic cells"""
    import pandas as pd
    total_score = []
    df = pd.DataFrame(adata.uns["ligands"]).loc[["GABA", "L-glutamic acid"]]
    for cluster in df:
        gaba = df[cluster].loc["GABA"] > 0
        glut = df[cluster].loc["L-glutamic acid"] > 0

        if gaba and glut:  # cell population expresses both ligands
            total_score.append(-1)
        elif not gaba and not glut:  # cell population expresses no ligands
            total_score.append(-1)
        else:  # cell population expresses only one ligand
            total_score.append(1)

    # normalize on total number of cell populations
    return(sum(total_score)/len(total_score))


def plot_segregation(adata, save=False, filename=None):
    """Plot gabaergic and glutamaterig cell populations"""
    import holoviews as hv
    from holoviews import opts
    import pandas as pd
    hv.extension("matplotlib")

    df = pd.DataFrame(adata.uns["ligands"]).loc[[
        "GABA", "L-glutamic acid"]].stack().reset_index()
    df.columns = ["ligand", "cluster", "value"]
    df = df.sort_values(by="cluster", axis=0)

    opts.defaults(
        opts.Bars(
            stacked=True,
            xrotation=90,
            legend_position="right",
            ylabel="Ligand score"
        )
    )
    bars = hv.Bars(df, kdims=["cluster", "ligand"])
    if save is True:
        hv.save(bars, filename)
    return bars


def optimize_segregation(adata, groupby, method="meanTH", start=0, stop=1, steps=10, organism="mmusculus", fraction=0.8, iterations=10, save=False, filename=None, return_df=False, verbouse=False, transformation="log1p", use_raw=False):
    """Calculate segregation score for a range of thresholds using specified method.
    
    available methods:
    meanTH
    medianTH
    NBT

    The estimation will be iterated over selected number of times, and an average segregation will 
    be plotted. error area is 99% CI.

    Returns: Seaborn lineplot of score over threshold

    Save file using save=True and a file name with .png or.svg extention
    """
    from collections import OrderedDict
    import scanpy as sc
    from random import randint, sample
    import seaborn as sns
    import pandas as pd
    import numpy as np
    import scConnect as cn
    from scConnect.tools import printProgressBar

    # filter out only genes nessesary to calculate gaba and glutamate presence
    mol_ligands = cn.database.get_molecule_ligands(organism=organism)
    seg_filter = [True if ligand in ["GABA", "L-glutamic acid"] else False for ligand in mol_ligands["ligand"]]
    mol_ligands = mol_ligands[seg_filter][["synthesis","transport","reuptake"]]
    gene_list = list()
    for genes in mol_ligands.values.flatten():
        if genes is not None:
            for gene in genes:
                    gene_list.append(gene)

    adata = adata[:, gene_list]
    
    th_score = list()
    space = np.linspace(start, stop, steps)
    tot_steps = len(space)*iterations
    cur_step = 1
    obs = adata.obs_names

    for th in space:
        for iteration in np.arange(0,iterations):

            # Create random subset of dataset for distribution calculation
            obs_s = sample(set(obs), int(len(obs)*fraction))
            temp = adata[obs_s]

            #Store method calls with relevant parameters
            methods = {
                "NBT": "betabinomialTrinarization(temp, groupby=groupby, f=th, verbouse=False)",
                "medianTH": "medianThreshold(temp, groupby=groupby, threshold=th, layer=None, use_raw=use_raw, transformation=transformation)",
                "meanTH": "meanThreshold(temp, groupby=groupby, threshold=th, layer=None, use_raw=use_raw, transformation=transformation)"
                }

            #Runs selected method with relevant local variabls
            temp = eval(methods[method])
            temp = cn.connect.ligands(temp, organism=organism, select_ligands=[
                                    "GABA", "L-glutamic acid"])
            score = segregation_score(temp)
            th_score.append((iteration, th, score))
            if verbouse:
                printProgressBar(cur_step, tot_steps,
                                prefix=f"Finding optimized threshold")
            cur_step += 1

    # store data in a pandas dataframe
    df = pd.DataFrame(data=th_score, columns=["Iteration", "Threshold", "Segregation"])

    # group by threshold to calculate best score and best threshold
    g_df = df.groupby(by="Threshold")["Segregation"].mean()

    # idxmax will return the first index with max score, which is what we want.
    max_segregation = g_df[g_df.idxmax()]
    best_th = g_df.idxmax()
    print(f"Maximun segregation of {max_segregation} was achived at threshold {best_th}")


    if return_df:
        return df

    f = sns.lineplot(x=df["Threshold"], y=df["Segregation"], ci=99, estimator="mean")
    f.axvline(best_th, ls="--", color="red")

    if save is True:
        fs = f.get_figure()
        fs.savefig(filename)
    return f