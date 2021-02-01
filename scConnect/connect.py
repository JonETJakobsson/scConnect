# Module deals with creation of ligand and receptor scores, and creation of scConnect tables etc.
import scConnect as cn
import scanpy as sc

version = cn.database.version
organism = cn.database.organism
# Scoring logic for ligands


def ligandScore(ligand, genes):
    """calculate ligand score for given ligand and gene set"""

    from scipy.stats.mstats import gmean
    import numpy as np

    if ligand.ligand_type == "peptide" and isinstance(ligand.preprogene, str):
        # check if multiple genes needs to be accounted for
        if isinstance(eval(ligand.preprogene), list):
            ligand_genes = list()
            for gene in eval(ligand.preprogene):
                try:
                    ligand_genes.append(genes[gene])
                except KeyError:
                    #print(f"{gene} not found")
                    ligand_genes.append(0.0)
            # use max, as there might be many orthologs genes for one original 
            # gene and not all have to be expressed
            try:
                ligand_score = max(ligand_genes)
            except ValueError:
                print(f"something is wrong with the list {ligand_genes}")
                ligand_score = 0.0

            return ligand_score


    elif ligand.ligand_type == "molecule":
        synthesis = ligand.synthesis
        transport = ligand.transport
        reuptake = ligand.reuptake
        excluded = ligand.excluded

        # get geometric mean of synthesis genes (all need to be present)
        if not isinstance(synthesis, str):
            # If no genes are needed, synthesis is set to nan
            synthesis = np.nan
        else:
            synthesis_expression = list()
            for gene in eval(synthesis):
                try:
                    synthesis_expression.append(genes[gene])
                except KeyError:
                    # If gene was not found append 0
                    #print(f"{gene} not found")
                    synthesis_expression.append(0.0)
            synthesis = gmean(synthesis_expression)

        # get maximum of vesicle transporters (only one is needed for molecule transport)
        if not isinstance(transport, str):
            # If no specific genes are needed, set transport to nan
            transport = np.nan
        else:
            transport_expression = list()
            for gene in eval(transport):
                try:
                    transport_expression.append(genes[gene])
                except KeyError:
                    # If gene was not found append 0
                    #print(f"{gene} not found")
                    transport_expression.append(0.0)
            transport = max(transport_expression)

        # Get maximum of reuptake genes (only one is needed)
        if not isinstance(reuptake, str):
            # If no specific genes are needed, set reuptake to nan
            reuptake = np.nan
        else:
            reuptake_expression = list()
            for gene in eval(reuptake):
                try:
                    reuptake_expression.append(genes[gene])
                except KeyError:
                    # If gene was not found append 0
                    #print(f"{gene} not found")
                    reuptake_expression.append(0.0)
            reuptake = max(reuptake_expression)

        # get maximum among exluding genes where any gene expression divert to other ligands
        if not isinstance(excluded, str):
            # If no specific genes are needed, set excluded to 0
            excluded = 0
        else:
            excluded_expression = list()
            for gene in eval(excluded):
                try:
                    excluded_expression.append(genes[gene])
                except KeyError:
                    # If gene was not found append 0
                    #print(f"{gene} not found")
                    excluded_expression.append(0.0)
            excluded = max(excluded_expression)

        # return geometric mean of synthesis, transport and reuptake multipled exclusion
        promoting_factor = gmean(([x for x in [synthesis, transport, reuptake] if str(x) != "nan"])) # genes driving ligand production, remove nan values
        
        if str(promoting_factor) == "nan": # capture cases where no promoting genes were present
            print(f"no promoting genes detected for {ligand.ligand}")
            return 0.0  # exit before running exclusion calculation
            
        ligand_score = promoting_factor - excluded # correct ligand expression based on the exclusion factor
        if ligand_score < 0: # ligand score should be 0 or positive
                ligand_score = 0.0
        return ligand_score

    # If genes are missing from ligand gene list
    else:
        print("Big error! ligand type is not defined!")
        return 0.0


def ligands(adata, organism=organism, select_ligands=None):
    """return a dataframe with ligand scores for each cluster.

    .. note::
        Needs a gene call dataframe under adata.uns.gene_call. 
        Use scConnect.genecall to create such dataframe

    organism defaults to mouse, to use genes for other organism select this here.

    use select_ligands to only asses given ligands 
    (used by optimize_segregation to only check for gaba and glutamate)

    Returns: Dict of ligand call for each cluster.
    """
    import scConnect as cn
    import pkg_resources
    import pandas as pd

    ligands = pd.read_csv(pkg_resources.resource_filename(
        __name__, (f"data/Gene_annotation/{version}/{organism}/ligands.csv")))
    if isinstance(select_ligands, list):
        select = [True if ligand in select_ligands else False for ligand in ligands.ligand]
        ligands = ligands[select]
    ligand_df = pd.DataFrame(index=ligands.ligand)
    
    for cluster, genes in adata.uns["gene_call"].items():
        cluster_scores = list()
        for ligand_data in ligands.iterrows():
            ligand = ligand_data[1]
            # fetch ligand score for specific ligand and gene set
            ligand_score = ligandScore(ligand, genes)
            cluster_scores.append(ligand_score)
        ligand_df[cluster] = cluster_scores

    adata.uns["ligands"] = ligand_df.to_dict()
    return adata


# Scoring logic for receptors
def receptorScore(receptor, genes):
    """calculate receptor score given receptor and gene set"""

    from scipy.stats.mstats import gmean

    gene_expression = list()
    for gene in eval(receptor.gene):
        try:
            gene_expression.append(genes[gene])
        except KeyError:
            # If gene was not found append 0
            #print(f"{gene} not found")
            gene_expression.append(0.0)
    # use max, as several genes might be found during ortholog search, 
    # not all might bee needed to create the receptor
    gene_expression = max(gene_expression)
    return gene_expression


def receptors(adata, organism=organism):
    """return a dataframe with receptor scores for each cluster.

    .. note::
        Needs a gene call dataframe under adata.uns.gene_call. 
        Use scConnect.genecall to create such dataframe.

    Returns: Dict of receptor call for each cluster.
    """
    import scConnect as cn
    import pkg_resources
    import pandas as pd

    receptors = pd.read_csv(pkg_resources.resource_filename(
        __name__, (f"data/Gene_annotation/{version}/{organism}/receptors.csv")))
    receptor_df = pd.DataFrame(index=receptors.receptor)

    for cluster, genes in adata.uns["gene_call"].items():
        cluster_scores = list()
        for receptor_data in receptors.iterrows():
            receptor = receptor_data[1]
            # fetch ligand score for specific ligand and gene set
            receptor_score = receptorScore(receptor, genes)
            cluster_scores.append(receptor_score)
        receptor_df[cluster] = cluster_scores
    adata.uns["receptors"] = receptor_df.to_dict()
    return adata


# Interaction logic

def interactions(emitter, target, self_reference=True, organism=organism, corr_pval=True):
    """return an edge list of interactions between clusters.
    If all connections are of interest, use the same data source for
    emitter and target.

    .. note::
        self_reference is only valid when emitter == target.

    .. note::
        edge_list is returned as a list, and not in a adata object. 
        This is since multiple adata objects can be passed in to the
        function, and whould lead to ambiguity of which object to append the edge_list to.
    
    Returns: List of edges between given emmitor and target clusters.
    """
    import pkg_resources
    import pandas as pd
    from itertools import product
    from scConnect.tools import printProgressBar

    interactions = pd.read_csv(pkg_resources.resource_filename(
        __name__, (f"data/Gene_annotation/{version}/interactions.csv")), index_col=[0, 1], sep=";")

    interactions.sort_index(axis="index", inplace=True)
    # Create a set of all possible index combinations.
    # This is used to test if ligand receptor combination is present.
    interaction_set = set(interactions.index)

    # An edge list should contain u, v and d,
    # where u is input node, v is output node
    # and d is a dictionary with edge attributes.
    edge_list = list()

    # get all clusters
    # NOTE: if the same cluster name is used in emitter and target datasets, they are
    # assumed to be the same cluster. Give your clusters uniqe names between your datasets.
    try:
        emitter_clusters = pd.DataFrame(emitter.uns["ligands"]).columns
        target_clusters = pd.DataFrame(target.uns["ligands"]).columns
    except KeyError:
        print(
            f"Please run connect.ligands() and connect.receptors() on your datasets first")
        return

    # Calculate total number of cluster combinations for the progress bar
    if self_reference is True:
        total_comb = len(list(product(emitter_clusters, target_clusters)))
    else:
        total_comb = len([(e, t) for (e, t) in product(
            emitter_clusters, target_clusters) if e != t])

    ligands = pd.DataFrame(emitter.uns["ligands"])
    receptors = pd.DataFrame(target.uns["receptors"])

    # load extra ligand and receptor statistics
    ligands_zscore = pd.DataFrame(emitter.uns["ligands_zscore"])
    receptors_zscore = pd.DataFrame(target.uns["receptors_zscore"])
    if corr_pval:
        ligands_pval = pd.DataFrame(emitter.uns["ligands_corr_pval"])
        receptors_pval = pd.DataFrame(target.uns["receptors_corr_pval"])
    else:
        ligands_pval = pd.DataFrame(emitter.uns["ligands_pval"])
        receptors_pval = pd.DataFrame(target.uns["receptors_pval"])

    # Fetch receptor and ligand information
    receptor_info = pd.read_csv(pkg_resources.resource_filename(
        __name__, (f"data/Gene_annotation/{version}/{organism}/receptors.csv")), index_col=1)
    receptor_info = receptor_info[["family", "gene"]]

    ligand_info = pd.read_csv(pkg_resources.resource_filename(
        __name__, (f"data/Gene_annotation/{version}/{organism}/ligands.csv")), index_col=1)
    ligand_info = ligand_info[["ligand_type", "comment"]]

    # Nested for-loop to get all combinations of
    # interactions between clusters.
    comb_tried = 0
    for emitter_cluster in emitter_clusters:
        for target_cluster in target_clusters:
            # Are we interested in self referencing information?
            # I leave that up to the user
            if emitter_cluster != target_cluster or self_reference == True:
                # Get only ligands and receptors expressed by the clusters
                # (speeds up itterative functions later)
                emitter_ligands = ligands[emitter_cluster][ligands[emitter_cluster] > 0]
                target_receptors = receptors[target_cluster][receptors[target_cluster] > 0]

                connections = get_connections(
                    emitter_ligands, 
                    target_receptors, 
                    interactions, 
                    interaction_set, 
                    receptor_info, 
                    ligand_info,
                    emitter_cluster,
                    target_cluster,
                    ligands_zscore,
                    ligands_pval,
                    receptors_zscore,
                    receptors_pval)

                if len(connections) > 0:
                    for connection in connections:
                        edge_list.append(connection)

                # Add the progress bar
                comb_tried += 1
                printProgressBar(
                    comb_tried, total_comb, prefix=f"finding connections between {len(emitter_clusters)} emitter clusters and {len(target_clusters)} target clusters")

    return edge_list


# get all connections based on Ligands and receptors, and provide score for interactions
# Also provide meta data as a dictionary for interaction
def scale(value, from_range=(0, 1), to_range=(10E-100, 1)): # mitagate log with 0
        value = to_range[0] + (to_range[1] - to_range[0]) * (value -from_range[0]) / (to_range[1] - to_range[0])
        return value

def get_connections(
    ligands,
    receptors, 
    interactions, 
    interaction_set, 
    receptor_info, 
    ligand_info,
    emitter_cluster,
    target_cluster, 
    ligands_zscore,
    ligands_pval,
    receptors_zscore,
    receptors_pval):
    """finds connections between ligands and receptors 
    and return a score and metadata for each interaction"""

    from scipy.stats.mstats import gmean
    import numpy as np
    # shorten the list of interactions to only contain relevant ligands.
    # This should speed up the algorithm
    ligand_filter = [True if ligand in ligands.keys() else False for ligand in interactions.index.get_level_values(0)]
    interactions = interactions.loc[ligand_filter]

    

    def interaction_specificity(l, r): # used to calculate interaction specificity score
        sig = -np.log10((l+r)/2)
        return sig

    connections = list()
    for ligand, l_score in ligands.iteritems():
        for receptor, r_score in receptors.iteritems():
            if (ligand, receptor) in interaction_set:
                interaction = interactions.loc[ligand, receptor]
                score = float(gmean((l_score, r_score)))
                ligand_pval = float(ligands_pval[emitter_cluster][ligand])
                receptor_pval = float(receptors_pval[target_cluster][receptor])
                specificity = float(interaction_specificity(ligand_pval, receptor_pval))
                log_score = float(np.log10(score + 1))
                importance = specificity * log_score

                connections.append((ligands.name, receptors.name, {
                    "score": float(score),
                    "log_score": log_score, # From here on, all values are +1ed and logaritmized with base of 10. # From here on, all values are +1ed and logaritmized with base of 10.
                    "ligand": str(ligand),
                    "ligand_zscore": float(ligands_zscore[emitter_cluster][ligand]),
                    "ligand_pval": ligand_pval,
                    "receptor": str(receptor),
                    "receptor_zscore": float(receptors_zscore[target_cluster][receptor]),
                    "receptor_pval": receptor_pval,
                    "interaction": f"{ligand} --> {receptor}",
                    "specificity": specificity,
                    "importance": importance,
                    "endogenous": f"{list(interaction.endogenous)}",
                    "action": f"{list(interaction.action)}",
                    "ligandspecies": f"{list(interaction.ligand_species)}",
                    "receptorspecies": f"{list(interaction.target_species)}",
                    "pubmedid": f"{list(interaction.pubmed_id)[:5]}",
                    "receptorfamily": str(receptor_info.loc[receptor]["family"]),
                    "receptorgene": str(receptor_info.loc[receptor]["gene"]),
                    "ligandtype": str(ligand_info.loc[ligand]["ligand_type"]),
                    "ligandcomment": str(ligand_info.loc[ligand]["comment"])}))
    return connections


def nodes(adatas):
    """
    Returns an list of nodes, attributes dictionary tuples. 
    Each tuple represent one node with an attribute dictionary: 
    *(cluster, dict(receptors: dict(receptor:score), ligands: dict(ligand:score) ))*
    
    """
    if not isinstance(adatas, list):
        adatas = [adatas, ]

    nodes = []

    for i, adata in enumerate(adatas):
        print(f"precessing adata #{i+1}")
        ligands_score = adata.uns["ligands"]
        ligands_zscore = adata.uns["ligands_zscore"]
        ligands_pval = adata.uns["ligands_pval"]
        ligands_corr_pval = adata.uns["ligands_corr_pval"]
        receptors_score = adata.uns["receptors"]
        receptors_zscore = adata.uns["receptors_zscore"]
        receptors_pval = adata.uns["receptors_pval"]
        receptors_corr_pval = adata.uns["receptors_corr_pval"]
        genes = adata.uns["gene_call"]
        clusters = ligands_score.keys()

        # Filter out ligands with positive score (remove non expressing ligands and receptors)
        for cluster in clusters:
            print(f"processing cluster {cluster}")
            cluster_ligands_score = {k: v for k,
                               v in ligands_score[cluster].items() if v > 0}
            cluster_receptors_score = {k: v for k,
                                 v in receptors_score[cluster].items() if v > 0}

            # Add all information to the node dicionary
            node = (cluster, {
                "ligands_score": cluster_ligands_score,
                "ligands_zscore": ligands_zscore[cluster],
                "ligands_pval": ligands_pval[cluster],
                "ligands_corr_pval": ligands_corr_pval[cluster],
                "receptors_score": cluster_receptors_score,
                "receptors_zscore": receptors_zscore[cluster],
                "receptors_pval": receptors_pval[cluster],
                "receptors_corr_pval": receptors_corr_pval[cluster],
                "genes": genes[cluster]})
            nodes.append(node)

    return nodes

# Statistic inference of ligand and receptor scores
# Here we shuffel the group annotations many times, calculate ligand and receptor scores
# and find the mean and standard deviation for each ligand/receptor score for each gorup.
# We can then calculate the z-score of the true ligand/receptor score, p-values and corrected p-values
# Data an be used to detect group specific expression of ligands and receptors.

def _ligand_receptor_call(adata, groupby, organism, transformation, return_df = True):
    import pandas as pd
    adata = cn.genecall.meanExpression(adata, groupby=groupby, normalization=False, use_raw=False, transformation=transformation)
    adata = cn.connect.ligands(adata, organism=organism)
    adata = cn.connect.receptors(adata, organism=organism)
    
    ligands = pd.DataFrame(adata.uns["ligands"])
    receptors = pd.DataFrame(adata.uns["receptors"])
    if return_df:
        return ligands, receptors

def _values_df(dfs):
    values_df = dfs[0].copy()
    
    for i in range(values_df.shape[0]):
        for j in range(values_df.shape[1]):
            values = list()
            for df in range(len(dfs)):
                values.append(dfs[df].iloc[i,j])
            values_df.iloc[i,j] = str(values)
    return values_df

            
def _mean_df(df):
    import numpy as np
    mean_df = df.copy()
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            mean_df.iloc[i,j] = np.mean(eval(df.iloc[i,j]))
    return mean_df
        
def _std_df(df):
    import numpy as np
    std_df = df.copy()
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            std_df.iloc[i,j] = np.std(eval(df.iloc[i,j]))
    return std_df
        

def _score_pv_df(mean, std, value, emperical, values, merge_dist):
    """Calculate z-scores and p-values for ligand and receptor 
    calls compared to random group designation
    
    returns: score_df and pval_df"""
    import numpy as np
    from scipy import stats

    assert (mean.shape == std.shape == value.shape), "dataframes are not of the same size, rerun ligand and receptor call"

    score_df = mean.copy()
    pval_df = mean.copy()
    warning = False # warning flag for if mean or std is 0 (mening no values were ever sampled to that group)
    faults = 0

    for i in range(score_df.shape[0]): # for each ligand and receptor
        if merge_dist == True:
            dist = list()
            for j in range(score_df.shape[1]):
                for val in eval(values.iloc[i,j]):
                    dist.append(val)
                
        for j in range(score_df.shape[1]): # for each celltype
            v = value.iloc[i,j]
            s = std.iloc[i,j]
            m = mean.iloc[i,j]
            if merge_dist == False:
                dist = eval(values.iloc[i,j])
        
            if s == 0: # sampeling never managed to include this ligand or receptor for this group
                z_score = 0.0
                pval = 1
                warning = True
                faults += 1
            else:
                z_score = (v-m)/s
                #pval =  float(stats.norm.sf(abs(z_score))*2) # Two tailed p-value
                pval =  float(stats.norm.sf(z_score)) # one tailed p-value

            if emperical == True: 
                # Calculate exact permutation pvalues emperically from the collected distribution
                # method from https://www.degruyter.com/view/journals/sagmb/9/1/article-sagmb.2010.9.1.1585.xml.xml
                # permutation without replacement (use full sequence)
                b = sum(dist > v)
                m = len(dist)
                pval = (b+1)/(m+1)

            
            score_df.iloc[i,j] = z_score
            pval_df.iloc[i,j] = pval
            
    if warning:
        total = score_df.shape[0] * score_df.shape[1]
        print(f"{faults/total*100} % of group metrices were 0. increase n to reduce this number")
    
    return score_df, pval_df

def _corrected_pvalue(pvalues, method="fdr_bh", scale_pval=False):
    """correct a dataframe of p-values to a dataframe of corrected p-values.
    
    Supports many different methods:
    bonferroni : one-step correction
    sidak : one-step correctio
    holm-sidak : step down method using Sidak adjustment
    holm : step-down method using Bonferroni adjustment
    simes-hochberg : step-up method (independent
    hommel : closed method based on Simes tests (non-negative)
    fdr_bh : Benjamini/Hochberg (non-negative)
    fdr_by : Benjamini/Yekutieli (negative)
    fdr_tsbh : two stage fdr correction (non-negative)
    fdr_tsbky : two stage fdr correction (non-negative)
    
    defaults to fdr_bh
    
    returns a pandas dataframe
    """
    import statsmodels.stats.multitest as mt
    import pandas as pd
    
    p_flat = pvalues.to_numpy().flatten()
    corr_p = mt.multipletests(p_flat, method=method)[1]
    
    corr_p = corr_p.reshape(pvalues.shape)
    corr_pval = pd.DataFrame(corr_p, columns=pvalues.columns, index=pvalues.index)
    
    # scale p values to remove abloslute 0 calls
    if scale_pval:
        corr_pval = scale(corr_pval)

    return corr_pval
    
def specificity(adata, n, groupby, organism=organism, return_values=False, transformation="log1p", emperical=True, merge_dist=False):
    """calculate statistics for the ligands and receptor scores.
    
    Compare the group ligand and receptor scores to the mean score of 
    that group after n number of permutations
    
    if emperical is True (default), calculates p-values emperically given the collected random distribution.
    p = (b+1)/(m+1) 
    where b is the number of permutated values higher than the observed
    and m is the number of permutations used (set this by the argument n)"""

    from random import shuffle
    import pandas as pd
    from scConnect.tools import printProgressBar
    _adata = adata.copy()
    groups = list(_adata.obs[groupby])
    
    ligand_dfs = list()
    receptor_dfs = list()
    
    # variable to store the setting in, can be used when saving the specificity
    settings = dict(
        groupby = groupby,
        permutations =  n,
        transformation = transformation,
        emperical = emperical,
        merge_dist = merge_dist
    )

    # Run normal ligand and receptor call without shuffel on original adata
    _ligand_receptor_call(adata, groupby=groupby, organism=organism, transformation=transformation, return_df=False)
    # shuffel group annotations n times and fetch ligand and receptor dataframes
    for i in range(n):
        printProgressBar(i+1, n, prefix=f"Shuffeling dataframe {i+1} out of {n}")
        shuffle(groups)
        _adata.obs[groupby] = groups
        ligand, receptor = _ligand_receptor_call(_adata, groupby=groupby, organism=organism, transformation=transformation)
        ligand_dfs.append(ligand)
        receptor_dfs.append(receptor)
    
    # Merge all dataframes to one datafram (with list of values for each element)
    ligand_values = _values_df(ligand_dfs)
    receptor_values = _values_df(receptor_dfs)
    
    # Calculate the mean values of the list in each element
    print("Calculating means...")
    ligand_mean = _mean_df(ligand_values)
    receptor_mean = _mean_df(receptor_values)
    
    # Calculate the standard deviation of the list in each element
    print("Calculating standard deviations...")
    ligand_std = _std_df(ligand_values)
    receptor_std = _std_df(receptor_values)
    
    # Calculate Z-scores, p-values and corrected p-values
    print("Calculating Z-score, p-values and corrected p-values...")
    ligand_value = pd.DataFrame(adata.uns["ligands"])
    ligand_score , ligand_pval = _score_pv_df(ligand_mean, ligand_std, ligand_value, emperical, ligand_values, merge_dist=merge_dist)
    ligand_corr_pval = _corrected_pvalue(ligand_pval, scale_pval=not emperical)
    
    receptor_value = pd.DataFrame(adata.uns["receptors"])
    receptor_score , receptor_pval = _score_pv_df(receptor_mean, receptor_std, receptor_value, emperical, receptor_values, merge_dist=merge_dist)
    receptor_corr_pval = _corrected_pvalue(receptor_pval, scale_pval=not emperical)
    
    adata.uns.update({"ligands_zscore": ligand_score.to_dict()})
    adata.uns.update({"receptors_zscore": receptor_score.to_dict()})
    adata.uns.update({"ligands_pval": ligand_pval.to_dict()})
    adata.uns.update({"receptors_pval": receptor_pval.to_dict()})
    adata.uns.update({"ligands_corr_pval": ligand_corr_pval.to_dict()})
    adata.uns.update({"receptors_corr_pval": receptor_corr_pval.to_dict()})
    adata.uns.update({"specificity_setting": settings})

    if return_values:
        return adata, ligand_values, receptor_values

    return adata

# Save and load specificity calculations (time consuming)

def save_specificity(adata, filename):
    """Saves data calculated by cn.connect.specificity to an excel file. 
    This file can later be loaded using cn.connect.load_specificity"""
    import pandas as pd
    
    keys = [
        'ligands',
        'receptors',
        'ligands_zscore',
        'receptors_zscore',
        'ligands_pval',
        'receptors_pval',
        'ligands_pval',
        'receptors_pval',
        'ligands_corr_pval',
        'receptors_corr_pval']
    
    xls = pd.ExcelWriter(filename)
    for key in keys:
        table = pd.DataFrame(adata.uns[key])
        table.to_excel(xls, sheet_name=key)
    s = pd.Series(adata.uns["specificity_setting"])
    s.to_excel(xls, sheet_name="specificity_setting")
    xls.close()

def load_specificity(adata, filename):
    """Loads previously calculated specificity to an andata object"""

    import pandas as pd
    keys = [
        'ligands',
        'receptors',
        'ligands_zscore',
        'receptors_zscore',
        'ligands_pval',
        'receptors_pval',
        'ligands_pval',
        'receptors_pval',
        'ligands_corr_pval',
        'receptors_corr_pval']

    for key in keys:
        data = pd.read_excel(filename, sheet_name=key, index_col=0)
        adata.uns[key] = {k:value.to_dict() for k, value in data.iteritems()}
    adata.uns["specificity_setting"] = pd.read_excel(filename, sheet_name="specificity_setting", index_col=0).to_dict()[0]
    