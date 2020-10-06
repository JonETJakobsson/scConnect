# Module deals with creation of ligand and receptor scores, and creation of scConnect tables etc.
import scConnect as cn
import scanpy as sc

version = cn.database.version
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
        utility = gmean(([x for x in [synthesis, transport, reuptake] if str(x) != "nan"])) # genes driving ligand production, remove nan values
        
        if str(utility) == "nan": # capture cases where no utility genes were present
            print(f"no utility genes detected for {ligand.ligand}")
            return 0.0  # exit before running exclusion calculation (division with 0)
            
        if utility == 0:  # if no utility genes are expressed, ligand score is 0
            return 0.0 # exit before running exclusion calculation (division with 0)

        exclusion_factor = ((utility - excluded) / utility) # correction factor for exclusion gene
        if exclusion_factor < 0: # exclusion should be 0 or positive
            exclusion_factor = 0.0
       
        ligand_score = utility * exclusion_factor # correct ligand expression based on the exclusion factor
        return ligand_score

    # If genes are missing from ligand gene list (TODO:this should not be the case!)
    else:
        print("Big error! ligand type is not defined!")
        return 0.0


def ligands(adata, organism="mmusculus", select_ligands=None):
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
    
    for cluster_data in adata.uns["gene_call"].iteritems():
        cluster = cluster_data[0]
        genes = cluster_data[1]
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


def receptors(adata, organism="mmusculus"):
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

    for cluster_data in adata.uns["gene_call"].iteritems():
        cluster = cluster_data[0]
        genes = cluster_data[1]
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

def interactions(emitter, target, self_reference=True, organism="mmusculus"):
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

    # Calculate total number of cluster combinations for the progress bar
    if self_reference is True:
        total_comb = len(list(product(emitter_clusters, target_clusters)))
    else:
        total_comb = len([(e, t) for (e, t) in product(
            emitter_clusters, target_clusters) if e != t])

    ligands = pd.DataFrame(emitter.uns["ligands"])
    receptors = pd.DataFrame(target.uns["receptors"])

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
                    emitter_ligands, target_receptors, interactions, interaction_set, receptor_info, ligand_info)
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

def get_connections(ligands, receptors, interactions, interaction_set, receptor_info, ligand_info):
    """finds connections between ligands and receptors 
    and return a score for each interaction"""

    from scipy.stats.mstats import gmean
    import numpy as np
    # shorten the list of interactions to only contain relevant ligands.
    # This should speed up the algorithm
    ligand_filter = [True if ligand in ligands.keys() else False for ligand in interactions.index.get_level_values(0)]
    interactions = interactions.loc[ligand_filter]

    connections = list()
    for ligand, l_score in ligands.iteritems():
        for receptor, r_score in receptors.iteritems():
            if (ligand, receptor) in interaction_set:
                interaction = interactions.loc[ligand, receptor]
                score = float(gmean((l_score, r_score)))
                connections.append((ligands.name, receptors.name, {
                    "score": score,
                    "log_score": np.log10(score + 1), # From here on, all values are +1ed and logaritmized with base of 10. # From here on, all values are +1ed and logaritmized with base of 10.
                    "ligand": ligand,
                    "receptor": receptor,
                    "interaction": f"{ligand} --> {receptor}",
                    "endogenous": f"{list(interaction.endogenous)}",
                    "action": f"{list(interaction.action)}",
                    "ligandspecies": f"{list(interaction.ligand_species)}",
                    "receptorspecies": f"{list(interaction.target_species)}",
                    "pubmedid": f"{list(interaction.pubmed_id)}",
                    "receptorfamily": receptor_info.loc[receptor]["family"],
                    "receptorgene": receptor_info.loc[receptor]["gene"],
                    "ligandtype": ligand_info.loc[ligand]["ligand_type"],
                    "ligandcomment": ligand_info.loc[ligand]["comment"]}))
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
        ligands = adata.uns["ligands"]
        receptors = adata.uns["receptors"]
        genes = adata.uns["gene_call"]
        clusters = ligands.keys()

        for cluster in clusters:
            print(f"processing cluster {cluster}")
            cluster_ligands = {k: v for k,
                               v in ligands[cluster].items() if v > 0}
            cluster_receptors = {k: v for k,
                                 v in receptors[cluster].items() if v > 0}
            node = (cluster, {"ligands": cluster_ligands,
                              "receptors": cluster_receptors,
                              "genes": genes[cluster].to_dict()},)
            nodes.append(node)

    return nodes
