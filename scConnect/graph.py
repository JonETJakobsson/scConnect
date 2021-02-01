import scConnect as cn
version = cn.database.version
organism = cn.database.organism

def build_graph(edge_list, node_list=False):
    """Construct a multi directional graph from an edgelist and nodelist and save as graph.
    Also return the graph object as a nx.MultiDiGraph.

    If no nodelist is passed, nodes will be implied from edgelist. 
    Create a nodelist using :py:func:`scConnect.connect.nodes`
    """
    import networkx as nx
    G = nx.MultiDiGraph()
    G.add_edges_from([(u, v, d) for u, v, d in edge_list])

    if isinstance(node_list, list):
        try:
            G.add_nodes_from(node_list)
        except:
            pass
    print(
        f"Graph has {G.number_of_edges()} interactions between {G.number_of_nodes()} clusters")

    nodes = list(G.nodes())
    for node in nodes:
        print(f"{node} has {G.degree(node)} interactions")

    # Calculate all edge and node metrics
    G = score(G)

    return G


# Simple function to abstact away some networkX API
def save(G, file_name):
    """Save the graph to .yaml format

    .. note::
        Yaml was the only format that supported utf-8 attributes and dictionary attributes
    """
    import networkx as nx
    nx.write_yaml(G, file_name, encoding="utf-8")


def load(filename):
    """load the graph from yaml format"""
    import networkx as nx
    G = nx.read_yaml(filename)
    return G
# Here follows functions to calculate edge atributes to help fins specifi interactions


def promiscuity(G):
    """
    Calculates edge attribute measure based on 
    the number of clusters a ligand is used to connect with.

    This measure is ligand specific.
    """
    ligand_count = dict()
    for n, nbrs in G.adj.items():
        cluster_list = dict()
        for nbr, edict in nbrs.items():
            # We do only want to check if a ligand is connected to the neighbor once.
            ligand_check = dict()
            for e, d in edict.items():
                if d["ligand"] not in ligand_check.keys():
                    if d["ligand"] in cluster_list.keys():
                        cluster_list[d["ligand"]] += 1
                    else:
                        cluster_list[d["ligand"]] = 1
                    # Lets stop counting this ligand
                    ligand_check[d["ligand"]] = "counted for this neighbor"
        ligand_count[n] = cluster_list
    # Adds normalized promiscuity value for each edge depending on ligand used.
    tot_n = len(G.nodes())
    for n, nbrs in G.adj.items():
        for nbr, edict in nbrs.items():
            for e, d in edict.items():
                G.edges[n, nbr, e]["promiscuity"] = ligand_count[n][d["ligand"]]/tot_n

    return G


def loyalty(G):
    """
    calculates edge attribute measure based on ratio of interaction
    to the neighbor over total interactions to all neighbors

    this measure is neighbour specific.
    """
    score_dict = dict()
    for n, nbrs in G.adj.items():
        scores = dict()
        for nbr, edict in nbrs.items():
            scores[nbr] = len(edict)
        # normalize on total interactions
        total = sum(scores.values())
        for nbr in scores.keys():
            scores[nbr] /= total
        score_dict[n] = scores

    # Adds the loyalty score to each edge depending on neigbour
    for n, nbrs in G.adj.items():
        for nbr, edict in nbrs.items():
            for e in edict.keys():
                G.edges[n, nbr, e]["loyalty"] = score_dict[n][nbr]
    return G


def weighted_score(G):
    """
    Calculate a weighted score based on original score,
    loyalty and promiscuity for each interaction.

    High loyalty promotes interactions between a set of neighbours
    High promiscuity punish interactions of specific ligands.
    """
    from math import log1p
    from numpy import log10
    for n, nbrs in G.adj.items():
        for nbr, edict in nbrs.items():
            for e, d in edict.items():
                # punnish interaction score low loyalty and for high promiscuity
                weighted_score = d["score"] * d["loyalty"] / d["promiscuity"]
                # Add w_score to the edge attribute
                G.edges[n, nbr, e]["weighted_score"] = weighted_score
                G.edges[n, nbr, e]["log_weighted_score"] = log10(weighted_score)
    return G


def centrality(G):
    """
    Calculate information centrality for all nodes in the graph

    returns a graph with added node attributes"""
    import networkx as nx

    cent = nx.degree_centrality(G)
    for node in G.nodes():
        G.nodes[node]["centrality"] = cent[node]

    return G


# Summary function for all graph measures (used to quickly run all)


def score(G):
    """
    Calculate all edge and node attribute measures for interactions in the graph.

    This needs to be run before flatten graph and is automatically
    run when constructing a graph with :py:func:`scConnect.graph.build_graph`
    """
    # Edge attributes
    G = promiscuity(G)
    G = loyalty(G)
    G = weighted_score(G)

    # Node attributes
    G = centrality(G)
    # To be added if needed
    return G

# The flatten graph aproach has been depreceated.
# Use conversion to adjacency matrix on filtered graphs insted

# filter the graph based on interactions containing ligands and receptors under a threshold p-value
def significant_interactions(G, a=0.05):
    """Return a filtered graph containing only sigificanly differentiated ligands and receptor interaction.
    
    a: alpha, selects scores with a pvalue under alpha.
    Return a MultiDiGraph
    """

    import networkx as nx
    G_filtered = nx.MultiDiGraph()
    for u, v, n, d in G.edges(data=True, keys=True):
        # Find interactions with significant p-values (pval < a)
        # Must be significantly upregulated (zscore > 0)
        if d["ligand_pval"] < a and d["receptor_pval"] < a and d["ligand_zscore"] > 0 and d["receptor_zscore"] > 0:
            G_filtered.add_edge(u, v , key=n ,**d)
            
    G_filtered.add_nodes_from(G.nodes(data=True))
    
    return G_filtered
    


def flatten_graph(G, weight="score", log=True):
    """
    Collapse graph and calculate edge weight between clusters as:
    weight = sum(weighted_score)

    log: log the sum of scores

    Returns: a nx.DiGraph

    .. note::
        Use flatten graph to export to a simpler format, for instance
        if working with gephi or other softwares.
        Note that edge atributes for specific interactions are
        stored as a list of dictionaries, but might be 
        unaccessible for most programs.

    """
    import networkx as nx
    from numpy import log1p, log10
    G_flat = nx.DiGraph()
    for n, nbrs in G.adj.items():
        for nbr, edict in nbrs.items():
            interaction_scores = list()  # collect all weighted interaction scores in this list
            for e, d in edict.items():
                interaction_scores.append(edict[e][weight])

            weight_sum = sum(interaction_scores)
            if log:
                weight_sum = log10(weight_sum +1)                

            # Add edge and pass in edict for possibility to access later
            G_flat.add_edge(n, nbr, weight=weight_sum,
                            interactions=list(edict.values()))

    # Forward node attributes
    G_flat.add_nodes_from(G.nodes(data=True))
    return G_flat


# Functions to access underlying interactions in a flattened graph
def get_interactions(G, n, nbr):
    """
    Find interactions between nodes in a flattened graph
    and returns a list sorted on score
    """
    edge = G[n][nbr]
    import pandas as pd
    interactions = pd.DataFrame(edge["interactions"])
    interactions.sort_values("score", ascending=False, inplace=True)

    return interactions


def get_emmitors(G, n):
    """
    Find all clusters emitting to selected cluster and
    returns a dictionary of emittor weight pairs
    """
    import networkx as nx
    emittors = {emittor: data["weight"] for (
        emittor, target, data) in G.edges(data=True) if target == n}
    return emittors


def get_targets(G, n):
    """
    Find all clusters targeted by selected cluster and
    returns a dictionary of target weight pairs
    """
    import networkx as nx
    targets = {target: data["weight"] for (
        emittor, target, data) in G.edges(data=True) if emittor == n}
    return targets


def flatten_graph_to_adjacency_matrix(G, weight="Weighted_score", log=False):
    """flattens graph over provided weight and returns a sorted adjacency matrix
    if log = True, performs a log10 transformation of the results
    """
    import networkx as nx
    import numpy as np
    import pandas as pd

    if len(G.edges()) == 0:
        matrix = np.zeros(shape=(len(G.nodes()), len(G.nodes())))
        matrix = pd.DataFrame(matrix, columns=G.nodes(), index=G.nodes())
    else:
        matrix = nx.to_pandas_adjacency(G, weight=weight)
        ind = matrix.index.values
        matrix = matrix.sort_index(ascending=True)
        matrix = matrix[matrix.index.values]
        matrix = pd.DataFrame(matrix, columns=ind, index=ind)

    if log:
        matrix = np.log10(matrix)
    return matrix


def split_graph(G):
    """splits graph(s) on interactions and return a dictionary of graphs with interaction as keys."""
    import networkx as nx

    # Find all interactions to split the graph on
    split_by = "interaction"
    split_list = list()
    for u, v, d in G.edges(data=True):
        split_list.append(d[split_by])
    split_set = set(split_list)
    G_splits = dict()
    for split in split_set:
        G_split = nx.from_edgelist(
            [(u, v, d) for u, v, d in G.edges(data=True) if d[split_by] == split],
            create_using=nx.MultiDiGraph())

        G_split.add_nodes_from(G.nodes(data=True))
        G_splits[split] = G_split

    return G_splits


def plot_adjacency_matrix(G, save=False, filename=None, weight="weighted_score", log=False, width=200, height=200, fontsize=10, fontscale=1.5):
    """Flattens and plots an adjecency matric for a given graph.
    If a split graph is provided, plots adjacency matrix for each interaction.

    Weight can be "weighted_score" or "score".
    Returns: Holoviews heatmap
    """
    import numpy as np
    import holoviews as hv
    from holoviews import opts
    import networkx as nx
    hv.extension('bokeh')

    # Create adjecency matrix and sort index and columns by name
    matrix = nx.to_pandas_adjacency(G, weight=weight)

    if log:
        matrix = np.log10(matrix)

    matrix = matrix.sort_index(ascending=True)
    matrix = matrix[matrix.index.values]
    matrix = matrix.stack()
    matrix = matrix.reset_index(name="values")
    matrix.columns = ["source", "target", "values"]

    # Plot the matrix as a heatmap
    heatmap = hv.HeatMap(matrix, kdims=["target", "source"])

    # Rotate X lables and show the color bar
    heatmap.opts(
        opts.HeatMap(tools=["hover"], xrotation=90, colorbar=True, width=width, height=height, fontsize=fontsize, fontscale=fontscale)
    )

    if save is True:
        hv.save(heatmap, filename)
    return heatmap


def corr_adjacency_matrix(G1, G2, weight="score"):
    """finds the pearson correlation between two adjecency matrix.

    Returns: Pearson correlation coeficient"""
    import networkx as nx
    import numpy as np

    # Create adjacency matric from Graph 1 and sort index and columns
    # by name and flatten the matrix to an array
    m1 = nx.to_pandas_adjacency(G1, weight=weight)
    m1 = m1.sort_index(ascending=True)
    m1 = m1[m1.index.values]
    m1 = m1.values.flatten()

    # Create adjacency matric from Graph 2 and sort index and columns
    # by name and flatten the matrix to an array
    m2 = nx.to_pandas_adjacency(G2, weight=weight)
    m2 = m2.sort_index(ascending=True)
    m2 = m2[m2.index.values]
    m2 = m2.values.flatten()

    # Find the pearson correlation coefficiant between the two arrays
    corr = np.corrcoef(m1, m2)[0][1]
    return corr


def plot_chord(G, notebook=False):
    """
    render an interactive chord plot of the graph.

    If notebook is False, starts an bokeh app in a browser window, 
    if True, renders the plot directly in the cell.
    """

    import scConnect as cn
    import holoviews as hv
    from holoviews import opts, dim
    import networkx as nx
    import pandas as pd
    import numpy as np

    # instantiate the bokeh renderer
    renderer = hv.renderer('bokeh')
    hv.extension("bokeh")
    hv.output(size=250)

    # set visuals
    opts.defaults(
        opts.Chord(
            node_cmap='Category20',
            edge_cmap='Category20',
            edge_color=dim("source"),
            labels='cluster',
            node_color=dim('cluster'),
            inspection_policy="edges",
            toolbar="above",
        )
    )

    # Creat the Dataset object to be passed to the Chord object (NOTE: aggregates the data
    # leaving only one edge per cluster pair)
    edges = nx.to_pandas_edgelist(G)
    links = hv.Dataset(edges, ["source", "target"], ["weighted_score", "loyalty"]).sort(
        by="source").aggregate(function=np.sum)
    nodes = hv.Dataset(list(G.nodes), 'cluster').sort(by="cluster")

    # Calculate values for threshold representing given percentiles
    percentiles = [0, 20, 40, 60, 80, 90, 95, 99]
    th_values = np.percentile(links.data["weighted_score"], percentiles)

    th = hv.Dimension(("th", "weighted scores threshold"),
                      default=th_values[0])

    # Filter data on threshold, and return a chord element
    def chord_graph(th):
        links_ = links.select(weighted_score=(th, None))
        chord = hv.Chord((links_, nodes))

        return chord

    # instantiate the dynamic map
    chord = hv.DynamicMap(chord_graph, kdims=[th]).redim.values(th=th_values)

    # Run the server if not in notebook
    if notebook == False:
        server = renderer.app(chord, show=True, new_window=True)

    if notebook == True:
        return chord


def plot_sankey(G, notebook=False):
    """
    render an interactive sankey plot of the graph.

    If notebook is False, starts an bokeh app in a browser window, 
    if True, renders the plot directly in the cell.
    """

    import scConnect as cn
    import holoviews as hv
    from holoviews import opts, dim
    import networkx as nx
    import pandas as pd
    import numpy as np

    # instantiate the bokeh renderer
    renderer = hv.renderer('bokeh')
    hv.extension("bokeh")
    hv.output(size=100)

    # set visuals
    opts.defaults(
        opts.Sankey(
            cmap='Category20',
            edge_cmap='Category20',
            edge_color=dim("receptorfamily"),
            labels='cluster',
            node_color=dim('cluster'),
            inspection_policy="edges",
            selection_policy="edges",
            colorbar=True,
            toolbar="above"
        ),
        opts.Bars(
            invert_axes=False,
            xrotation=70,
            toolbar="above"
        )
    )

    # Create values to be used for filtering of the graph
    edges = nx.to_pandas_edgelist(G)
    percentiles = [0, 20, 40, 60, 80, 90, 95, 99]
    th_values = np.percentile(edges["weighted_score"], percentiles)
    nodes_list = list(G.nodes())

    node = hv.Dimension(("node", "Focus cluster"), default=nodes_list[0])
    th = hv.Dimension(('th', 'Weighted score threshold'), default=th_values[0])

    # Filter data on node and threshold, and return a sankey element
    def sankey_graph(node, th):
        # Find all interactions where node is target or source node
        G_s = nx.MultiDiGraph()
        for n, nbrs in G.adj.items():
            for nbr, edict in nbrs.items():
                if n == node:
                    for e, d in edict.items():
                        # append dash after the target node
                        G_s.add_edge(n, nbr + "_", **d)
                if nbr == node:
                    for e, d in edict.items():
                        # append dash before the source node
                        G_s.add_edge("_" + n, nbr, **d)
        # create the dataset used to build the sankey graph.
        # Sort values on weight to get ordered representation on plot.
        edges = nx.to_pandas_edgelist(G_s)
        links = hv.Dataset(edges, ["source", "target"],
                           ["weighted_score", "interaction", "receptorfamily"]).sort("weighted_score")
        nodes = hv.Dataset(list(G_s.nodes), 'cluster')
        sankey = hv.Sankey((links, nodes)).select(weighted_score=(th, None))

        # calculate bars
        ligands = hv.Dataset(edges, ["ligand", "source"], ["score"]).select(
            source=node).aggregate(function=np.mean).sort("score", reverse=True)
        receptors = hv.Dataset(edges, ["receptor", "target"], ["score"]).select(
            target=node).aggregate(function=np.mean).sort("score", reverse=True)
        bars = (hv.Bars(ligands, "ligand") +
                hv.Bars(receptors, "receptor")).cols(2)

        # calculate table
        ligands = hv.Dataset((G.node[node]["ligands"]), "ligand", "score").sort(
            "score", reverse=True)
        receptors = hv.Dataset((G.node[node]["receptors"]), "receptor", "score").sort(
            "score", reverse=True)
        table = hv.Layout(hv.Table(ligands) + hv.Table(receptors)).cols(2)

        return (bars + table + sankey).cols(2)

    sankey = hv.DynamicMap(sankey_graph, kdims=[node, th]).redim.values(
        node=nodes_list, th=th_values)

    layout = sankey
    # Run the server if not in notebook
    if notebook == False:
        server = renderer.app(layout, show=True, new_window=True)

    if notebook == True:
        return layout




def compare_interactions_df(G, node_a=str, node_b=str, method="ratio"):
    
    import networkx as nx
    import pandas as pd
    import numpy as np
    
    Gs = split_graph(G)
    interaction_dict_i = dict()
    for node in G.nodes():
        ratio_dict = dict()
        for interaction in Gs:
            
            if method == "ratio":
                df = nx.to_pandas_adjacency(Gs[interaction], weight="log_score")
                a = df[node_a][node]
                b = df[node_b][node]
                ratio = a - b   # calculate the ratio of interaction scores (using log values log(a)-log(b) = log(a/b))
                                # This bypass divition with 0 issues as in a/b, b might be 0

            if method == "difference":
                df = nx.to_pandas_adjacency(Gs[interaction], weight="score")
                a = df[node_a][node]
                b = df[node_b][node]
                ratio = a - b  # Calulate the difference in interaction score

            if method == "specificity":
                df = nx.to_pandas_adjacency(Gs[interaction], weight="specificity")
                a = df[node_a][node]
                b = df[node_b][node]
                ratio = a - b  # Calulate the difference in specificity

            ratio_dict[interaction] = ratio
        interaction_dict_i[node] = ratio_dict
    df_i = pd.DataFrame(interaction_dict_i)


    interaction_dict_o = dict()
    for node in G.nodes():
        ratio_dict = dict()
        for interaction in Gs:
            if method == "ratio":
                df = nx.to_pandas_adjacency(Gs[interaction], weight="log_score")
                a = df[node][node_a]
                b = df[node][node_b]
                ratio = a - b   # calculate the ratio of interaction scores (using log values log(a)-log(b) = log(a/b))
                                # This bypass divition with 0 issues as in a/b, b might be 0
                
            if method == "difference":
                df = nx.to_pandas_adjacency(Gs[interaction], weight="score")
                a = df[node][node_a]
                b = df[node][node_b]
                ratio = a - b  # Calulate the difference in interaction score
                
            ratio_dict[interaction] = ratio
        interaction_dict_o[node] = ratio_dict
    df_o = pd.DataFrame(interaction_dict_o)
    return df_i, df_o




def compare_interactions_plot(G, node_a=str, node_b=str, th=float, figsize_i=(5, 15), figsize_o=(5,15), colormap="nipy_spectral", method="ratio", save=None):
    """Compares interactions scores between node a <--> x and node b <--> x. returns the logaritmized ratio.
    
    if method = "ratio": log(a/b) = log(a)-log(b), "th": times difference ie. 4x
                "difference": values = a-b, "th": difference ie. 200
     

    returns two plots, first incomming interactions to a an b, then outgoing interaction from a and b. 
    also a lut and color pallet for the receptor families
    """
    import seaborn as sns
    import numpy as np
    import scConnect as cn
    import matplotlib.pyplot as plt
    
    df_i, df_o = compare_interactions_df(G, node_a, node_b, method)

    if method == "ratio": # convert a ratio th to the corresponding log values
        th = np.log10(th)


    # Filter interactions
    interaction_filter_i = [True if max(abs(values)) > th else False for i, values in df_i.iterrows()]
    df_i = df_i[interaction_filter_i]
    interaction_filter_o = [True if max(abs(values)) > th else False for i, values in df_o.iterrows()]
    df_o = df_o[interaction_filter_o]

    # Find receptors used by both incomming and outgoing interactions
    receptors_i = [interaction.split("--> ")[1] for interaction in df_i.index]
    receptors_o = [interaction.split("--> ")[1] for interaction in df_o.index]
    receptors = set(receptors_i + receptors_o)
    
    # find receptor families for all selected interactions to color rows
    targets = cn.database.get_data("targets")
    targets = targets.set_index("Target name")

    # create a dictionary of receptors and families and a correspinding lut
    families_dict = targets.loc[receptors]["Family name"].to_dict() # selects only families for which we have receptors in the dataset
    unique_families = set(families_dict.values())
    lut = dict(zip(unique_families, sns.color_palette(colormap, len(unique_families))))

    if df_i.shape[0] > 1:
        # Make color mapping for incomming interactions
        row_colors_i = list()
        for receptor in receptors_i:
            family = families_dict[receptor]
            color = lut[family]
            row_colors_i.append(color)

        # plot a clustered heatmap with receptor families as cluster row color for incomming interactions
        print("incomming interactions")

    
    
        sns.clustermap(df_i, center=0, cmap="bwr", figsize=figsize_i, row_colors=row_colors_i)
        if save != None:
            path_i = f"{save}{node_a}_to_{node_b}_incomming.pdf"
            plt.savefig(path_i)

    if df_o.shape[0] > 1:
        # Make color mapping for outgoing interactions
        row_colors_o = list()
        for receptor in receptors_o:
            family = families_dict[receptor]
            color = lut[family]
            row_colors_o.append(color)

        # plot a clustered heatmap with receptor families as cluster row color
        print("Outgoing interactions")
        sns.clustermap(df_o, center=0, cmap="bwr", figsize=figsize_o, row_colors=row_colors_o)
        if save != None:
            path_o = f"{save}{node_a}_to_{node_b}_outgoing.pdf"
            plt.savefig(path_o)

    print(f"receptor families in order: {list(lut.keys())}")
    sns.palplot(lut.values())

    return None


def plot_interaction_heatmap(G, th=0.5, method='ward', metric='euclidean', z_score=None, standard_scale=None, figsize=None, save=None):
    """plots a heatmap of logaritmized interaction scores for each Source-->Target combination. 
    
    
    G: multi-directional graph
    th: filter out all interactions where the maximum interaction score is lower than the threshold (log values)
    method: passed to seaborn.clustermap
    metric: passed to seaborn.clustermap
    z_score: passed to seaborn.clustermap. calculates z-scores for each value x in the matrix as z=(x-mean)/std.
    standard_scale: passed to seaborn.clustermap
    figsize: passed to seaborn.clustermap. this is automatically set based on size of dataset, but can be explicity set using a tuple.
    save: filename with file extension

    return: seaborn.matrix.ClusterGrid object
    """

    import networkx as nx
    import pandas as pd
    import seaborn as sns
    import itertools

    # split graph on interactions
    print("split graph...")
    Gs = split_graph(G)
    
    ams = dict()
    print("flatten interactions...")
    for interaction in Gs.keys():
        g = Gs[interaction]
        am = nx.to_pandas_adjacency(g, weight="log_score") #make adjacency matrix with scores
        am = pd.Series(am.values.flatten()) # flatten the array
        ams[interaction] = am # add flatten adjacency to dict
    
    print("recovering index...")
    #create index
    am = nx.to_pandas_adjacency(Gs[list(Gs.keys())[0]])
    ix = list(itertools.product(am.index.values, repeat=2))
    ix = [f"{pops[0]}-->{pops[1]}" for pops in ix]

    # Create a df with the flatten arrays
    print("build dataframe...")
    df = pd.DataFrame(ams)
    df.index = ix
    filter_interactions = df.aggregate(max, axis=0) > th
    df = df.loc[:,filter_interactions.values]

    print("plotting...")
    if figsize == None:
        figsize = (len(df.index)/3, len(df.columns)/3)
    f = sns.clustermap(df.T, figsize=figsize, method=method, metric=metric, z_score=z_score, standard_scale=standard_scale)
    f.savefig(save)
    
    return f


# Fetch edge list dataframe for plotting functions
def edge_list(
    G, 
    filter_by="specificity", 
    th=1, 
    cellphonedb_interactions=False, 
    interaction_list=None,
    version = version,
    organism = organism,
    all_connections=True):
    """retrieves all interactions from the graphs and returns a pandas dataframe suitable for plotting"""

    import networkx as nx
    import pandas as pd
    import numpy as np
    import plotly.express as px
    import pkg_resources
    
    
    df = nx.to_pandas_edgelist(G)
    
    # change name of scConnects interaction to better compare to cellPhoneDB
    receptor_info = pd.read_csv(pkg_resources.resource_filename("scConnect", (f"data/Gene_annotation/{version}/{organism}/receptors.csv")), index_col=1)
    ligand_info = pd.read_csv(pkg_resources.resource_filename("scConnect", (f"data/Gene_annotation/{version}/{organism}/ligands.csv")), index_col=1)

    df["ligand_gene"] = [eval(ligand_info.loc[ligand]["preprogene"])[0] if isinstance(ligand_info.loc[ligand]["preprogene"], str) else ligand for ligand in df["ligand"]]
    df["receptor_gene"] = [eval(receptor_info.loc[receptor]["gene"])[0] if isinstance(receptor_info.loc[receptor]["gene"], str) else receptor for receptor in df["receptor"]]

    df["cellphonedb_interactions"] = [f"{l}_{r}" for l, r in list(zip(df["ligand_gene"], df["receptor_gene"]))]

    #Create a new connection annotation merging source and target info
    df["connection"] = [str(connection[0]+"|"+connection[1]) for connection in zip(df["source"], df["target"])]

    #Collect all connection for specific interactions
    if all_connections:
        # select all significant interaction (significance > 1)
        interactions = df[df[filter_by] >= th]["interaction"].unique()
        # keep all interactions that were significant between any connection
        df = df[[True if interaction in interactions else False for interaction in df["interaction"]]]

        # only sort if all_connections == True
         # sort based on connection and interactoin names to get order in the plot
        if cellphonedb_interactions:
            df.sort_values(by=[ "cellphonedb_interactions", "connection"], key=np.vectorize(str.lower), ascending=False, inplace=True)
        else:
            df.sort_values(by=["connection", "interaction"], key=np.vectorize(str.lower), ascending=False, inplace=True)
    
    else:
        # if only selected interactions, sort by filter_by key
        df= df[df[filter_by] >= th]
        df.sort_values(by=[filter_by], ascending=False, inplace=True)
   
    return df

# plotting function using the edge list data
def dotplot(
    G=None, 
    df=None, 
    filter_by="specificity", 
    th=1,
    organism=organism,
    version=version,
    all_connections=True,
    cellphonedb_interactions=False, 
    height_scale=40, 
    width_scale=40, 
    cmap=None):

    import plotly.express as px
    import pandas as pd
    import numpy as np

    if cmap == None:
        cmap = px.colors.sequential.Viridis_r
    
    if isinstance(df, pd.DataFrame):
        df = df
        # Filter by feature and threshold
        interactions = df[df[filter_by] >= th]["interaction"].unique()
        # keep all interactions that were significant between any connection
        df = df[[True if interaction in interactions else False for interaction in df["interaction"]]]
    else:
        df = edge_list(G, filter_by, th, organism=organism, version=version, all_connections=all_connections)
        
    # set width and hight of figure
    width = len(df["connection"].unique())*width_scale
    height = len(df["interaction"].unique())*height_scale
    if cellphonedb_interactions:
        y = "cellphonedb_interactions"
    else:
        y = "interaction"
    
    # log (specificity +1)for plotting sized
    df["log_specificity"] = np.log10(df["specificity"]+1)
    
    # find and sort all connections
    connections = list(set(df["connection"]))
    connections.sort()

    
    
    plot = px.scatter(df, x="connection", y=y, 
        size="log_specificity", 
        color="log_score",
        opacity=0.8,
        height=height, width=width,
        size_max=10,
        template="plotly_white",
        render_mode="svg",
        title=f"Interactions with {filter_by} higher than {th}",
        color_continuous_scale=cmap,
        hover_data=["specificity"],
        category_orders = {"connection": connections})
    
    plot.update_layout(coloraxis_colorbar=dict(
        title="Log score",
        thicknessmode="pixels", thickness = 20,
        lenmode="pixels", len=200,
        yanchor="top", y=1
        )
    )
    return plot
    

# functions to retrieve dataframes of ligand and receptor
# information for further analysis
def get_ligand_df(G, color_map=False, corr_pval=True):
    """fetch all ligand information in a graph and return a DataFrame in long format"""
    import pandas as pd
    import numpy as np

    nodes = G.nodes(data=True)
    df = pd.DataFrame(columns=["ligand", "score", "node"])
    for node, data in G.nodes.items():
        temp_df = pd.DataFrame.from_dict(data["ligands_score"], orient="index", columns=["score"])
        temp_df["z_score"] = pd.DataFrame.from_dict(data["ligands_zscore"], orient="index", columns=["z-score"])
        if corr_pval:
            temp_df["pvalue"] = pd.DataFrame.from_dict(data["ligands_corr_pval"], orient="index", columns=["pvalue"])
        else:
             temp_df["pvalue"] = pd.DataFrame.from_dict(data["ligands_pval"], orient="index", columns=["pvalue"])
        temp_df["specificity"] = -np.log10(temp_df["pvalue"])
        temp_df["log_score"] = np.log10(temp_df["score"]+1)
        temp_df.rename_axis(index="ligand", inplace=True)
        temp_df["node"] = node
        temp_df.reset_index(inplace=True)
        df = df.append(temp_df, ignore_index=True,)
        
    color = dict(G.nodes(data="color"))
    for k, v in color.items():
        color[k] = "rgb" + str(tuple(v))
        
    if color_map == True:
        return df, color
    
    return df

def get_receptor_df(G, color_map=False, corr_pval=True):
    """fetch all ligand information in a graph and return a DataFrame in long format"""
    import pandas as pd
    import numpy as np

    nodes = G.nodes(data=True)
    df = pd.DataFrame(columns=["receptor", "score", "node"])
    for node, data in G.nodes.items():
        temp_df = pd.DataFrame.from_dict(data["receptors_score"], orient="index", columns=["score"])
        temp_df["z_score"] = pd.DataFrame.from_dict(data["receptors_zscore"], orient="index", columns=["z-score"])
        if corr_pval:
            temp_df["pvalue"] = pd.DataFrame.from_dict(data["receptors_corr_pval"], orient="index", columns=["pvalue"])
        else:
            temp_df["pvalue"] = pd.DataFrame.from_dict(data["receptors_pval"], orient="index", columns=["pvalue"])
        temp_df["specificity"] = -np.log10(temp_df["pvalue"])
        temp_df["log_score"] = np.log10(temp_df["score"]+1)
        temp_df.rename_axis(index="receptor", inplace=True)
        temp_df["node"] = node
        temp_df.reset_index(inplace=True)
        df = df.append(temp_df, ignore_index=True,)
        
    color = dict(G.nodes(data="color"))
    for k, v in color.items():
        color[k] = "rgb" + str(tuple(v))
        
    if color_map == True:
        return df, color
    
    return df