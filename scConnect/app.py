
# Dash app to investigate graphs


def graph(G, mode=None):
    """
    G: a multidirectional graph
    mode: inline to show inside the jupyter nodebook, default is None
    """

    import dash
    from jupyter_dash import JupyterDash
    import dash_cytoscape as cyto
    from dash.dependencies import Output, Input
    import dash_html_components as html
    import dash_core_components as dcc
    import dash_table
    import networkx as nx
    import scConnect as cn
    import plotly.graph_objs as go
    import pandas as pd
    import numpy as np
    import json
    import matplotlib
    import matplotlib.pyplot as plt

    cyto.load_extra_layouts()

    JupyterDash.infer_jupyter_proxy_config()

    app = JupyterDash(__name__)

    server = app.server
    # Add a modified index string to change the title to scConnect
    app.index_string = '''
        <!DOCTYPE html>
        <html>
            <head>
                {%metas%}
                <title>scConnect</title>
                {%favicon%}
                {%css%}
            </head>
            <body>
                {%app_entry%}
                <footer>
                    {%config%}
                    {%scripts%}
                    {%renderer%}
            </body>
        </html>
        '''
    # Add colors to each node
    nodes = pd.Categorical(G.nodes())
    # make a list of RGBA tuples, one for each node
    colors = plt.cm.tab20c(nodes.codes/len(nodes.codes), bytes=True)
    # zip node to color
    color_map_nodes = dict(zip(nodes, colors))

    # add these colors to original graph
    for node, color in color_map_nodes.items():
        G.nodes[node]["color"] = color[0:3]  # Save only RGB

    # Add colors to edges(source node color) for  G
    for u, v, k in G.edges(keys=True):
        G.edges[u,v,k]["color"] = color_map_nodes[u][0:3]

    #load graph into used formes
    G_flat = cn.graph.flatten_graph(G, weight="score", log=True)
    
    # Add colors to edges(source node color) for G_flat
    for u, v,  in G_flat.edges():
        G_flat.edges[u,v]["color"] = color_map_nodes[u][0:3]


    G_split = cn.graph.split_graph(G)
      
    # find and sort all found interactions
    interactions = list(G_split.keys())
    interactions.sort()

    
    G_cyto = nx.cytoscape_data(G_flat)

    # get min and max weight for all edges for flat and normal graph
    #weights = [d["weight"] for u, v, d in G_flat.edges(data=True)]
    scores = [d["score"] for u, v, d in G.edges(data=True)] + [d["weighted_score"] for u, v, d in G.edges(data=True)]
    cent = [d["centrality"] for n, d in G.nodes(data=True)]


    # prepare data for network graph
    nodes = G_cyto["elements"]["nodes"]
    elements = []

    # collect all available genes
    genes = list(nodes[0]["data"]["genes"].keys())


    # Styling parameters
    font_size = 20

    # Style for network graph
    default_stylesheet = [
        {
            'selector': 'node',
            'style': {
                'background-color': 'data(color)',
                'label': 'data(id)',
                'shape': 'ellipse',
                'opacity': 1,
                'font-size': f'{font_size}',
                'font-weight': 'bold',
                'text-wrap': 'wrap',
                'text-max-width': "100px",
                'text-opacity': 1,
                'text-outline-color': "white",
                'text-outline-opacity': 1,
                'text-outline-width': 2
            }
        },

        {
            'selector': 'node:selected',
            'style': {
                'background-color': 'data(color)',
                'label': 'data(id)',
                'shape': 'ellipse',
                'opacity': 1,
                'border-color': "black",
                'border-width': "5"
            }
        },
        
        {
            'selector': 'edge',
            'style': {
                'line-color': 'data(color)',
                "opacity": 0.7,
                "curve-style": "unbundled-bezier",
                "width": "data(weight)",
                "target-arrow-shape": "vee",
                "target-arrow-color": "black",
                'z-index': 1,
                'font-size': f'{font_size}'
            }
        },

        {
            'selector': 'edge:selected',
            'style': {
                'line-color': 'red',
                'line-style': "dashed",
                'opacity': 1,
                'z-index': 10,
            }
        }
    ]
    app.layout = html.Div(className="wrapper", children=[  # wrapper

        html.Div(className="header", children=[  # header
            html.Img(src="assets/logo.png", alt="scConnect logo"),
            
            html.Div(className="graph-info", id="graph-stat", children=[
                html.H3(f'Loaded graph with {len(G.nodes())} nodes and {len(G.edges())} edges')
            ])
        ]),

        html.Div(className="network-settings", children=[  # network settings
            html.H2("Network settings"),

            html.Label("Interactions"),
            dcc.Dropdown(
                id="network-interaction",
                options=[{'label': "all interactions", 'value': "all"}]+
                [{'label': interaction, 'value': interaction} for interaction in interactions],
                value="all"
            ),

            html.Label("Graph Layout"),
            dcc.Dropdown(
                id="network-layout",
                options=[
                    {'label': name.capitalize(), 'value': name}
                    for name in ['grid', 'random', 'circle', 'cose', 'concentric', 'breadthfirst', 'cose-bilkent', 'cola', 'euler', 'spread', 'dagre', 'klay']
                ],
                value="circle",
                clearable=False),

            html.Label("Weight Filter", style={"paddingBottom": 500, "paddingTop": 500}),
            dcc.Slider( # min, max and value are set dynamically via a callback
                id="network-filter",
                step=0.001,
                updatemode="drag",
                tooltip={
                    "always_visible": True,
                    "placement": "right"
                },
                ),

            html.Label("Node size"),
            dcc.RangeSlider(
                id="node-size", 
                value=[10,50], 
                min=0, 
                max=100,
                updatemode="drag"), 

            html.Label("Select gene"),
            dcc.Dropdown(
                id="gene_dropdown", 
                options=[{"label": gene, "value": gene} for gene in genes],
                clearable=True,
                placeholder="Color by gene expression",
                ),

            #Store node colors "hidden" for gene expresison
            html.Div(id="node-colors", style={"display": "none"}, children=[
                ""
            ]),

            html.Div(id="min-max", children=[

            ])

        ]),  # end network settings
        html.Div(id="network-graph", className="network-graph", children=[  # network graph
            cyto.Cytoscape(
                id="cyto-graph",
                style={
                    'width': '100%',
                    'height': '80vh'},
                stylesheet=default_stylesheet,
                elements=elements,
                autoRefreshLayout=True,
                zoomingEnabled=False)
        ]),  # end network graph

        html.Div(className="sankey-settings", children=[  # network settings
            html.H2("Connectivity Settings"),

            html.Label("Weight Filter"),
            dcc.Slider(
                id="sankey-filter", 
                min=min(scores),
                max=max(scores), 
                value=0.75, 
                step=0.001,
                updatemode="drag",
                tooltip={
                    "always_visible": True,
                    "placement": "right"
                }),
            
            html.Label("Toggle weighted"),
            dcc.RadioItems(id="sankey-toggle", options=[
                {"label": "Score", "value": "score"},
                {"label": "Weighted score", "value": "weighted_score"},
                {"label": "Log score", "value": "log_score"}
            ], value="score")
            
        ]),  # end network settings

        html.Div(className="sankey", id="sankey", children=[  # sankey graph
            dcc.Graph(id="sankey-graph")
        ]),  # end sankey graph

        html.Div(className="interaction-list", children=[  # interaction list

            html.Div(id="selection", children=[

                html.H3(id="edge-info"),

                dash_table.DataTable(id="edge-selection",
                                     style_table={
                                         "overflowX": "scroll",
                                         "overflowY": "scroll",
                                         "height": "50vh",
                                         "width": "95%"

                                     },
                                     style_cell={
                                         "minWidth": "0px",
                                         "overflow": "hidden"
                                     },
                                     sort_action="native",
                                     fixed_rows={'headers': True, 'data': 0}
                                     )
            ])
        ]),  # end interaction list

        html.Div(className="L-R-settings", children=[ # ligand and receptor settings (search)
            html.H2("Ligands and receptors"),
            html.Label("Search for ligands and receptors:"),
            dcc.Input(id="filter_l_r", type="search", value="", placeholder="Search"),
        ]),

        html.Div(className="L-R-scores", children=[  # ligand and receptor lists
            dcc.Tabs([
                dcc.Tab(label="Ligands", children=[
                    dcc.Graph(
                        id="ligand-graph",
                        config=dict(
                            autosizable=True,
                            responsive=True)
                        )
                ]),
                dcc.Tab(label="Receptors", children=[
                    dcc.Graph(
                        id="receptor-graph",
                        config=dict(
                            autosizable=True,
                            responsive=True)
                    )
                ])
            ])
        ]) # end ligand receptor list
    ])  # end wrapper



    #Instantiate the graph and produce the bounderies for filters
    @app.callback([
        Output("cyto-graph", "elements"),
        Output("network-filter", "min"),
        Output("network-filter", "max"),
        Output("network-filter", "value")
        ],
        [Input("network-interaction", "value")])
    def make_graph(interaction):
        
        if interaction == "all": # if no interaction is selected, use full graph
            G_cyto = nx.cytoscape_data(G_flat)
            weights = [d["weight"] for u, v, d in G_flat.edges(data=True)]

            # prepare data for network graph
            nodes = G_cyto["elements"]["nodes"]
            edges = G_cyto["elements"]["edges"]
            elements = nodes+edges

            return elements, min(weights), max(weights), np.mean(weights)


        else: # an interaction is selected, select only that interaction
            G_split_flat = cn.graph.flatten_graph(G_split[interaction], weight="score", log=True)
            # Add colors to edges(source node color) for G_split_flat
            for u, v,  in G_split_flat.edges():
                G_split_flat.edges[u,v]["color"] = color_map_nodes[u][0:3]

            G_cyto = nx.cytoscape_data(G_split_flat)
            weights = [d["weight"] for u, v, d in G_split_flat.edges(data=True)]

            # prepare data for network graph
            nodes = G_cyto["elements"]["nodes"]
            edges = G_cyto["elements"]["edges"]
            elements = nodes+edges

            return elements, min(weights), max(weights), np.mean(weights)





    # Change layout of network graph
    @app.callback(Output("cyto-graph", "layout"),
                  [Input("network-layout", "value")])
    def update_network_layout(layout):
        return {
            "name": layout,
            "automate": True,
            "fit": True
        }



    # Choose gene to color nodes by
    @app.callback(
        [Output("node-colors", "children"),
        Output("min-max", "children")],
        [
        Input("gene_dropdown", "value")
        ]
    )
    def calculate_colors(gene):
        if gene is None:
            return [None, ""]
        #get all gene expression values for selected gene
        gene_data = {celltype["data"]["id"]: celltype["data"]["genes"][gene] for celltype in nodes}
        
        min_value = min(gene_data.values())
        max_value = max(gene_data.values())

        #package min max expression information to a list that will be returned
        expression = html.Ul(children=[
            html.Li(f"minimum gene expression: {min_value}"), 
            html.Li(f"maximum gene expression: {max_value}")
            ]
        )
        
        cmap = matplotlib.cm.get_cmap("coolwarm")
        
        color_dict= dict()
        for k, v in gene_data.items():
            color_dict[k] = {"rgb": cmap(v, bytes=True)[0:3], "expression": v}

        color = pd.Series(color_dict)

        return color.to_json(), expression


    # Select visible edges of network graph depending on filter value
    # node color depending on selected gene
    # width of edges
    @app.callback(
        Output("cyto-graph", "stylesheet"),
        [
        Input("network-filter", "value"),
        Input("network-filter", "min"),
        Input("network-filter", "max"),
        Input("node-size", "value"),
        Input("node-colors", "children")
        ]
    )
    def style_network_graph(th, min_weight, max_weight, size, colors):
        
        # create a filter for edges
        filter_style = [{
            "selector": f"edge[weight < {th}]",
            "style": {
                "display": "none"
            }
        },
        {
            "selector": "node",
            "style": {
                'height': f'mapData(centrality, {min(cent)}, {max(cent)}, {size[0]}, {size[1]})',
                'width': f'mapData(centrality, {min(cent)}, {max(cent)}, {size[0]}, {size[1]})'
            }
        }]
        

        # create a color style for nodes based on gene expression
        if isinstance(colors, str):
            colors = pd.read_json(colors, typ="series", convert_dates=False)
            color_style = [{'selector': f'node[id = "{str(index)}"]', 'style': {'background-color': f'rgb{tuple(colors[index]["rgb"])}'}} for index in colors.index]
            filter_style += color_style
        else:
            color_style = {
                "selector": "node",
                "style": {'background-color': 'BFD7B5'}
            }


        # Map edges width to a set min and max value (scale for visibility)    
        edge_style = [{
            "selector": "edge",
            "style": {
                "width": f"mapData(weight, {min_weight}, {max_weight}, 1, 10)"
            }
        }]

        return default_stylesheet + filter_style + edge_style




    # Produce a table of all edge data from tapped edge
    @app.callback([
        Output("edge-info", "children"),
        Output("edge-selection", "columns"),
        Output("edge-selection", "data")
    ],
        [Input("cyto-graph", "tapEdgeData")])
    def update_data(edge):
        import pandas as pd

        info = f"Interactions from source: {edge['source']} to target: {edge['target']}. Weight: {edge['weight']}"

        columns = [{"name": i, "id": i} for i in [
            "interaction", "receptorfamily", "score", "log_score", "weighted_score", "pubmed_id", "target_species", "action"]]

        interactions = pd.DataFrame(edge["interactions"])[
            ["interaction", "receptorfamily", "score", "log_score", "weighted_score", "pubmed_id", "target_species", "action"]]

        interactions.sort_values(by="score", ascending=False, inplace=True)
        records = interactions.to_dict("records")

        return [info, columns, records]




    # Produce ligand and receptor graphs based on tapped node
    @app.callback([
        Output("ligand-graph", "figure"),
        Output("receptor-graph", "figure")
        ],
        [
        Input("cyto-graph", "tapNodeData"),
        Input("filter_l_r", "value")
        ]
    )
    def plot_l_r_expression(node, filter_text):

        if isinstance(node, dict):
            
            ligands = pd.DataFrame.from_dict(node['ligands'], orient="index", columns=["log_scores"])
            ligands = ligands.sort_values("log_scores", axis=0, ascending=True)
            ligands = np.log10(ligands + 1) # Turn ligand score into log10 +1
            if filter_text != "":
                ligands = ligands.filter(like=filter_text, axis = 0)

            ligand_fig = go.Figure(
                data=go.Bar(y=list(ligands.index), x=list(ligands["log_scores"]), orientation="h"),
                layout=go.Layout(title=f"Ligands: {node['name']}", showlegend=False, xaxis={'title': "log(Ligand Score)"}, autosize=True, height=800)
                )

            
            receptors = pd.DataFrame.from_dict(node['receptors'], orient="index", columns=["log_scores"])
            receptors = receptors.sort_values("log_scores", axis=0, ascending=True)
            receptors = np.log10(receptors +1) # Turn receptor scores into log10 +1
            if filter_text != "":
                receptors = receptors.filter(like=filter_text, axis = 0)

            receptors_fig = go.Figure(
                data=go.Bar(y=list(receptors.index), x=list(receptors["log_scores"]), orientation="h"),
                layout=go.Layout(title=f"Receptors: {node['name']}", showlegend=False, xaxis={'title': "log(Receptor Score)"}, autosize=True, height=800)
                )
            
            
        else:
            figure=dict()
        
            
        return [ligand_fig, receptors_fig]





    # Builds a sankey graph based on the tapped node
    @app.callback(
        Output("sankey-graph", "figure"),
        [
            Input("cyto-graph", "tapNodeData"),
            Input("sankey-filter", "value"),
            Input("sankey-toggle", "value")
        ]
    )
    def build_sankey_graph(node, th, score):
        node = node["id"]
        # Find all interactions where node is target or source node
        G_s = nx.MultiDiGraph()
        for n, nbrs in G.adj.items():
            for nbr, edict in nbrs.items():
                if n == node:
                    for e, d in edict.items():
                        if d[score] > th:
                            # append dash after the target node
                            G_s.add_edge(n, " Post " + nbr, **d)
                if nbr == node:
                    for e, d in edict.items():
                        if d[score] > th:
                            # append dash before the source node
                            G_s.add_edge("Pre " + n, nbr, **d)
        
        edges = nx.to_pandas_edgelist(G_s)
        if len(edges) < 1:
            fig = dict()
            return fig
        # add same color scheme as network graph
        for node_s in G_s.nodes():
            if " Post" in node_s:
                original_node = str(node_s).split(sep=" Post")[1]
            elif "Pre " in node_s:
                original_node = str(node_s).split(sep="Pre ")[1]
            else:
                original_node = str(node_s)

            new_color = color_map_nodes[original_node.strip()]
            G_s.nodes[node_s]["color"] = new_color

        nodes = G_s.nodes()

        node_map = {cluster: id for id, cluster in enumerate(list(nodes))}

        sankey = go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(
                    color="black",
                    width=0.5
                ),
                label=list(nodes),
                color=[f'rgb{tuple(d["color"][0:3])}' for n, d in G_s.nodes(data=True)]

            ),
            link=dict(
                source=list(edges["source"].map(node_map)),
                target=list(edges["target"].map(node_map)),
                value=list(edges[score]),
                label=edges["interaction"]
            )
        )

        data = [sankey]

        layout = go.Layout(
            autosize=True,
            title=f"Interactions: {node}",
            font=dict(
                size=font_size
            )
        )

        fig = go.Figure(data=data, layout=layout)

        return fig

    
    
    
    
    # Run server
    app.run_server(mode=mode, debug=False)



if __name__ == "__main__":

    import scConnect as cn
    print("Loading graph file....")
    G = cn.graph.load("assets/graph.yaml")
    print("graph file loaded.")
    G = cn.graph.score(G)
    #G = "test"
    graph(G)
