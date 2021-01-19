
# Dash app to investigate graphs


def graph(G, mode="external", **kwargs):
    """
    G: a multidirectional graph

    kwargs are passed to the Jupyter_Dash.run_server() function. Some usefull arguments are:
        mode: "inline" to run app inside the jupyter nodebook, default is external 
        debug: True or False, Usefull to catch errors during development.
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
    import plotly.io as pio
    import pandas as pd
    import numpy as np
    import json
    import matplotlib
    import matplotlib.pyplot as plt
    pio.templates.default = "plotly_white"

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
        G.edges[u, v, k]["color"] = color_map_nodes[u][0:3]

    # load graph into used formes
    def G_to_flat(G, weight):
        G_flat = cn.graph.flatten_graph(G, weight=weight, log=True)

        # Add colors to edges(source node color) for G_flat
        for u, v, in G_flat.edges():
            G_flat.edges[u, v]["color"] = color_map_nodes[u][0:3]
        return G_flat

    # produce full graph variante to extract metadata
    G_flat =G_to_flat(G, weight="score")
    G_split = cn.graph.split_graph(G)

    # find and sort all found interactions
    interactions = list(G_split.keys())
    interactions.sort()

    G_cyto = nx.cytoscape_data(G_flat)

    # get min and max weight for all edges for flat and normal graph
    #weights = [d["weight"] for u, v, d in G_flat.edges(data=True)]
    scores = [d["score"] for u, v, d in G.edges(
        data=True)]
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
                html.H3(
                    f'Loaded graph with {len(G.nodes())} nodes and {len(G.edges())} edges')
            ])
        ]),

        html.Div(className="network-settings", children=[  # network settings
            html.H2("Network settings", style={"text-align": "center"}),

            html.Label("Interactions"),
            dcc.Dropdown(
                id="network-interaction",
                options=[{'label': "all interactions", 'value': "all"}] +
                [{'label': interaction, 'value': interaction}
                    for interaction in interactions],
                value="all"
            ),
            # select if only significant ligands and receptors should be shown
            html.Label("Graph weight:"),
            dcc.RadioItems(id="weight-select", 
                options=[
                    {"label": "Score", "value": "score"},
                    {"label": "Log score", "value": "log_score"},
                    {"label": "Specificity", "value": "specificity"},
                    {"label": "Importance", "value": "importance"}],
                value="importance",
                labelStyle={
                    'display': 'block',
                    "margin-left": "50px"
                },
                style={"padding": "10px", "margin": "auto"}
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

            html.Label("Weight Filter", style={
                       "paddingBottom": 500, "paddingTop": 500}),
            dcc.Slider(  # min, max and value are set dynamically via a callback
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
                value=[10, 50],
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
            
            # Store node colors "hidden" for gene expresison
            html.Div(id="node-colors", style={"display": "none"}, children=[
                ""
            ]),

            html.Div(id="min-max", children=[

            ]),
            # Click to download image of network graph
            html.Button(children="Download current view",  id="download-network-graph", style={"margin":"10px"}) 
            
        ]),  # end network settings
        html.Div(id="network-graph", className="network-graph", children=[  # network graph
            html.H2("Network graph", style={"text-align": "center"}),
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
            html.H2("Sankey Settings", style={"text-align": "center"}),

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
                {"label": "Log score", "value": "log_score"},
                {"label": "Specificity", "value": "specificity"},
                {"label": "Importance", "value": "importance"}
            ], value="importance",
            labelStyle={"display": "block"})

        ]),  # end network settings

        html.Div(className="sankey", id="sankey", children=[  # sankey graph
            html.H2("Sankey graph", style={"text-align": "center"}),
            dcc.Graph(id="sankey-graph")
        ]),  # end sankey graph

        html.Div(className="interaction-list", children=[  # interaction list

            html.Div(id="selection", children=[
                html.H2("Interactions", style={"text-align": "center"}),
                html.H3(id="edge-info", style={"text-align": "center"}),
            
                dcc.Graph(id="interaction-scatter"),

                html.Div(id="interaction-selection", style={"display": "none"}, children=[
                ""
                ])
            ]),
            html.Div(children=[
                dash_table.DataTable(
                    id="edge-selection",
                    page_size=20,
                    style_table={
                        "overflowX": "scroll",
                        "overflowY": "scroll",
                        "height": "50vh",
                        "width": "95%"
                    },

                    style_cell_conditional=[
                        {
                            "if": {"column_id": "interaction"},
                            "textAlign": "left"
                        },
                        {
                            "if": {"column_id": "receptorfamily"},
                            "textAlign": "left"
                        },
                        {
                            "if": {"column_id": "pubmedid"},
                            "textAlign": "left"
                        }
                    ],

                    style_header={
                        "fontWeight": "bold",
                        "maxWidth": "200px",
                        "minWidth": "70px"
                    },

                    style_data={
                        "maxWidth": "200px",
                        "minWidth": "70px",
                        "textOverflow": "ellipsis"
                    },
                    sort_action="native",

                    fixed_rows={'headers': True, 'data': 0}
                )
            ])
        ]),  # end interaction list

        html.Div(className="L-R-scores", children=[  # ligand and receptor lists
            html.H2("Ligand and receptors", style={"text-align": "center"}),
            html.Div(children=[
                html.H3(id="selected-node", style={"text-align": "center"}, children=["Select a node in the notwork graph"]),
                html.Label("Search for ligands and receptors:", style={"margin-right": "10px"}),
                dcc.Input(id="filter_l_r", type="search", value="", placeholder="Search")
            ]),
           
           
            dcc.Tabs([
                dcc.Tab(label="Ligands", children=[
                    dcc.Graph(
                        id="ligand-graph",
                        config=dict(
                            autosizable=True,
                            responsive=True)
                    ),
                    dash_table.DataTable(
                    id="ligand-table",
                    page_size=20,
                    style_table={
                        "overflowX": "scroll",
                        "overflowY": "scroll",
                        "height": "50vh",
                        "width": "95%"
                    },

                    style_cell_conditional=[
                        {
                            "if": {"column_id": "Ligand"},
                            "textAlign": "left"
                        }
                    ],

                    style_header={
                        "fontWeight": "bold",
                        "maxWidth": "200px",
                        "minWidth": "70px"
                    },

                    style_data={
                        "maxWidth": "200px",
                        "minWidth": "70px",
                        "textOverflow": "ellipsis"
                    },
                    sort_action="native",

                    fixed_rows={'headers': True, 'data': 0})
                ]),
                dcc.Tab(label="Receptors", children=[
                    dcc.Graph(
                        id="receptor-graph",
                        config=dict(
                            autosizable=True,
                            responsive=True)
                    ),
                    dash_table.DataTable(
                    id="receptor-table",
                    page_size=20,
                    style_table={
                        "overflowX": "scroll",
                        "overflowY": "scroll",
                        "height": "50vh",
                        "width": "95%"
                    },

                    style_cell_conditional=[
                        {
                            "if": {"column_id": "Receptor"},
                            "textAlign": "left"
                        }
                    ],

                    style_header={
                        "fontWeight": "bold",
                        "maxWidth": "200px",
                        "minWidth": "70px"
                    },

                    style_data={
                        "maxWidth": "200px",
                        "minWidth": "70px",
                        "textOverflow": "ellipsis"
                    },
                    sort_action="native",

                    fixed_rows={'headers': True, 'data': 0})
                ])
            ])
        ])  # end ligand receptor list
    ])  # end wrapper

    # Instantiate the graph and produce the bounderies for filters
    @app.callback([
        Output("cyto-graph", "elements"),
        Output("network-filter", "min"),
        Output("network-filter", "max"),
        Output("network-filter", "value")
        ],
        [
        Input("network-interaction", "value"),
        Input("weight-select", "value")])
    def make_graph(interaction, score):
        G_flat = G_to_flat(G, score)
            

        if interaction == "all":  # if no interaction is selected, use full graph
            G_cyto = nx.cytoscape_data(G_flat)
            weights = [d["weight"] for u, v, d in G_flat.edges(data=True)]

            # prepare data for network graph
            nodes = G_cyto["elements"]["nodes"]
            edges = G_cyto["elements"]["edges"]
            elements = nodes+edges

            return elements, min(weights), max(weights), np.mean(weights)

        else:  # an interaction is selected, select only that interaction
            G_split = cn.graph.split_graph(G)
            G_split_flat = G_to_flat(G_split[interaction], score)
            G_cyto = nx.cytoscape_data(G_split_flat)
            weights = [d["weight"]
                       for u, v, d in G_split_flat.edges(data=True)]

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
        # get all gene expression values for selected gene
        gene_data = {celltype["data"]["id"]: celltype["data"]
                     ["genes"][gene] for celltype in nodes}

        min_value = min(gene_data.values())
        max_value = max(gene_data.values())

        # package min max expression information to a list that will be returned
        expression = html.Ul(children=[
            html.Li(f"minimum gene expression: {min_value}"),
            html.Li(f"maximum gene expression: {max_value}")
        ]
        )

        cmap = matplotlib.cm.get_cmap("coolwarm")

        color_dict = dict()
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
            color_style = [{'selector': f'node[id = "{str(index)}"]', 'style': {
                'background-color': f'rgb{tuple(colors[index]["rgb"])}'}} for index in colors.index]
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

    # download an image of current network graph view
    @app.callback(
        Output("cyto-graph", "generateImage"), 
        Input("download-network-graph", "n_clicks"))
    def download_networkgraph_image(get_request):
        
        if get_request == None:
            return dict()

        return {
            "type": "svg",
            "action": "download"
        }


    # Produce a table of all edge data from tapped edge
    @app.callback([
        Output("edge-info", "children"),
        Output("edge-selection", "columns"),
        Output("edge-selection", "data")
    ],
        [Input("cyto-graph", "tapEdgeData"),
        Input("interaction-selection", "children")])
    def update_data(edge, selection):
        import pandas as pd
        import json

        # check if an edge has really been clicked, return default otherwise
        if edge is None:
            return ["", None, None]

        info = f"Interactions from {edge['source']} to {edge['target']}."

        # map visible names for columns with columns in edge[interaction]
        columns = [
            {
                "name": "Interaction",
                "id": "interaction"
            },
            {
                "name": "Receptor Family",
                "id": "receptorfamily"
            },
            {
                "name": "Score",
                "id": "score"
            },
            {
                "name": "Log10(score)",
                "id": "log_score"
            },
            {
                "name": "Specificity",
                "id": "specificity"
            },
            {
                "name": "Importance",
                "id": "importance"
            },
            {
                "name": "Ligand z-score",
                "id": "ligand_zscore"
            },
            {
                "name": "Ligand p-value",
                "id": "ligand_pval"
            },
            {
                "name": "Receptor z-score",
                "id": "receptor_zscore"
            },
            {
                "name": "Receptor p-value",
                "id": "receptor_pval"
            },
            {
                "name": "PubMed ID",
                "id": "pubmedid"
            }
        ]

        interactions = pd.DataFrame(edge["interactions"])[
            ["interaction", "receptorfamily", "score", "log_score", "specificity", "importance", "ligand_zscore",
             "ligand_pval", "receptor_zscore", "receptor_pval", "pubmedid"]]

        # Sort values based on score
        interactions.sort_values(by="score", ascending=False, inplace=True)

        # round values for scores to two decimals
        interactions[[
            "score",
            "log_score",
            "specificity",
            "importance",
            "ligand_zscore",
            "receptor_zscore"]] = interactions[[
                "score",
                "log_score",
                "specificity",
                "importance",
                "ligand_zscore",
                "receptor_zscore"]].round(decimals=2)

        interactions[["ligand_pval", "receptor_pval"]] = interactions[[
            "ligand_pval", "receptor_pval"]].round(decimals=4)

        # if selection from interaction graph, filter dataframe
        if selection != "":
            selection = json.loads(selection)
            interactions = interactions.loc[interactions["interaction"].isin(selection)]

        records = interactions.to_dict("records")

        return [info, columns, records]


    @app.callback([
        Output("interaction-scatter", "figure")
    ],
        [Input("cyto-graph", "tapEdgeData")])
    def interaction_scatter_plot(edge):
        import plotly.express as px

        fig = go.Figure()
        if not isinstance(edge, dict):
            return [fig, ]

        interactions = pd.DataFrame(edge["interactions"])[
                    ["interaction", "receptorfamily", "score", "log_score", "ligand_zscore",
                    "ligand_pval", "receptor_zscore", "receptor_pval", "specificity", "importance", "pubmedid"]]
        
        # add 10% to the min and max value to not clip the datapoint
        range_x = (-max(interactions["log_score"])*0.1, max(interactions["log_score"])*1.1)
        range_y = (-max(interactions["specificity"])*0.1, max(interactions["specificity"])*1.1)
        #interactions["specificity"] = np.log10( interactions["specificity"])

        fig = px.scatter(interactions, 
                x="log_score",
                range_x=range_x,
                y="specificity",
                range_y=range_y,
                color="importance",
                hover_name="interaction",
                hover_data=["ligand_pval", "receptor_pval", "score","specificity", "receptorfamily"],
                color_continuous_scale=px.colors.sequential.Viridis_r,
                labels={
                    "ligand_zscore": "Ligand Z-score",
                    "receptor_zscore": "Receptor Z-score",
                    "log_score": "log(Interaction score)",
                    "score": "Interaction score",
                    "specificity": "Specificity",
                    "importance": "Importance",
                    "receptorfamily": "Receptor family",
                    "pubmedid": "PubMed ID",
                    "ligand_pval": "Ligand p-value",
                    "receptor_pval": "Receptor p-value"
                    }
                )
        return [fig,]

    @app.callback(
        Output("interaction-selection", "children"),
        [
        Input("interaction-scatter", "selectedData")
        ]
        )
    def interaction_select(selected_data):
        import json
        if isinstance(selected_data, dict):
            interactions = [point["hovertext"] for point in selected_data["points"]]
        else:
            return ""
        return json.dumps(interactions)



    # Produce ligand and receptor graphs based on tapped node

    @app.callback([
        Output("ligand-graph", "figure"),
        Output("receptor-graph", "figure"),
        Output("selected-node", "children")
    ],
        [
        Input("cyto-graph", "tapNodeData"),
        Input("filter_l_r", "value")
    ]
    )
    def plot_l_r_expression(node, filter_text):

        # set output variables to empty figures
        ligand_fig = go.Figure()
        receptor_fig = go.Figure()
        node_id = "Select a node in the network graph"

        if isinstance(node, dict):
            import plotly.express as px

            node_id = node["id"]
            
            ligands_score = pd.DataFrame.from_dict(node["ligands_score"], orient="index", columns=["Score"])
            ligands_zscore = np.log2(pd.DataFrame.from_dict(node["ligands_zscore"], orient="index", columns=["Z-score"]))
            ligands_corr_pval = pd.DataFrame.from_dict(node["ligands_corr_pval"], orient="index", columns=["p-value"])
            ligands_merge = ligands_score.merge(ligands_zscore, how="left", left_index=True, right_index=True)
            ligands_merge = ligands_merge.merge(ligands_corr_pval, how="left", left_index=True, right_index=True)
            ligands_merge["log(score + 1)"] = np.log10(ligands_merge["Score"]+1)
            ligands_merge["Significant"] = [True if p_val < 0.05 else False for p_val in ligands_merge["p-value"]]
            ligands_merge["-log(p-value)"] = -np.log10(ligands_merge["p-value"])
            
            if filter_text != "":
                ligands_merge = ligands_merge.filter(like=filter_text, axis=0)

            ligand_fig = px.scatter(
                ligands_merge, 
                x="log(score + 1)", 
                y="-log(p-value)",
                color="Significant",
                hover_name=ligands_merge.index,
                hover_data=["Score", "Z-score", "p-value"])

            receptors_score = pd.DataFrame.from_dict(node["receptors_score"], orient="index", columns=["Score"])
            receptors_zscore = np.log2(pd.DataFrame.from_dict(node["receptors_zscore"], orient="index", columns=["Z-score"]))
            receptors_corr_pval = pd.DataFrame.from_dict(node["receptors_corr_pval"], orient="index", columns=["p-value"])
            receptors_merge = receptors_score.merge(receptors_zscore, how="left", left_index=True, right_index=True)
            receptors_merge = receptors_merge.merge(receptors_corr_pval, how="left", left_index=True, right_index=True)
            receptors_merge["log(score + 1)"] = np.log10(receptors_merge["Score"]+1)
            receptors_merge["Significant"] = [True if p_val < 0.05 else False for p_val in receptors_merge["p-value"]]
            receptors_merge["-log(p-value)"] = -np.log10(receptors_merge["p-value"])
            
            if filter_text != "":
                receptors_merge = receptors_merge.filter(like=filter_text, axis=0)

            receptor_fig = px.scatter(
                receptors_merge, 
                x="log(score + 1)", 
                y="-log(p-value)",
                color="Significant",
                hover_name=receptors_merge.index,
                hover_data=["Score", "Z-score", "p-value"])
       
        return [ligand_fig, receptor_fig, node_id]

    # Builds a sankey graph based on the tapped node (store in global G_s)
    G_s = nx.MultiDiGraph() #variable holding sankey graph
    @app.callback([
        Output("sankey-filter", "min"),
        Output("sankey-filter", "max"),
        Output("sankey-filter", "value")
        ],
        [
        Input("cyto-graph", "tapNodeData"),
        Input("sankey-toggle", "value")])
    def build_sankey_graph(node, score):
        import numpy as np
        # If no node has been selected, dont try to build graph
        if node is None:
            return (0 ,0 ,0)

        node = node["id"]
        # Find all interactions where node is target or source node
        nonlocal G_s
        G_s = nx.MultiDiGraph() # reset content
        weight = list() # list to store all weights (used to set min and max for the filter)
        for n, nbrs in G.adj.items(): # graph has been modified by network graph before
            for nbr, edict in nbrs.items():
                if n == node:
                    for e, d in edict.items():
                        G_s.add_edge(n, " Post " + nbr, **d)
                        weight.append(d[score])
                if nbr == node:
                    for e, d in edict.items():
                        G_s.add_edge("Pre " + n, nbr, **d)
                        weight.append(d[score])

        if len(weight) == 0:
            weight = [0,1]
        if score == "specificity": 
            # set default start value to specificity value for ligand and receptor 
            # p-value of (0.05 and 0.05)/2 = 1.3
            return (min(weight), max(weight), 1.3)
        return (min(weight), max(weight), np.mean(weight))
    
    @app.callback(
        Output("sankey-graph", "figure"),
        [Input("sankey-filter", "value"),
        Input("sankey-toggle", "value"),
        Input("cyto-graph", "tapNodeData")])
      
    def filter_sankey_graph(th, score, node):

        if node:
            node = node["id"]

        _G_s = nx.MultiDiGraph()
        for u, v, n, d in G_s.edges(data=True, keys=True):
            if d[score] > th:
                _G_s.add_edge(u, v, n, **d)
        _G_s.add_nodes_from(G_s.nodes(data=True))
        
        edges = nx.to_pandas_edgelist(_G_s)
        if len(edges) < 1:
            fig = dict()
            return fig
        # add same color scheme as network graph
        for node_s in _G_s.nodes():
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
                color=[f'rgb{tuple(d["color"][0:3])}' for n,
                       d in G_s.nodes(data=True)]

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

    @app.callback([Output("ligand-table", "columns"), Output("ligand-table", "data")], [Input("ligand-graph", "figure"), Input("ligand-graph", "selectedData")])
    def select_ligands(figure, selected):
        import json
        ligands = []
        score= []
        zscore = []
        pval = []
        
        

        for group in figure["data"]:
            for ligand in group["hovertext"]:
                ligands.append(ligand)
            for data in group["customdata"]:
                score.append(data[0])
                zscore.append(data[1])
                pval.append(data[2])

        
        df = pd.DataFrame({"Ligand": ligands, "Score": score, "Z-score": zscore, "P-value": pval})
        df.index = df["Ligand"]
        df.sort_values(by="Score", ascending=False, inplace=True)

        if isinstance(selected, dict):
            filt=[]
            for point in selected["points"]:
                filt.append(point["hovertext"])
            df = df.loc[filt]

        
        columns = [
            {"name": "Ligand", "id": "Ligand"},
            {"name": "Score", "id": "Score"},
            {"name": "Z-score", "id": "Z-score"},
            {"name": "P-value", "id": "P-value"}]

        data = df.to_dict("records")

        return columns, data

    @app.callback([Output("receptor-table", "columns"), Output("receptor-table", "data")], [Input("receptor-graph", "figure"), Input("receptor-graph", "selectedData")])
    def select_ligands(figure, selected):
        import json
        receptors = []
        score= []
        zscore = []
        pval = []
        
        

        for group in figure["data"]:
            for receptor in group["hovertext"]:
                receptors.append(receptor)
            for data in group["customdata"]:
                score.append(data[0])
                zscore.append(data[1])
                pval.append(data[2])

        
        df = pd.DataFrame({"Receptor": receptors, "Score": score, "Z-score": zscore, "P-value": pval})
        df.index = df["Receptor"]
        df.sort_values(by="Score", ascending=False, inplace=True)


        if isinstance(selected, dict):
            filt=[]
            for point in selected["points"]:
                filt.append(point["hovertext"])
            df = df.loc[filt]

        
        columns = [
            {"name": "Receptor", "id": "Receptor"},
            {"name": "Score", "id": "Score"},
            {"name": "Z-score", "id": "Z-score"},
            {"name": "P-value", "id": "P-value"}]

        data = df.to_dict("records")

        return columns, data

    # Run server
    app.run_server(**kwargs)


if __name__ == "__main__":

    import scConnect as cn
    print("Loading graph file....")
    G = cn.graph.load("assets/graph.yaml")
    print("graph file loaded.")
    G = cn.graph.score(G)
    #G = "test"
    graph(G)
