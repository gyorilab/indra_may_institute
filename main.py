import argparse
import pandas
import requests
from protmapper import uniprot_client
import networkx as nx
import plotly.graph_objects as go

DATASET_NAME = "SMALL_MOLECULE_INTERVENTION"
DATASET_VALUES = ["BREAST_CANCER", "SMALL_MOLECULE_INTERVENTION"]
LABELS = ["DMSO-DbET6"]
P_VALUE_CONSTRAINT = 0.05


def validate_settings():
    if DATASET_NAME not in DATASET_VALUES:
        raise Exception("ADJUST YOUR DATASET NAME")


def run_indra_demo():
    args = get_arguments()
    filename = get_filename(args)
    res = query_indra(filename)
    print(res.json())


def construct_networkx_graph(res):
    output = res.json()
    G = nx.DiGraph()
    for entry in output:
        if entry['data']['stmt_type'] == 'Complex':
            G.add_node(entry['source_id'])
            G.add_node(entry['target_id'])
            if G.has_edge(entry['target_id'], entry['source_id']):
                if (entry['data']['evidence_count']
                        <= G.get_edge_data(entry['target_id'], entry['source_id'])['evidence']):
                    # Accounts for bidirectional edges for now
                    continue
            G.add_edge(
                entry['source_id'],
                entry['target_id'],
                evidence=entry['data']['evidence_count'],
                belief=entry['data']['belief'],
                type=entry['data']['stmt_type']
            )

    # TODO: Set positions based on gene set clusters
    Z = nx.random_geometric_graph(len(G.nodes), 0.125)
    for index, node in enumerate(G.nodes):
        nx.set_node_attributes(G, {node: Z.nodes[index]['pos']}, name='pos')
    return G


def create_plotly_graph(G):
    mnode_x, mnode_y, mnode_txt = [], [], []
    arrow_list = []
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']

        mnode_x.extend([(x0 + x1) / 2])
        mnode_y.extend([(y0 + y1) / 2])
        mnode_txt.extend([f'HGNC:{edge[0]}->HGNC:{edge[1]} evidence count: {G.edges[edge]['evidence']}'])

        arrow = go.layout.Annotation(dict(
            x=x0,
            y=y0,
            xref="x", yref="y",
            showarrow=True,
            axref="x", ayref='y',
            ax=x1,
            ay=y1,
            arrowhead=3,
            arrowwidth=1.5,
            arrowcolor='rgb(255,51,0)', )
        )

        arrow_list.append(arrow)

    mnode_trace = go.Scatter(
        x=mnode_x, y=mnode_y,
        mode="markers",
        showlegend=False,
        hovertext=mnode_txt,
        hovertemplate="Edge %{hovertext}<extra></extra>",
        marker=go.scatter.Marker(opacity=0)
    )

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        hoverinfo='text',
        text=['HGNC:' + node for node in list(G.nodes())],
        textposition="bottom center",
        marker=dict(
            showscale=True,
            # colorscale options
            # 'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            # 'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            # 'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))

    node_adjacencies = []
    node_text = []
    for node, adjacencies in enumerate(G.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append('# of connections into node: ' + str(len(adjacencies[1])))

    node_trace.marker.color = node_adjacencies
    node_trace.hovertext = node_text

    fig = go.Figure(data=[mnode_trace, node_trace],
                    layout=go.Layout(
                        title='<br>Network graph made with Python',
                        font=dict(
                            family="Courier New, monospace",
                            size=18,
                            color="RebeccaPurple"
                        ),
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    fig.update_layout(annotations=arrow_list)
    fig.show()


def search_uniprot(mnemonic_id):
    protein_id = uniprot_client.get_id_from_mnemonic(mnemonic_id)
    if protein_id:
        return uniprot_client.get_hgnc_id(protein_id)
    else:
        return None


def construct_df(filename):
    pandas_df = pandas.read_csv(filename)
    if DATASET_NAME == "BREAST_CANCER":
        pandas_df[['sp', 'Protein', 'Gene']] = pandas_df['Protein'].str.split('|', expand=True)
        pandas_df = pandas_df.loc[pandas_df['adj.pvalue'] < P_VALUE_CONSTRAINT]
        pandas_df['HGNC'] = pandas_df['Protein'].apply(lambda protein_id: uniprot_client.get_hgnc_id(protein_id))
    elif DATASET_NAME == "SMALL_MOLECULE_INTERVENTION":
        pandas_df = pandas_df.loc[pandas_df['adj.pvalue'] < P_VALUE_CONSTRAINT]
        pandas_df = pandas_df.loc[pandas_df['issue'].isnull()]
        pandas_df = pandas_df.loc[pandas_df['Label'].isin(LABELS)]
        pandas_df['HGNC'] = pandas_df['Protein'].apply(
            lambda mnemonic_id: search_uniprot(mnemonic_id)
        )
        pandas_df = pandas_df.loc[pandas_df['HGNC'].notnull()]
    return pandas_df


def query_indra(filename):
    pandas_df = construct_df(filename)
    groundings = create_groundings(pandas_df)
    res = requests.post(
        'https://discovery.indra.bio/api/indra_subnetwork_relations',
        json={'nodes': list(groundings.values())}
    )
    return res


def create_groundings(df):
    hgnc_ids = set()
    groundings = {}
    for index, row in df.iterrows():
        hgnc_ids.add(row['HGNC'])
    for index, entry in enumerate(hgnc_ids):
        groundings[index] = 'HGNC', entry
    return groundings


def get_filename(args):
    return args["filename"]


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename")
    return vars(parser.parse_args())


if __name__ == "__main__":
    run_indra_demo()
