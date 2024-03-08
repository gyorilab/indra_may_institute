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


def create_plotly_graph(res):
	output = res.json()
	G = nx.DiGraph()
	for entry in output:
		G.add_node(entry['source_id'])
		G.add_node(entry['target_id'])
		G.add_edge(entry['source_id'], entry['target_id'])

	# Setting positions at random
	Z = nx.random_geometric_graph(len(G.nodes), 0.125)
	for index, node in enumerate(G.nodes):
		nx.set_node_attributes(G, {node: Z.nodes[index]['pos']}, name='pos')

	# Need to add position....

	edge_x = []
	edge_y = []
	for edge in G.edges():
		x0, y0 = G.nodes[edge[0]]['pos']
		x1, y1 = G.nodes[edge[1]]['pos']
		edge_x.append(x0)
		edge_x.append(x1)
		edge_x.append(None)
		edge_y.append(y0)
		edge_y.append(y1)
		edge_y.append(None)

	edge_trace = go.Scatter(
		x=edge_x, y=edge_y,
		line=dict(width=0.5, color='#888'),
		hoverinfo='none',
		mode='lines')

	node_x = []
	node_y = []
	for node in G.nodes():
		x, y = G.nodes[node]['pos']
		node_x.append(x)
		node_y.append(y)

	node_trace = go.Scatter(
		x=node_x, y=node_y,
		mode='markers',
		hoverinfo='text',
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
		node_text.append('# of connections: ' + str(len(adjacencies[1])))

	node_trace.marker.color = node_adjacencies
	node_trace.text = node_text

	fig = go.Figure(data=[edge_trace, node_trace],
					layout=go.Layout(
						title='<br>Network graph made with Python',
						titlefont_size=16,
						showlegend=False,
						hovermode='closest',
						margin=dict(b=20, l=5, r=5, t=40),
						annotations=[dict(
							text="Python code: <a href='https://plotly.com/ipython-notebooks/network-graphs/'> https://plotly.com/ipython-notebooks/network-graphs/</a>",
							showarrow=False,
							xref="paper", yref="paper",
							x=0.005, y=-0.002)],
						xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
						yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
					)
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
