import argparse
import pandas
import requests
from protmapper import uniprot_client


def run_indra_demo():
	args = get_arguments()
	filename = get_filename(args)
	return query_indra(filename)


def construct_df(filename):
	pandas_df = pandas.read_csv(filename)
	pandas_df[['sp', 'Protein', 'Gene']] = pandas_df['Protein'].str.split('|', expand=True)
	pandas_df = pandas_df.loc[pandas_df['adj.pvalue'] < 0.01]
	pandas_df['HGNC'] = pandas_df['Protein'].apply(lambda protein_id: uniprot_client.get_hgnc_id(protein_id))
	return pandas_df


def query_indra(filename):
	pandas_df = construct_df(filename)
	print(len(pandas_df))
	groundings = create_groundings(pandas_df)
	res = requests.post(
		'https://discovery.indra.bio/api/indra_subnetwork_relations',
		json={'nodes': list(groundings.values())}
	)
	print(res.json())
	return res


def create_groundings(df):
	groundings = {}
	for index, row in df.iterrows():
		groundings[index] = 'HGNC', row['HGNC']
	return groundings


def get_filename(args):
	return args["filename"]


def get_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filename")
	return vars(parser.parse_args())


if __name__ == "__main__":
	run_indra_demo()
