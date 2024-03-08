import argparse
import pandas
import requests
from protmapper import uniprot_client


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
