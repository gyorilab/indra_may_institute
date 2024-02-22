import argparse
import pandas


def run_indra_demo():
	args = get_arguments()
	filename = get_filename(args)
	pandas_df = pandas.read_csv(filename)
	print(pandas_df)


def get_filename(args):
	return args["filename"]


def get_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filename")
	return vars(parser.parse_args())


if __name__ == "__main__":
	run_indra_demo()
