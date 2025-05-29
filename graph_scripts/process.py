import argparse
import csv
import gzip
import os
import time

from joblib import Parallel, delayed
from tqdm import tqdm


def main():
    parser = argparse.ArgumentParser()
    # parser.add_argument('--edges_dir')
    parser.add_argument('--edges_file')
    parser.add_argument('--output_dir')
    # parser.add_argument('--atom_mapping', action='store_true')
    parser.add_argument('--prop_file', default='synthon_properties.csv.gz', required=True)
    parser.add_argument('--n_descriptors', type=int, default=5)
    args = parser.parse_args()

    print('Loading properties file')
    prop_dict = {}
    with gzip.GzipFile(filename=args.prop_file) as file:
        for line in tqdm(file):
            d = line.decode().strip().split(',')
            prop_dict[d[0]] = d[1:]

    print('Read properties file')
    print(len, 'synthons in property dictionary')
    process_one_file(args.edges_file, prop_dict, args.output_dir, args.n_descriptors)
    # edges_files = [os.path.join(args.edges_dir, i) for i in os.listdir(args.edges_dir) if 'csv.gz' in i]
    # Parallel(n_jobs=len(edges_files), backend="multiprocessing")(
    #     delayed(process_one_file)(edges_file, prop_dict, args.output_dir)
    #     for edges_file in edges_files
    # )


def process_one_file(edges_file, prop_dict, output_dir, n_descriptors):
    fz = gzip.open(os.path.join(output_dir, os.path.basename(edges_file)), 'wt')
    writer = csv.writer(fz, delimiter=',')

    print('Processing edges file')
    start = time.time()
    with gzip.GzipFile(filename=edges_file) as file:
        for i, line in enumerate(tqdm(file)):
            e = line.decode().strip().split(',')
            # add the properties
            if e[2] not in prop_dict.keys():
                e = e + [None]*n_descriptors
            else:
                e = e + prop_dict[e[2]]
            writer.writerow(e)

    fz.close()
    end = time.time()
    print('Time taken:', end-start)
    print('Written new file', os.path.basename(edges_file))


if __name__ == "__main__":
    main()
