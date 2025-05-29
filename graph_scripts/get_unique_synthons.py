from argparse import ArgumentParser

import gzip
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm


def main():
    parser = ArgumentParser()
    parser.add_argument('--csv_file')
    parser.add_argument('--output_csv')
    args = parser.parse_args()

    substructures = set()

    with gzip.GzipFile(filename=args.csv_file) as file:
        for line in tqdm(file):
            line = line.decode().strip()
            synthon = line.split(',')[2].split('|')[1]
            if synthon.count('Xe') == 1:
                substructures.add(synthon)

    print(len(substructures), 'unique substructures')
    df = pd.DataFrame({list(substructures)[0]: list(substructures)[1:]})
    df.to_csv(args.output_csv, index=False, compression='gzip')

if __name__ == "__main__":
    main()
