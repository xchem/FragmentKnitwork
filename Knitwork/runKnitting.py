import os
import shutil
import time
from argparse import ArgumentParser

from FragmentKnitwork.Knitwork.queries import impure_expansion, pure_expansion
from FragmentKnitwork.utils import knitworkConfig as config
import FragmentKnitwork.utils.Config as base_config

from FragmentKnitwork.utils.quilterUtils import split_fragment_pair_string
from FragmentKnitwork.utils.knitworkUtils import calc_pharm_fp, load_sigFactory
from FragmentKnitwork.utils.utils import dump_json, load_json
from joblib import Parallel, delayed
from rdkit import Chem

import logging
logger = logging.getLogger('FragmentKnitwork')

def pure_merge(substructure_pair: list, working_dir: str, idx=None, total=None, write_to_file=True, limit=None):
    """
    Run the pure merging query (incorporate an exact substructure into a merge).

    :param substructure_pair: pair of substructures to use, first is the initial seed node
    :param working_dir: working directory to save results
    :param idx: index of the sub pair being run
    :param total: total number of sub pairs being run
    :param write_to_file: whether to save results
    :param return_results: whether to return results
    :return: the query data
    """
    output_fname = os.path.join(working_dir, f"{substructure_pair}_pure_merge.json")
    if os.path.exists(output_fname):
        logger.warning(f'File already exists for this query: {substructure_pair}')
        return load_json(output_fname)

    subnode, synthon = substructure_pair.split('*')[0], substructure_pair.split('*')[1]
    logger.title(f'Running substructure pair {substructure_pair}: {idx + 1}/{total}')

    start = time.time()

    expansions, cmpd_ids = pure_expansion(subnode, synthon, limit=limit)
    logger.var(substructure_pair, f"{len(expansions)} expansions")
    substructure_pair_data = {'expansions': expansions,
                              'cmpd_ids': cmpd_ids}

    end = time.time()
    time_taken = round(end - start, 2)
    time_data = {'time': time_taken,
                 'limit': limit}

    if write_to_file:
        logger.writing(output_fname)
        dump_json(substructure_pair_data, output_fname)
        logger.writing(os.path.join(working_dir, f"{substructure_pair}_pure_merge_timing.json"))
        dump_json(time_data, os.path.join(working_dir, f"{substructure_pair}_pure_merge_timing.json"))

    return substructure_pair_data


def impure_merge(substructure_pair: list, working_dir: str, desc, idx=None, total=None, write_to_file=True, limit=None):
    """
    Run the impure (bioisosteric) merging query, where a substructure can be replaced by an analogue.

    :param substructure_pair: pair of substructures to use, first is the initial seed node
    :param working_dir: working directory to save results
    :param idx: index of the sub pair being run
    :param total: total number of sub pairs being run
    :param write_to_file: whether to save results
    :param return_results: whether to return results
    :return: the query data
    """
    output_fname = os.path.join(working_dir, f"{substructure_pair}_{desc}_impure_merge.json")
    if os.path.exists(output_fname):
        logger.warning(f'File already exists for this query: {substructure_pair}')
        return load_json(output_fname)

    if desc == 'prop_pharmfp':
        sigFactory = load_sigFactory()

    subnode, synthon = substructure_pair.split('*')[0], substructure_pair.split('*')[1]
    logger.title(f'Running substructure pair {substructure_pair}: {idx + 1}/{total}')

    # get the descriptor for the query substructure to calculate similarity against
    # TODO: right now the code only supports the pre-computed prop_pharmfp descriptor
    if desc == 'prop_pharmfp':
        vector = calc_pharm_fp(Chem.MolFromSmiles(synthon), sigFactory, asStr=False)

    start = time.time()

    expansions, repl_synthons, similarities, cmpd_ids = impure_expansion(subnode, vector, synthon, limit=limit)
    logger.var(substructure_pair, f"{len(expansions)} expansions")
    substructure_pair_data = {'expansions': expansions,
                              'used_synthons': repl_synthons,  # the similar substructures actually incorporated
                              'similarities': similarities,
                              'cmpd_ids': cmpd_ids}

    end = time.time()
    time_taken = round(end-start, 2)
    time_data = {'time': time_taken,
                 'limit': limit}

    if write_to_file:
        logger.writing(output_fname)
        dump_json(substructure_pair_data, output_fname)
        logger.writing(os.path.join(working_dir, f"{substructure_pair}_pure_merge_timing.json"))
        dump_json(time_data, os.path.join(working_dir, f"{substructure_pair}_{desc}_impure_merge_timing.json"))

    return substructure_pair_data


def runKnitting(substructure_pair_file, n_parallel, target, working_dir, output_dir, limit, substructure_dir, descriptor='prop_pharmfp',
                pure_search=False, prolif_prioritization=True, max_prioritize=None, r_group_search=None, r_group_data_file=None,
                equiv_synthon_data_file=None, fragalysis_dir=base_config.FRAGALYSIS_DATA_DIR):
    """

    :param substructure_pair_file: the substructure pair file generated by runEnumeration.py (containing the substructure pairs to run querying on)
    :param n_parallel: number of queries to run in parallel
    :param target: name of target (for Fragalysis-formatted files)
    :param working_dir: to save intermed files
    :param output_dir: to save final results files
    :param limit: limit to the number of molecules to retrieve from the database
    :param descriptor: only support prop_pharmfp right now (for impure/bioisosteric merging)
    :param substructure_dir: the substructure dir created by runEnumeration.py
    :param pure_search: whether to run a pure merge search (not impure/bioisosteric)
    :param prolif_prioritization: whether to prioritize bioisosteric results according to whether substructures may make an interaction
    :param max_prioritize: the max number of compounds to prioritize
    :param r_group_search: whether to run r group expansion on retrieved compounds (will be added as extra compounds in lsit)
    :param r_group_data_file: where r group data is saved (generated by runEnumeration.py) so we know which r groups to query with
    :param equiv_synthon_data_file: used for r group expanding
    :return:
    """
    if prolif_prioritization:
        from FragmentKnitwork.Knitwork.prioritize import prioritize_data
        logger.info('Running prioritization using ProLIF...')
        prolif_dir = os.path.join(output_dir, 'prolif_prioritized')
        if not os.path.exists(prolif_dir):
            os.mkdir(prolif_dir)
        prolif_output, prolif_working = os.path.join(prolif_dir, 'output'), os.path.join(prolif_dir, 'working')
        if not os.path.exists(prolif_output):
            os.mkdir(prolif_output)
        if not os.path.exists(prolif_working):
            os.mkdir(prolif_working)

    if r_group_search:
        r_group_dir = os.path.join(output_dir, 'r_group_expanded')
        if not os.path.exists(r_group_dir):
            os.mkdir(r_group_dir)

    data = load_json(substructure_pair_file)

    # get the unique queries from the enumerated pair data (as the same substructure pairs may exist for multiple
    # fragment pairs
    total_queries = 0
    substructure_pairs = set()
    for fragment_pair in data:
        sub_pairs = data[fragment_pair]
        for _sub_pair in sub_pairs:
            total_queries += 1
            sub_pair = '*'.join(_sub_pair)
            substructure_pairs.add(sub_pair)

    substructure_pairs = list(substructure_pairs)
    num_pairs = len(substructure_pairs)
    idxs = list(range(num_pairs))

    logger.var('#total queries', total_queries)
    logger.var('#queries removing redundancy', num_pairs)

    if pure_search:
        results = Parallel(n_jobs=n_parallel, backend="multiprocessing")(
            delayed(pure_merge)(sub_pair, working_dir, idx, num_pairs, True, limit)
            for sub_pair, idx in zip(substructure_pairs, idxs)
        )
    else:  # if impure merge search
        results = Parallel(n_jobs=n_parallel, backend="multiprocessing")(
            delayed(impure_merge)(sub_pair, working_dir, descriptor, idx, num_pairs, True, limit=limit)
            for sub_pair, idx in zip(substructure_pairs, idxs)
        )

    logger.success('Queries run')
    logger.header('Processing data')
    # save the data into individual files for each FRAGMENT with the substructure pair as keys (joined by an asterisk)
    # and the data as the values
    output_data = {}
    for sub_pair, res in zip(substructure_pairs, results):
        output_data[sub_pair] = res

    n = len(data)

    for i,fragment_pair in enumerate(data):
        fragment_data = {}
        logger.header(f"{i}/{n}")
        logger.var('fragment_pair', fragment_pair)
        if pure_search:
            output_fname = os.path.join(output_dir, f"{fragment_pair}_pure_merge.json")
        else:
            output_fname = os.path.join(output_dir, f"{fragment_pair}_{descriptor}_impure_merge.json")
        fragment_sub_pairs = data[fragment_pair]
        for _sub_pair in fragment_sub_pairs:
            sub_pair = '*'.join(_sub_pair)
            fragment_data[sub_pair] = output_data[sub_pair]
        dump_json(fragment_data, output_fname)

        if prolif_prioritization:

            # second_fragment = fragment_pair.split('-')[1]
            second_fragment = split_fragment_pair_string(fragment_pair)[1]

            logger.header('Prolif prioritization starting')
            try:
            	priori_data = prioritize_data(fragment_data, second_fragment, target, substructure_dir,
                          	prolif_working, n_parallel, os.path.join(prolif_output, f"{fragment_pair}_{descriptor}_impure_merge.json"),
                          	max_prioritize, fragalysis_dir=fragalysis_dir)
            	logger.success('Prolif prioritization run')
            except Exception as e:
                logger.error(f"prolif_priotization error: {e}")
                logger.error("Skipping fragment pair")
                continue

        if r_group_search:
            logger.header('R group expansion starting')
            from FragmentKnitwork.Knitwork.RGroupExpansion import write_r_group_expansions
            r_group_data = load_json(r_group_data_file)
            equiv_synthon_data = load_json(equiv_synthon_data_file)

            if prolif_prioritization:
                fragment_data = priori_data

            write_r_group_expansions(fragment_data, fragment_pair, r_group_data, equiv_synthon_data,
                                     os.path.join(r_group_dir, f"{fragment_pair}_{descriptor}_impure_merge.json"))
            logger.success('R group expansion done')

    if prolif_prioritization:
        shutil.rmtree(prolif_working)

    logger.writing(f'{output_dir}/knitwork.tgz')
    if pure_search:
        os.system(f'cd {output_dir}; tar -czf knitwork_pure.tgz *_merge.json')
    else:
        os.system(f'cd {output_dir}; tar -czf knitwork_impure.tgz *_merge.json')

    logger.success('Files processed')



def main():
    """
    Runs all pure merge or impure merge querying, writes into a working directory all the timings and results for
    individual substructure pairs, and then processes it so it is in a format for filtering per fragment pair
    (the same substructure pairs may exist for multiple fragment pairs so we only want to run these queries once to
    remove redundancy - but they will be filtered separately)
    """
    parser = ArgumentParser()
    parser.add_argument('--substructure_pair_file', default=None,
                        help="substructure_pairs.json file generated by runEnumeration.py, if not provided, checks within substructure_dir")
    parser.add_argument('--descriptor', required=False, default=config.DESCRIPTOR_NAME, choices=['prop_pharmfp'],
                        help="descriptor property for bioisosteric merging search, currently only supports prop_pharmfp!")
    parser.add_argument('--pure_search', action='store_true', help="flag for if you want to search for only pure merges")
    parser.add_argument('--working_dir', required=False, default=config.WORKING_DIR, help="where to save intermediate files")
    parser.add_argument('--output_dir', required=False, default=config.OUTPUT_DIR, help="where to save output files")
    parser.add_argument('--n_parallel', required=False, type=int, help="number of parallel queries to run")
    parser.add_argument('--run_parallel', action='store_true', help="whether to run queries in parallel")
    parser.add_argument('--clear_working_dir', action='store_true', help="whether to clear the contents of the working dir")
    parser.add_argument('--limit', default=None, help="limit to the number of compounds to retrieve per query")
    parser.add_argument('--prolif_prioritization', action='store_true', help="whether to prioritize the substructures used for those that make an interactions")
    parser.add_argument('--target', required=False, default=None, help="name of the target (in fragalysis format)")
    parser.add_argument('--substructure_dir', required=False, default=None, help="location of the enumeration dir creater by runEnumeration.py")
    parser.add_argument('--max_prioritize', type=int, required=False, default=None, help="if prioritizing, the max number of compounds to prioritize")
    parser.add_argument('--r_group_search', required=False, default=None, help="whether you want to run an r group search on the compounds automatically using R groups from the original fragments")
    parser.add_argument('-fd', '--fragalysis_dir', required=False, default=base_config.FRAGALYSIS_DATA_DIR, help="location of fragalysis structures for prolif prioritization")
    args = parser.parse_args()

    if not os.path.exists(args.working_dir):
        os.mkdir(args.working_dir)

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    if not args.substructure_pair_file:
        substructure_file = os.path.join(args.substructure_dir, 'substructure_pairs.json')
    else:
        substructure_file = args.substructure_pair_file

    # if r group expansions are being run
    r_group_data_file = None  # created by runEnumeration.py (if flag is on)
    if args.r_group_search:
        r_group_data_file = os.path.join(args.substructure_dir, 'r_group_expansions.json')

    runKnitting(substructure_file, args.n_parallel, args.target, args.working_dir, args.output_dir, args.limit,
                args.substructure_dir, args.descriptor, args.pure_search, args.prolif_prioritization, args.max_prioritize,
                args.r_group_search, r_group_data_file, fragalysis_dir=args.fragalysis_dir)

    if args.clear_working_dir:
        shutil.rmtree(args.working_dir)


if __name__ == "__main__":
    main()
