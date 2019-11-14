#!/usr/bin/env ost
""" Check if the domains cover residues scored by lDDT-BS.

The aim is to find out if can we use some residues of the domain as "functional"
sites, and do we still have enough residues covered by lDDT-BS?
"""
from __future__ import print_function
import argparse
import json
import string

import pandas


def parse_args():
    """Parse options."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose',
                        dest='verbosity',
                        action='append_const',
                        const=1,
                        default=[ost.GetVerbosityLevel()],
                        help="increase verbosity. Can be used multiple times.")
    parser.add_argument('-q', '--quiet',
                        dest='verbosity',
                        action='append_const',
                        const=-1,
                        help="decrease verbosity. Can be used multiple times.")
    parser.add_argument('-o', '--out',
                        default='-',
                        help="Output file (CSV format).")
    parser.add_argument('-m', '--models-path',
                        default="cathsm_models",
                        help="Path to the folder containing the CATH-SM "
                             + "models.")
    parser.add_argument('-s', '--sm-models-path',
                        default="sm_models",
                        help="Path to the folder containing the SWISS-MODEL "
                             + "models.")
    parser.add_argument('-t', '--targets-path',
                        default="cameo_targets",
                        help=("Path to the folder containing the CAMEO target"
                              + " data."))
    parser.add_argument('-c', '--crh',
                        default="cameo_models.allhits.crh",
                        help=("The ranked hits file."))
    parser.add_argument('targets',
                        nargs="+",
                        help=("Name(s) of the target(s) to process, such as "
                              + "20180825_00000040_2."))

    args = parser.parse_args()

    # Set verbosity
    args.verbosity = sum(args.verbosity)
    ost.PushVerbosityLevel(args.verbosity)

    return args


def read_crh(filename):
    """ Read the Ranked Hits from CATH Resolve Hits file. """
    crh = pandas.read_table(filename, sep=" ", comment="#", header=None)   

    column_names = ['target', 'hit', 'normalised_qmean4','query_range']
    if len(crh.columns) == 5:
        crh.columns = column_names + ['query_range_duplicated']
    else:
        crh.columns = column_names

    return crh


def get_query_range(crh, target, hit):
    """ Get the query range that was used in this hit.  Return a tuple (start,
    end). """
    range_str = crh[(crh['target'] == target) & (crh['hit'] == hit)]['query_range']
    assert len(range_str) == 1  # a Series with 1 row
    range_list = range_str.iloc[0].split("-")
    return tuple(range_list)


def iter_hits(crh, target):
    """ Return a generator that yields the hits ranked according to the
    normalized qmean4 value."""
    hits = crh[crh['target'] == target]
    hits = hits.sort_values(by=['normalised_qmean4'], ascending=False)
    for hit in hits['hit']:
        yield hit


def split_comma(csv):
    """ Split a string of comma-delimited values with optional spaces into a list"""
    if csv is None or csv == "":
        return []
    else:
        return [item.strip() for item in csv.split(",")]


def assess_target(target, args, crh):
    """ Assess coverage of the binding site"""
    ost.LogDebug(target)
    if target == "20180811_00000008_1":
        ost.LogWarning("Skipping 20180811_00000008_1 (5Y2D) missing from SWISS-MODEL data.")
        return []

    results = []

    json_path = os.path.join(args.targets_path, target, "lddtha_bs.json")
    with open(json_path) as jfile:
        jdata = json.load(jfile)

    # Get all the hits
    for rank, hit in enumerate(iter_hits(crh, target), start=1):
        ost.LogDebug(hit)
        query_range = get_query_range(crh, target, hit)

        if jdata['results']['score_lddtha_bs'] is None:
            ost.LogScript("No ligand for %s, %s" % (target, hit))
            continue

        best_coverage = 0

        # Get all the ligands
        for ligand_name, ligand_data in jdata['results']['score_lddtha_bs'].items():
            ost.LogDebug(ligand_name)
            trg_res = split_comma(ligand_data['trg_bs_res'])
            if len(trg_res) == 0:
                ost.LogWarning("Empty target binding site in %s %s, ignoring." %(target, ligand_name))
                continue
            mdl_res = split_comma(ligand_data['mdl_bs_res'])
            # How much of BS is covered by domain?
            reschains = [res.split('.')[0] for res in mdl_res]
            # Only allow 3 letter codes with 3 uppercase letters
            # We only expect the 20 standard residues
            # This is to make sure we don't accidentally have residues named like 8OG
            trg_resnames = [res.split('.')[1].strip(string.digits) for res in trg_res]
            mdl_resnames = [res.split('.')[1].strip(string.digits) for res in mdl_res]
            assert all([len(resname) == 3 for resname in trg_resnames])
            assert all([len(resname) == 3 for resname in mdl_resnames])

            # Make sure we only have 1 chain
            assert all(reschain == reschains[0] for reschain in reschains)

            # Now look at the target
            trg_resnums = [int(res.split('.')[1].strip(string.ascii_uppercase)) for res in trg_res]
            domain_start, domain_stop = [int(x) for x in query_range]
            bs_covered = 0
            bs_size = 0
            for trg_resnum in trg_resnums:
                bs_size += 1
                if trg_resnum >= domain_start and trg_resnum <= domain_stop:
                    bs_covered += 1

            new_coverage = float(bs_covered) / float(bs_size)
            if new_coverage > best_coverage:
                best_coverage = new_coverage

            results.append({
                'target': target,
                'hit': hit,
                'hit_rank': rank,
                'hit_from': query_range[0],
                'hit_to': query_range[1],
                'ligand_name': ligand_name,
                'bs_size': bs_size,
                'bs_covered': bs_covered,
                'coverage': new_coverage,
            })

        if best_coverage == 0:
            ost.LogScript("No coverage for %s, %s" % (target, hit))
        else:
            ost.LogScript("Coverage for %s, %s: %s" %(target, hit, best_coverage))

    return results


def _main():
    args = parse_args()
    crh = read_crh(args.crh)

    results = []
    for target in args.targets:
        results.extend(assess_target(target, args, crh))

    table = pandas.DataFrame.from_dict(results)
    if args.out == '-':
        args.out = sys.stdout
    table.to_csv(args.out)


if __name__ == "__main__":
    _main()
