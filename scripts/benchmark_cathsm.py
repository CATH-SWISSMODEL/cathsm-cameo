#!/usr/bin/env ost
""" This script generates a CSV file with lDDT values. """
from __future__ import print_function
import argparse
from contextlib import contextmanager
from glob import glob
import json
import os
import pandas
import subprocess
import tempfile
import sys

import ost


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


@contextmanager
def temporary_file(**kwargs):
    """Context manager to create temporary file and close and delete it after
    use."""
    tf = tempfile.NamedTemporaryFile(**kwargs)
    try:
        ost.LogVerbose("Making temporary file %s" % tf.name)
        yield tf
    finally:
        ost.LogVerbose("Closing temporary file %s" % tf.name)
        tf.close()
        if not kwargs.get("delete", True):
            ost.LogVerbose("Deleting temporary file %s" % tf.name)
            try:
                os.remove(tf.name)
            except Exception as e:
                ost.LogVerbose("Couldn;t remove file %s: %s" % (tf.name,
                                                                str(e)))


def process(json_file, args, cc):
    """
    :param cc: the chain cluster index :class:`~sm.smtl.ChainClusters`
    """
    pass


def read_crh(filename):
    """ Read the Ranked Hits from CATH Resolve Hits file. """
    crh = pandas.read_table(filename, sep=" ",
                            names=['target', 'hit',
                                   'normalised_qmean4',
                                   'query_range'])
    return crh


def iter_hits(crh, target):
    """ Return a generator that yields the hits ranked according to the
    normalized qmean4 value."""
    hits = crh[crh['target'] == target]
    hits = hits.sort_values(by=['normalised_qmean4'], ascending=False)
    for hit in hits['hit']:
        yield hit


def iter_sm_models(base, target):
    """ Return a generator that yields the SWISS-MODEL models for the target
    in the order they were returned to CAMEO.
    """
    models_glob = os.path.join(base, target, "part-*.pdb")
    model_filenames = sorted(glob(models_glob))
    for model_filename in model_filenames:
        yield ost.io.LoadPDB(model_filename)


def get_query_range(crh, target, hit):
    """ Get the query range that was used in this hit.  Return a tuple (start,
    end). """
    range_str = crh[(crh['target'] == target) & (crh['hit'] == hit)]['query_range']
    assert len(range_str) == 1  # a Series with 1 row
    range_list = range_str.iloc[0].split("-")
    return tuple(range_list)


def get_cathsm_model(base, target, hit):
    """ Get the CATH-SM model as an ost Entity. The model is re-numbered
    according to the query range in order to be directly comparable to the
    target, as required by CAMEO."""
    model_path = os.path.join(base, target, "%s.%s.pdb" % (target, hit))
    json_path = os.path.join(base, target, "%s.%s.json" % (target, hit))
    metadata = json.load(open(json_path))
    model = ost.io.LoadPDB(model_path)
    # Model is numbered based on query_sequence. We need to renumber it
    # so that it corresponds to the CAMEO sequence. In this case we need the
    # first residue of the model to have the number of the start of query_range.
    new_start = int(metadata['query_range'].split('-')[0])
    renumber_model(model, new_start)
    return model


def get_sm_model(base, target):
    """ Get the structure of the top-ranked SWISS-MODEL model. """
    model_glob = os.path.join(base, target, "part-*.pdb")
    models_path = sorted(glob(model_glob))
    top_model = models_path[0]
    return ost.io.LoadPDB(top_model)


def renumber_model(model, new_start):
    """ Renumber the model in place, using the new start number. """
    peptide_chains = [
        chain for chain in model.chains if chain.name not in set(["-", "_"])]
    assert len(peptide_chains) == 1  #
    target_chain = peptide_chains[0]
    assert target_chain.IsValid()
    editor = model.EditICS()
    ost.LogInfo("Renumbering target from %s" % new_start)
    editor.RenumberChain(target_chain, new_start, True)
    editor.UpdateXCS()


def get_target_structures(base, target):
    """ Get the structures (BioUnits) of this target. Returns a list of
    targets."""
    target_glob = os.path.join(base, target, "bu_target_*.pdb")
    targets = []
    for biounit in glob(target_glob):
        targets.append(ost.io.LoadPDB(biounit))
    return targets


def get_target_structure(base, target):
    """ Get the target structure of this target."""
    target_filename = os.path.join(base, target, "target.pdb")
    return ost.io.LoadPDB(target_filename)


def get_target_alignment(base, target):
    """ Get the length of the target sequence."""
    target_aln_path = os.path.join(base, target, "target.aln")
    seqlist = ost.io.LoadSequenceList(target_aln_path, format='fasta')
    assert len(seqlist) == 2
    aln = ost.seq.CreateAlignment()
    for seq in seqlist:
        aln.AddSequence(seq)
    return aln


def cut_targets(targets, query_range):
    """ Cut the targets to the query range. Return a list of targets restricted
    to the query range."""
    target_views = []
    for target in targets:
        target_views.append(cut_target(target))
    return target_views


def cut_target(target, query_range):
    """ Cut the target to the query range. Return a view of the target with the
    range selected."""
    sele_str = 'rnum>=%s and rnum<=%s' % query_range
    sele_str += ' and cname!="_" and cname!="-"'  # no ligand or water
    return target.Select(sele_str)


def compare_structures(model, reference):
    """ Compare the structures and return an lDDT value. """
    with temporary_file(delete=False) as model_file, \
            temporary_file(delete=False) as target_file, \
            temporary_file(delete=False) as output_file:
        # We'll open/close the file later by name so we don't want python
        # to mess up
        model_file.close()
        target_file.close()
        output_file.close()
        ost.io.SavePDB(model, model_file.name)
        ost.io.SavePDB(reference, target_file.name)
        cmd = ['/usr/bin/env', 'ost', 'compare-structures',
               '-r', target_file.name,
               '-m', model_file.name,
               '-o', output_file.name,
               '--lddt',
               '-v', str(ost.GetVerbosityLevel()),
               '--residue-number-alignment'  # do not re-align
               ]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            ost.LogError("ost compare-structures returned %s, error stream included below " % proc.returncode)
            ost.LogError(stderr)
            if reference.atom_count == 0:
                # Reference is empty, probably because it covered an other
                # portion of the target sequence that doesn't include our
                # domain of interest. We keep track of this hit with missing
                # values (None)
                ost.LogError("The reference structure appears to be empty. "
                             + "This probably means it didn't cover this "
                             + "domain. Returning missing values.")
                return {'weighted_lddt': None,
                        'oligo_lddt': None,
                        'single_chain_lddt': None,
                        'coverage': None,
                        }
            else:
                raise RuntimeError("ost compare-structure returned %s" % proc.returncode)
        compare_data = json.load(open(output_file.name))
        try:
            compare_results = compare_data['result'][
                os.path.basename(model_file.name)][
                os.path.basename(target_file.name)]
        except KeyError:
            ost.LogError("ost compare-structures retured a success code but "
                         + "the expected data is not available. This is "
                         + "probably caused by an acceptable error such as "
                         + "a failure in chain mapping. Returning missing "
                         + "values. STDERR output included below: ")
            ost.LogError(stderr)
            return {'weighted_lddt': None,
                    'oligo_lddt': None,
                    'single_chain_lddt': None,
                    'coverage': None,
                    }

        # What is the coverage of the model?
        coverages = []
        for aln_str in compare_results['info']['mapping']['alignments']:
            aln = ost.io.AlignmentFromString(str(aln_str), 'fasta')
            coverages.append(aln.GetCoverage(1))
            assert aln.GetCount() == 2

        lddt_results = compare_results['lddt']
        single_chain_lddts = []
        for sc_lddt in  lddt_results['single_chain_lddt']:
            assert sc_lddt['status'] == 'SUCCESS'
            single_chain_lddts.append(sc_lddt['global_score'])

        assert lddt_results['oligo_lddt']['status'] == 'SUCCESS'
        assert lddt_results['weighted_lddt']['status'] == 'SUCCESS'

        return {'weighted_lddt': lddt_results['weighted_lddt']['global_score'],
                'oligo_lddt': lddt_results['oligo_lddt']['global_score'],
                'single_chain_lddt': max(single_chain_lddts),
                'coverage': max(coverages),
                }


def assess_target(target, args, crh):
    """ Assess SM and CATH-SM models against a target. """
    results = []
    ost.LogScript(target)
    target_structure = get_target_structure(args.targets_path, target)
    for rank, hit in enumerate(iter_hits(crh, target), start=1):
        ost.LogScript(hit)
        query_range = get_query_range(crh, target, hit)
        cathsm_model = get_cathsm_model(args.models_path, target, hit)
        target_for_hit = cut_target(target_structure, query_range)

        # Compare
        lddt_cathsm = compare_structures(cathsm_model, target_for_hit)

        for sm_rank, sm_model in enumerate(
                iter_sm_models(args.sm_models_path, target),
                start=1):
            ost.LogScript("SM model %s" % sm_rank)
            lddt_sm = compare_structures(sm_model, target_for_hit)
            results.append({
                'target': target,
                'target_length':
                    len(get_target_alignment(args.targets_path, target)),
                'hit': hit,
                'hit_rank': rank,
                'hit_from': query_range[0],
                'hit_to': query_range[1],
                'sm_oligo_lddt': lddt_sm['oligo_lddt'],
                'cathsm_oligo_lddt': lddt_cathsm['oligo_lddt'],
                'sm_weighted_lddt': lddt_sm['weighted_lddt'],
                'cathsm_weighted_lddt': lddt_cathsm['weighted_lddt'],
                'sm_single_chain_lddt': lddt_sm['single_chain_lddt'],
                'cathsm_single_chain_lddt': lddt_cathsm['single_chain_lddt'],
                'sm_coverage': lddt_sm['coverage'],
                'cathsm_coverage': lddt_cathsm['coverage'],
                'sm_rank': sm_rank,
            })
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
