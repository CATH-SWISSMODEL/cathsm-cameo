#!/usr/bin/env python3

"""
CLI tool to run some post-processing on 3D models

Currently this collates the hits for each CAMEO model and
generates a file containing the query range and a normalised
QMEAN4 score. This file can be used as the input to 
cath-resolve-hits (CRH).

Note: cath-resolve-hits (CRH) expects a strictly positive score,
so the score reported in the output is QMEAN4 + 100.
"""

# core
import argparse
import logging
import json
import os
import re
import string


# logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S')
LOG = logging.getLogger(__name__)

parser = argparse.ArgumentParser()
parser.add_argument('--in', '-i', dest='in_dir', type=str, required=True,
                    help='directory containing the models')

class SubmitModelData:
    def __init__(self, *, auth_asym_id, pdb_id, query_range, target_sequence, template_seqres_offset, template_sequence):
        self.auth_asym_id = auth_asym_id
        self.pdb_id = pdb_id
        self.query_range = query_range
        self.target_sequence = target_sequence
        self.template_seqres_offset = template_seqres_offset
        self.template_sequence = template_sequence


def qmean4_from_pdb(pdb_file):
    qmean4 = None
    with open(pdb_file) as pdb_io:
        for line in pdb_io:
            match = re.match(r'REMARK\s+\d+\s+QMN4\s+([\-\d.]+)', line)
            if match:
                qmean4 = match.group(1)
                break
    return qmean4


if __name__ == '__main__':
    args = parser.parse_args()
    
    models = {}
    for (dirpath, dirnames, filenames) in os.walk(args.in_dir):

        pdb_files = {f[:-4]: os.path.join(dirpath, f) for f in filenames if f.endswith('.pdb')}

        for model_and_hit_id, pdb_file in pdb_files.items():
            json_file = re.sub('\.pdb$', '.json', pdb_file)

            model_id, hit_id = model_and_hit_id.split('.')

            with open(json_file) as json_io:
                model_data = json.load(json_io)
                model = SubmitModelData(**model_data)

            qmean4 = qmean4_from_pdb(pdb_file)

            if model_id not in models:
                models[model_id] = []

            models[model_id].extend([{
                'hit': hit_id,
                'qmean4': qmean4,
                **model.__dict__,
            }])
    
    for model_id in models:
        hits = models[model_id]
        for hit in hits:
            print("{} {} {} {}".format(
                model_id, 
                hit['hit'], 
                str(100 + float(hit['qmean4'])), 
                hit['query_range']))
