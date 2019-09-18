# CATH-SM / CAMEO benchmarking

This repo contains benchmarking results of the CATH-SM modeling pipeline
on CAMEO targets from August 2018 to January 2019.

The CAMEO targets were submitted to a version of the that was slightly
altered to allow more than one hit per region, so this has been resolved
to the single "best" hit using `QMEAN4` and `cath-resolve-hits`.

The SWISS-MODEL models are the corresponding models that have been submitted
to CAMEO.

The evaluation was performed in a way similar to CAMEO with the
[OpenStructure `compare-structure` action](https://www.openstructure.org/docs/1.10/actions/#comparing-two-structures).
The assessment is similar to what CAMEO does, but  because we are specifically
interested in domain modeling we cut the target structures to the domain
of interest (as determined by FunFams) in order to score only the domain part
of the target structure.
The normal CAMEO evaluation, which weights lDDT results by coverage,
would negatively penalize long targets and/or targets with short
domains.


## CATH-SM model overview

There are 752 predicted 3D models from 203 CAMEO target sequences. Each
model has two files:

* JSON file containing the data sent to SWISS-MODEL
* PDB file sent back from SWISS-MODEL

eg

```sh
cameo_models/20180811_00000087_1/20180811_00000087_1.hit1.json
cameo_models/20180811_00000087_1/20180811_00000087_1.hit1.pdb
```

A post-processing script was used to collate these data files and
generate a summary of all hits for all targets.

```sh
python3 ./scripts/process_models.py \
    --in cameo_models/ > cameo_models.allhits.crh
```

The format of this file is:

```sh
target  hit  normalised_qmean4   query_range
20180811_00000087_1 hit1 94.91 10-270
20180811_00000087_1 hit2 97.77 87-139
```

Note: "normalised" QMEAN4 score = QMEAN4 + 100 as the next stage needs
this number to be strictly positive.

We can now use `cath-resolve-hits` to select the best hit per region
and provide a pretty visualisation of the decision making process.

```sh
cath-resolve-hits \
    --html-output-to-file cameo_models.html \
    cameo_models.allhits.crh > cameo_models.besthits.crh
```

## Benchmarking

The benchmarking data (lDDT values) in `benchmark_lddt.csv` was generated by
the following [OpenStructure](https://openstructure.org/) script:

```sh
all_targets=$(ls cathsm_models)
ost ./scripts/benchmark_cathsm.py -o benchmark_lddt.csv $all_targets
```

Further analysis of the data and plots were performed in the
`Plots_lDDT.ipynb` Jupyter notebook.
