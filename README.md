# CATH-SM / CAMEO models

This repo contains the results of running a set of CAMEO targets
through the CATH-SM modelling pipeline. The pipeline was altered
to allow more than one hit per region, so this has been resolved
to the single "best" hit using `QMEAN4` and `cath-resolve-hits`.

## Overview

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
