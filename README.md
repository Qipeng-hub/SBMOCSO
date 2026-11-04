# SBMOCSO

Sparse binary multi-objective competitive swarm optimizer for bi-objective gene selection on ophthalmic
transcriptomic data, and the code of the experiments in the accompanying paper.

## Requirements

Python 3.10 or later with the packages in `requirements.txt`:

```
pip install -r requirements.txt
```

## Layout

```
sbmocso/                  algorithms and problem
  core.py                 population container and budgeted evaluation
  operators.py            non-dominated sorting, crowding distance, selection and variation operators
  problem.py              bi-objective gene selection problem (fraction of genes, balanced 3-NN error)
  metrics.py              non-dominated filtering and two-dimensional hypervolume
  algorithms/             SBMOCSO and the eight compared algorithms
scripts/
  prepare_data.py         builds data/processed from the GEO files
  run_one.py              one algorithm, one dataset, one split, one budget
  run_all.py              all runs of the paper
  analysis.py             statistical tests and tables
  figures.py              result figures
  gene_check.py           database check of the selected genes
  gene_literature.py      literature search of the most frequently selected genes
data/
  processed/              the six datasets (2000 most variable genes each)
  gene_check/             Open Targets scores (data release 26.09) and OpenAlex search results
```

## Data

`data/processed` contains the six datasets used in the paper. They can be rebuilt from the public GEO
series GSE29801, GSE22138, GSE44295, GSE58331 and GSE160306: place the files listed in
`scripts/prepare_data.py` in `data/raw` and run

```
python scripts/prepare_data.py
```

## Usage

One run of SBMOCSO on the macular retina AMD dataset, split 0, 10,000 evaluations:

```
python scripts/run_one.py SBMOCSO AMD-RET-GSE29801 0 10000
```

Variants are selected with `key=value` options, for example
`init=uniform`, `full_copy=true`, `mutate=false` or `credit=true`.
Each run writes `results/<tag>/<algorithm>/<dataset>/r<run>.npz` with the final population, its
training objectives, test errors, hypervolumes and the hypervolume curve.

Using SBMOCSO on another problem:

```python
from sbmocso.algorithms import sbmocso
from sbmocso.core import MOTask
from sbmocso.problem import GeneSelection

problem = GeneSelection("UM-GSE22138", run=0)
task = MOTask(problem, N=100, max_fes=10000, seed=1)
result = sbmocso.run(task)          # result.decs, result.objs: final population
```

Any problem object with the attributes `M`, `D`, `lower`, `upper`, `encoding` and the methods
`cal_obj(dec)` and `cal_con(dec)` can be passed to `MOTask`.

## Reproducing the paper

```
python scripts/run_all.py --workers 32     # 8820 runs: main comparison, budgets, ablation, credit variant
python scripts/analysis.py                  # output/tables, output/summary.json
python scripts/figures.py                   # output/figures (Figures 3-9)
python scripts/gene_check.py counts --gene-info Homo_sapiens.gene_info.gz
python scripts/gene_check.py test           # Table S3, part A
python scripts/gene_check.py figure         # Figure 10
```

`python scripts/run_all.py --list` prints the command line of every run, which can be distributed over
a cluster. Every algorithm uses the population size 100 and the default parameters listed in the paper.
The train/test split of run `r` and the random seed of each algorithm are fixed in
`sbmocso/problem.py` and `scripts/run_one.py`, so the runs are reproducible.
