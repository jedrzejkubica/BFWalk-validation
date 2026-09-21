# BFWalk-validation

This repository contains scripts for the validation of **[BFWalk](https://github.com/jedrzejkubica/BFWalk)** as described in the submitted manuscript [DOI TBA]. We performed leave-one-out cross-validation (LOO CV) and tissue enrichment validation to compare the performance of BFWalk with RWR (as implemented in MultiXrank[^1]) and NetCore[^2].


### Step 1. Run LOO CV for BFWalk

We assume that BFWalk is installed as described in [BFWalk](https://github.com/jedrzejkubica/BFWalk) and that input data (interactome, seeds) is prepared as described in [BFWalk Interactome](https://github.com/jedrzejkubica/BFWalk/tree/main/Interactome). We provide the interactome and seeds for three studied phenotypes (dyschromatopsia, hypertrophic cardiomyopathy and chronic kidney disease) in [data/](data/). Scores (`scores_LOO.tsv`) and ranks (`ranks_LOO.tsv`) for the left-out genes will be saved in `--out`, logs will be in `log_LOO.txt`.

```
( for PHENO in DYSCHROM HYPCARD CKD ; do
    mkdir -p BFWalk-output/${PHENO}/ ;
done )
```

```
( for PHENO in DYSCHROM HYPCARD CKD ; do
    python run_BFWalk/leave_one_out.py \
        --network data/interactome_human_evidence2_OR_direct1.sif \
        --seeds data/causal_proteins_${PHENO}.txt \
        --out BFWalk-output/${PHENO}/ \
        2>BFWalk-output/${PHENO}/log_LOO.txt ;
done )
```


### Step 2. Run LOO CV for MultiXrank

Create a Python environment and install MultiXrank as described in [https://github.com/anthbapt/multixrank](https://github.com/anthbapt/multixrank).

A default config file is provided at [run_multixrank/default/config.yml](run_multixrank/default/config.yml). All output files will be saved in `multixrank-output/`.

```
( for PHENO in DYSCHROM HYPCARD CKD ; do
    mkdir -p multixrank-output/${PHENO}/ ;
done )
```

The MultiXrank scoring script takes the interactome SIF file and BFWalk-derived left-out genes' ranks as input. The interactome SIF file from BFWalk will be automatically converted to a TSV file (with two columns: node1, node2) as required by MultiXrank. The scoring script produces scores for all genes in the interactome, they will be saved in `multixrank-output/multiplex_1.tsv`.

```
( for PHENO in DYSCHROM HYPCARD CKD ; do
    python run_multixrank/run_multixrank.py \
        --network data/interactome_human_evidence2_OR_direct1.sif \
        --BFWalk_ranks BFWalk-output/${PHENO}/ranks_LOO.tsv \
        --config run_multixrank/default/config.yml \
        --out multixrank-output/${PHENO}/ ;
done )
```

Then, run leave-one-out with the same input files as in the previous step. One subdirectory for each left-out gene will created for temporary files.

```
( for PHENO in DYSCHROM HYPCARD CKD ; do
    python run_multixrank/run_leave_one_out.py \
        --network data/interactome_human_evidence2_OR_direct1.sif \
        --BFWalk_ranks BFWalk-output/${PHENO}/ranks_LOO.tsv \
        --config run_multixrank/default/config.yml \
        --out multixrank-output/${PHENO}/ ;
done )
```


### Step 3. Run LOO CV for NetCore

Create a Python environment and install NetCore in `run_netcore/` following instructions at [https://github.molgen.mpg.de/barel/NetCore](https://github.molgen.mpg.de/barel/NetCore) (note that NetCore requires Python 3.7).

The interactome TSV file created for MultiXrank will be reused here. All output files will be saved in `netcore-output/`.

```
mkdir netcore-output/
```

First, NetCore requires running edge permutations on the interactome. This script takes the interactome TSV file and produces a subdirectory `permutations/` in `netcore-output/`, it will be used as input in the scoring step, log will be in `log_permut.txt`.

```
python run_netcore/run_permutations.py \
    --interactome multixrank-output/${PHENO}/interactome_human_evidence2_OR_direct1.tsv \
    --output-path netcore-output/permutations/ \
    2>netcore-output/log_permut.txt
```

Then, run the NetCore scoring script with an interactome TSV file (`-e`),the seeds file (`-c`) and the permutations directory (`-pd`). It saves scores for all genes in the network in `netcore-output/random_walk_weights.txt` (`-o`). Both stdout and stderr will be captured in `log.txt`.

```
( for PHENO in DYSCHROM HYPCARD CKD ; do
    mkdir -p netcore-output/${PHENO}/output/ ;
    python run_netcore/NetCore/netcore/netcore.py \
            -e multixrank-output/${PHENO}/interactome_human_evidence2_OR_direct1.tsv \
            -s data/causal_proteins_${PHENO}.txt \
            -pd netcore-output/permutations/interactome_human_edge_permutations/ \
            -o netcore-output/${PHENO}/ \
            1>netcore-output/${PHENO}/log.txt \
            2>&1 ;
done )
```

Then, run the leave-one-out bash script, it follows this usage:

```
./run_leave_one_out.sh <interactome> <seeds_dir> <permutation_dir> <out_dir>
```

It takes arguments in the following order: interactome TSV file (produced for MultiXrank), the temporary files produced by MultiXrank (in `multixrank-output/tmp/`), the permutations subdirectory and the output directory. One subdirectory will be created for each left-out gene in the output directory, each containing `random_walk_weights.txt` with genes scores.

```
( for PHENO in DYSCHROM HYPCARD CKD ; do
    run_netcore/run_leave_one_out.sh \
        multixrank-output/${PHENO}/interactome_human_evidence2_OR_direct1.tsv \
        multixrank-output/${PHENO}/tmp/ \
        netcore-output/permutations/interactome_human_edge_permutations/ \
        netcore-output/${PHENO}/ ;
done )
```


### Part 2. Perform the analyses

This part covers three analyses to compare BFWalk, MultiXrank and NetCore.

```
mkdir validation_logs ;
mkdir figures ;
```

[validation_CDF.py](validation_CDF.py) calculates and plots a cumulative distribution function (CDF) for left-out gene ranks. CDF curves show the proportion of left-out genes recovered at or above rank x, for every rank x. The area under the curve (AUC) is calculated and shown in the legend.

```
python validation_CDF.py --help
```

```
( for PHENO in DYSCHROM HYPCARD CKD ; do
    python validation_CDF.py \
       --network data/interactome_human_evidence2_OR_direct1.sif \
       --BFWalk_ranks BFWalk-output/${PHENO}/ranks_LOO.tsv \
       --multixrank_ranks multixrank-output/${PHENO}/ranks_LOO.tsv \
       --netcore_LOO_dir netcore-output/${PHENO}/ \
       --cdf figures/${PHENO}_CDF.png \
       --pheno ${PHENO} \
       2>>validation_logs/log_CDF ;
done )

```

[validation_TE.py](validation_TE.py) compares the ratio of the highest-scoring genes enriched in the tissue of interest with the ratio of all genes. It checks whether the highest-scoring genes are more enriched in the tissue than would be expected by chance.

> [!NOTE]
> For the tissue enrichment validation we downloaded the GTEx data (release V8) from the EBI Expression Atlas.

```
python validation_TE.py --help
```

```
( for PHENO in DYSCHROM HYPCARD CKD ; do
    if [ "$PHENO" == "CKD" ]; then
       TISSUE="cortex of kidney"
    else
       TISSUE="heart left ventricle"
    fi

    python validation_TE.py \
        --network data/interactome_human_evidence2_OR_direct1.sif \
        --uniprot data/uniprot_parsed.tsv \
        --gtex data/E-GTEX-8-query-results.tpmss.tsv \
        --tissue="${TISSUE}" \
        --BFWalk_scores BFWalk-output/${PHENO}/scores.tsv \
        --multixrank_scores multixrank-output/${PHENO}/multiplex_1.tsv \
        --netcore_scores netcore-output/${PHENO}/output/random_walk_weights.txt \
        --matrix figures/${PHENO}_comparison_matrix.png \
        2>>validation_logs/log_TE ;
done )
```

[validation_ranksVsDeg.py](validation_ranksVsDeg.py) examines the relationship between the degrees of the left-out genes and the differences in their ranks between BFWalk and MultiXrank, and between BFWalk and NetCore.

```
python validation_ranksVsDeg.py --help
```

```
python validation_rankVsDeg.py \
    --network data/interactome_human_evidence2_OR_direct1.sif \
    --phenotypes HYPCARD DYSCHROM CKD \
    --BFWalk_out_dir BFWalk-output/ \
    --multixrank_out_dir multixrank-output/ \
    --netcore_out_dir netcore-output/ \
    --rankVsDeg figures/ \
    2>>validation_logs/log_rankVsDeg
```


### Dependencies

For validation we used Python 3.9 with the following libraries:
- numpy 1.23
- networkx 3.2
- matplotlib 3.4
- scipy 1.13


### References

[^1]: Baptista, A., Gonzalez, A., & Baudot, A. (2022). Universal multilayer network exploration by random walk with restart. Communications Physics, 5(1), 1–9. https://doi.org/10.1038/s42005-022-00937-9
[^2]: Gal Barel, Ralf Herwig, NetCore: a network propagation approach using node coreness, Nucleic Acids Research, Volume 48, Issue 17, 25 September 2020, Page e98, https://doi.org/10.1093/nar/gkaa639