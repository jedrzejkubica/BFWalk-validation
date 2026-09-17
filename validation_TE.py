###########################################################################################
# Copyright (C) Jędrzej Kubica, Nicolas Thierry-Mieg, 2026
#
# This file was written by Jędrzej Kubica and Nicolas Thierry-Mieg
# (CNRS, France) Nicolas.Thierry-Mieg@univ-grenoble-alpes.fr
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <https://www.gnu.org/licenses/>.
############################################################################################

import os
import sys
import argparse
import pathlib
import logging

import networkx
import scipy
import matplotlib.pyplot

sys.path.append("/home/kubicaj/Software/BFWalk")
import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)

matplotlib.pyplot.set_loglevel('warning')


def parse_uniprot(uniprot_file):
    """
    Parse a TSV file from uniprot_parser.py with columns:
    - Uniprot Primary AC
    - Uniprot Secondary AC(s) (comma-separated)
    - tax ID
    - gene name
    - synonyms (comma-separated)

    Returns:
      - uniprot2_gene: dict with key=Uniprot Primary AC, value=gene name
    """
    uniprot2_gene = {}

    with open(uniprot_file, 'r') as f:
        # skip header
        header = f.readline()
        if not header.startswith("PrimaryAC\t"):
            raise Exception(f"uniprot file {uniprot_file} is headerless? expecting headers but got {header}")

        for line in f:
            line_split = line.rstrip("\n").split("\t")

            # if some records are incomplete, die
            if(len(line_split) != 5):
                raise Exception(f"Bad line in the uniprot file {uniprot_file}, not 5 tab-separated fields")

            (primaryAC, secondaryACs, taxID, gene, synonyms) = line_split
            # only keep human proteins
            if taxID != "9606":
                continue

            # make sure there is a gene name
            if gene == "":
                continue
            uniprot2_gene[primaryAC] = gene

    return(uniprot2_gene)


def parse_expression_data(expression_data, tissue):
    """
    Parses a TSV file created by the Ensembl Expression Atlas:
    Gene Name is the second column, followed by the TPMs in each tissue.
    Provided tissue must match one of the tissues listed in the header of expression_data.

    arguments:
    - expression_data: path to the TSV file from Expression Atlas
    - tissue: string, tissue of interest.

    returns:
    - gene2enrichment: dict, key=gene name, value=tissue enrichment, ie TPM in tissue divided
      by average TPM in all tissues.
    """
    gene2enrichment = {}

    try:
        f = open(expression_data, 'r')
    except Exception as e:
        raise Exception("cannot open provided expression file")

    line = f.readline()
    if not line.startswith("# Expression Atlas"):
        raise Exception("Expression file problem")
        
    # skip header comments
    while line.startswith("#"):
        line = f.readline()
    # parse header
    split_line = line.rstrip().split('\t')
    if split_line[1] != 'Gene Name':
        raise Exception("Expression file header: not gene name column")
    tissueCol = -1
    for i in range(2, len(split_line)):
        if split_line[i] == tissue:
            tissueCol = i
            break
    if tissueCol == -1:
        raise Exception("Cannot find tissue in Expression file header")
    nbTissues = len(split_line) - 2

    # parse data lines
    for line in f:
        split_line = line.rstrip("\n").split("\t")
        gene = split_line[1]
        expTiss = split_line[tissueCol]
        # use 0 for missing data
        if expTiss == '':
            expTiss = 0.0
        else:
            expTiss = float(expTiss)
        expAverage = 0.0
        for i in range(2, len(split_line)):
            if split_line[i] != '':
                expAverage += float(split_line[i])
        expAverage /= nbTissues
        tissue_ratio = expTiss / expAverage
        gene2enrichment[gene] = tissue_ratio
    
    return(gene2enrichment)


def parse_BFWalk_scores(scores_file):
    """
    Parses a TSV ranks file (with header) with two columns: NODE, SCORE

    Returns:
    - node2score: dict with key=node, value=score
    """
    node2score = {}
    with open(scores_file, 'r') as f:
        header = f.readline()
        if not header.startswith("NODE\t"):
            raise Exception("Scores file problem, wrong format")

        for line in f:
            split_line = line.rstrip().split("\t")
            if len(split_line) != 2:
                raise Exception(f"Invalid line in scores file: {line}")

            (node, score) = split_line
            try:
                score = float(score)
            except ValueError:
                raise Exception(f"Invalid score value in scores file: {score} is not a float")
            node2score[node] = score

    return(node2score)


def parse_multixrank_scores(scores_file):
    """
    Parses a TSV ranks file (with header) with three columns: multiplex, node, score

    Returns:
    - node2score: dict with key=node, value=score
    """
    node2score = {}
    with open(scores_file, 'r') as f:
        header = f.readline()
        if not header.startswith("multiplex\t"):
            raise Exception("Scores file problem, wrong format")

        for line in f:
            split_line = line.rstrip().split("\t")
            if len(split_line) != 3:
                raise Exception(f"Invalid line in scores file: {line}")

            (multiplex, node, score) = split_line
            try:
                score = float(score)
            except ValueError:
                raise Exception(f"Invalid score value in scores file: {score} is not a float")
            node2score[node] = score

    return(node2score)


def parse_netcore_scores(scores_file):
    """
    Parses a TSV ranks file (with header) with four columns: node_index, node, prop_weight, pvalue

    Returns:
    - node2score: dict with key=node, value=score
    """
    node2score = {}
    with open(scores_file, 'r') as f:
        header = f.readline()
        if not header.startswith("node_index\t"):
            raise Exception("Scores file problem, wrong format")

        for line in f:
            split_line = line.rstrip().split("\t")
            if len(split_line) != 4:
                raise Exception(f"Invalid line in scores file: {line}")

            (idx, node, score, pvalue) = split_line
            try:
                score = float(score)
            except ValueError:
                raise Exception(f"Invalid score value in scores file: {score} is not a float")
            node2score[node] = score

    return(node2score)


def top_percent(list, x):
    """
    Finds the top x% of values in a list
    
    returns: 
    - top: list with top x% of values
    - rest: list with the rest of values
    """
    num_elements = round(len(list) * x/100)

    top = list[:num_elements]
    rest = list[num_elements:]
    
    return(top, rest)


def get_highest_scoring_nodes(node2score, network, protein2enrichment, threshold=10.0):
    """
    Get the top x% of highest-scoring genes (including causal) based on the provided
    node2score dictionary.

    arguments:
    - node2score: dict with key=node, value=score
    - network: networkx graph of the interactome
    - protein2enrichment: dict with key=protein, value=tissue enrichment
    - threshold: float, percentage threshold for selecting highest-scoring genes

    returns:
    - highest_scoring: list of highest-scoring genes (top x%)
    - low_scoring: list of low-scoring genes (the rest)
    """

    # # sort scores by descending order
    # results_sorted = sorted(node2score.keys(), key=lambda item: node2score[item], reverse=True)
    # scores_sorted = [node2score[node] for node in results_sorted]

    node2score_copy = node2score.copy()

    # get top x% of highest-scoring genes (including causal)
    # due to obsolete GTEx data, some genes are removed from the tissue enrichment analysis
    for protein in network:
        if protein not in protein2enrichment:
            del node2score_copy[protein]

    results_sorted = sorted(node2score_copy.keys(), key=lambda item: node2score_copy[item], reverse=True)
    highest_scoring, low_scoring = top_percent(results_sorted, threshold)

    return(highest_scoring, low_scoring)


def intersection(list_1, list_2):
    """
    Returns a list with common elements for two lists
    """
    return list(set(list_1) & set(list_2))


def format_pvalue(pvalue):
    """
    Returns a p-value in the scientific notation, eg. p < 0.05*
    """
    p = '{:.2e}'.format(pvalue)

    if pvalue <= 0.0001:
        return "p < 10^-4****"
    elif pvalue <= 0.001:
        return "p < 10^-3***"
    elif pvalue <= 0.01:
        return "p < 0.01**"
    elif pvalue <= 0.05:
        return "p < 0.05*"    
    return f"p = {p}"


def contingency_matrix(highest_scoring, low_scoring, tissue_enriched, non_tissue_enriched):
    """
    Compares the ratio of the highest scoring genes enriched in the expected tissues
    with the ratio of all genes enriched in that tissue using Fisher's exact test.
    Checks if the ratio of the highest scoring genes enriched in tissue is greater 
    than the ratio of all genes enriched in the tissue

    Note: alternative='greater' is the probability that a random table has x >= a, 
    where contingency_matrix = [[a, b], [c, d]]

    arguments:
    - highest_scoring: list of highest-scoring genes (top x%)
    - low_scoring: list of low-scoring genes (the rest)
    - tissue_enriched: list of genes enriched in the tissue
    - non_tissue_enriched: list of genes non-enriched in the tissue
    """
    x1 = len(intersection(highest_scoring, tissue_enriched))
    x2 = len(intersection(low_scoring, tissue_enriched))
    x3 = len(intersection(highest_scoring, non_tissue_enriched))
    x4 = len(intersection(low_scoring, non_tissue_enriched))

    contingency_matrix = [[x1, x2], [x3, x4]]

    _, pvalue = scipy.stats.fisher_exact(contingency_matrix, alternative='greater')

    return(format_pvalue(pvalue))


def comparison_matrix_row(method, highest_scoring, tissue_enriched, non_tissue_enriched, pvalue):
    """
    Creates a row (list) for comparison matrix containing:
    - method's name,
    - the number of the highest-scoring genes that are tissue-enriched,
    - the number of the highest-scoring genes that are non-tissue-enriched,
    - percentage of the highest-scoring tissue-enriched genes,
    - p-value for tissue enrichment
    """
    x1 = len(intersection(highest_scoring, tissue_enriched))
    x2 = len(intersection(highest_scoring, non_tissue_enriched))
    perc = round(x1 / (x1 + x2), 2)

    comparison_row = [x1, x2, perc, pvalue]

    return(comparison_row)


############################
########### MAIN ###########
############################

def main(network_file, uniprot_file, gtex_file, tissue, BFWalk_scores_file,
         multixrank_scores_file=None, netcore_scores_file=None,
         comparison_matrix_path="comparison_matrix.png",
         highest_scoring_threshold=10.0, enrichement_threshold=10.0,
         weighted=False, directed=False):

    if not (multixrank_scores_file or netcore_scores_file):
        logger.warning("Provide MultiXrank and/or NetCore scores files")
        raise Exception("No MultiXrank and/or NetCore scores files provided")

    logger.info(f"Parsing network {network_file}")
    (edge_list, node2idx, idx2node) = data_parser.parse_network(network_file, weighted, directed)

    # construct a networkx Graph from the edge list
    edge_list_no_weights = []
    for edge in edge_list:
        edge_list_no_weights.append((idx2node[edge[0]], idx2node[edge[1]]))
    interactome = networkx.from_edgelist(edge_list_no_weights)

    logger.info("Parsing Uniprot")
    uniprot2gene = parse_uniprot(uniprot_file)

    logger.info("Parsing GTEx")
    gene2enrichment = parse_expression_data(gtex_file, tissue)

    protein2enrichment = {}
    for protein in interactome:
        # due to obsolete GTEx data, some genes are removed from the tissue enrichment analysis
        if protein in uniprot2gene:
            gene = uniprot2gene[protein]
            if gene in gene2enrichment:
                protein2enrichment[protein] = gene2enrichment[gene]

    # get top x% of genes enriched in the tissue
    tissue_expr_sorted = sorted(protein2enrichment.keys(), key=lambda item: protein2enrichment[item], reverse=True)
    tissue_enriched, non_tissue_enriched = top_percent(tissue_expr_sorted, enrichement_threshold)
    logger.info(f"Selected {len(tissue_enriched)} tissue-enriched genes (removed {len(interactome) - len(protein2enrichment)} due to obsolete GTEx)")

    comparison_matrix = []
    row_labels = []
    logger.info("Parsing BFWalk scores")
    BFWalk_node2score = parse_BFWalk_scores(BFWalk_scores_file)
    (BFWalk_highest_scoring, BFWalk_low_scoring) = get_highest_scoring_nodes(BFWalk_node2score, interactome, protein2enrichment, highest_scoring_threshold)
    logger.info(f"BFWalk: Selected {len(BFWalk_highest_scoring)} highest-scoring genes")
    BFWalk_pvalue = contingency_matrix(BFWalk_highest_scoring, BFWalk_low_scoring, tissue_enriched, non_tissue_enriched)
    logger.info(f"BFWalk enrichment: {BFWalk_pvalue}")
    BFWalk_comparison_row = comparison_matrix_row("BFWalk", BFWalk_highest_scoring, tissue_enriched, non_tissue_enriched, BFWalk_pvalue)
    comparison_matrix.append(BFWalk_comparison_row)
    row_labels.append("BFWalk")

    if multixrank_scores_file:
        logger.info("Parsing MultiXrank scores")
        multixrank_node2score = parse_multixrank_scores(multixrank_scores_file)
        (multixrank_highest_scoring, multixrank_low_scoring) = get_highest_scoring_nodes(multixrank_node2score, interactome, protein2enrichment, highest_scoring_threshold)
        logger.info(f"MultiXrank: Selected {len(multixrank_highest_scoring)} highest-scoring genes")

        multixrank_pvalue = contingency_matrix(multixrank_highest_scoring, multixrank_low_scoring, tissue_enriched, non_tissue_enriched)
        logger.info(f"MultiXrank enrichment: {multixrank_pvalue}")

        multixrank_comparison_row = comparison_matrix_row("MultiXrank", multixrank_highest_scoring, tissue_enriched, non_tissue_enriched, multixrank_pvalue)
        comparison_matrix.append(multixrank_comparison_row)
        row_labels.append("MultiXrank")

    if netcore_scores_file:
        logger.info("Parsing NetCore scores")
        netcore_node2score = parse_netcore_scores(netcore_scores_file)
        # in netcore, some gene scores are missing, because they are not linked
        # to the main network component, we assign them scores=0,
        # so that they are included in the analysis
        for protein in node2idx:
            if protein not in netcore_node2score:
                netcore_node2score[protein] = 0.0
        (netcore_highest_scoring, netcore_low_scoring) = get_highest_scoring_nodes(netcore_node2score, interactome, protein2enrichment, highest_scoring_threshold)
        logger.info(f"NetCore: Selected {len(netcore_highest_scoring)} highest-scoring genes")
        netcore_pvalue = contingency_matrix(netcore_highest_scoring, netcore_low_scoring, tissue_enriched, non_tissue_enriched)
        logger.info(f"NetCore enrichment: {netcore_pvalue}")

        netcore_comparison_row = comparison_matrix_row("NetCore", netcore_highest_scoring, tissue_enriched, non_tissue_enriched, netcore_pvalue)
        comparison_matrix.append(netcore_comparison_row)
        row_labels.append("NetCore")

    # plot the comparison table
    fig, ax = matplotlib.pyplot.subplots(figsize=(10, 8))
    table = ax.table(comparison_matrix, 
                    rowLabels=row_labels,
                    colLabels=['tissue-enriched', 'non-tissue-enriched', "ratio", "p-value"],
                    loc="center")
    ax.axis('tight')
    ax.axis('off')
    table.set_fontsize(10)
    table.scale(1,2)

    comparison_matrix_path.parent.mkdir(parents=True, exist_ok=True)  # Path.parent of a bare filename returns Path("."), and mkdir on "."
    matplotlib.pyplot.savefig(comparison_matrix_path, dpi=500, bbox_inches='tight')


if __name__ == "__main__":
    script_name = os.path.basename(sys.argv[0])
    # configure logging, sub-modules will inherit this config
    logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    # set up logger: we want script name rather than 'root'
    logger = logging.getLogger(script_name)

    parser = argparse.ArgumentParser(
        prog=script_name,
        description="""
        Validation of the tissue enrichment of the highest scoring genes,
        using the BFWalk, MultiXrank and NetCore scores.
        """)
    parser.add_argument('--network',
                        help="Path to the network SIF file",
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--uniprot',
                        help="Path to the UniProt file created with uniprot_parser.py, with columns: " \
                        "Uniprot Primary AC, Uniprot Secondary AC(s), tax ID, gene name, synonyms",
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--gtex',
                        help="path to the TSV file with expression data from Expression Atlas",
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--tissue',
                        help="name of tissue of interest (must be one of the tissues in GTEx file)",
                        type=str,
                        required=True)
    parser.add_argument('--BFWalk_scores',
                        help="Path to the BFWalk scores file (TSV with header, columns: NODE, SCORE)",
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--multixrank_scores',
                        help="Path to the MultiXrank scores file (TSV with header, columns: multiplex, node, score)",
                        type=pathlib.Path,
                        required=False,
                        default=None)
    parser.add_argument('--netcore_scores',
                        help="Path to the NetCore scores file (TSV with header, columns: node_index, node, prop_weight, pvalue)",
                        type=pathlib.Path,
                        required=False,
                        default=None)
    parser.add_argument('--matrix',
                        help="Path to the output figure for the comparison matrix (default: comparison_matrix.png)",
                        type=pathlib.Path,
                        required=False,
                        default="comparison_matrix.png")
    parser.add_argument('--highest_scoring_threshold',
                        help="Threshold for the highest scoring nodes (default: 10.0%)",
                        type=float,
                        required=False,
                        default=10.0)
    parser.add_argument('--enrichement_threshold',
                        help="Threshold for tissue enrichment (default: 10.0%)",
                        type=float,
                        required=False,
                        default=10.0)
    parser.add_argument('--weighted',
                        help="Whether the network is weighted (default: False)",
                        action='store_true',
                        required=False)
    parser.add_argument('--directed',
                        help="Whether the network is directed (default: False)",
                        action='store_true',
                        required=False)

    args = parser.parse_args()

    try:
        main(args.network,
             args.uniprot,
             args.gtex,
             args.tissue,
             args.BFWalk_scores,
             multixrank_scores_file=args.multixrank_scores,
             netcore_scores_file=args.netcore_scores,
             comparison_matrix_path=args.matrix,
             highest_scoring_threshold=args.highest_scoring_threshold,
             enrichement_threshold=args.enrichement_threshold,
             weighted=args.weighted,
             directed=args.directed)

    except Exception as e:
        # details on the issue should be in the exception name, print to stderr and die
        sys.stderr.write(f"ERROR in {script_name} : {repr(e)}\n")
        sys.exit(1)
