############################################################################################
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

import scipy
import statistics
import numpy
import networkx
import matplotlib.pyplot

sys.path.append("/home/kubicaj/Software/BFWalk")
import data_parser

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)

matplotlib.pyplot.set_loglevel('warning')


def parse_ranks(ranks_file):
    """
    Parses a TSV ranks file (with header) with two columns: NODE, RANK

    Returns:
    - node2rank: dict with key=node, value=rank
    """
    node2rank = {}
    
    with open(ranks_file, 'r') as f:
        header = f.readline()
        if not header.startswith("NODE\t"):
            raise Exception("Ranks file problem, wrong format")

        for line in f:
            split_line = line.rstrip().split("\t")
            if len(split_line) != 2:
                raise Exception(f"Invalid line in ranks file: {line}")

            (node, rank) = split_line
            try:
                rank = int(rank)
            except ValueError:
                raise Exception(f"Invalid rank value in ranks file: {rank} is not an integer")
            node2rank[node] = rank

    return(node2rank)


def netcore_scores_to_ranks(left_out, netcore_LOO_dir, n_nodes):
    """
    Iterates over the list of left-out nodes,
    for each left-out node opens a leave-one-out file (with nodes sorted by descending score,
    so the first node has rank=1, the second node has rank=2, etc),
    finds the rank of the left-out node;
    for left-out nodes for which scores are missing (not in the main network component), we assign scores=0.

    arguments:
    - left_out: list of left-out nodes
    - netcore_LOO_dir: directory with NetCore scores files for left-out nodes,
        eg. for each left-out node scores for all nodes in the interactome are in:
        netcore-output/{PHENOTYPE}/output_{LEFT-OUT-NODE}/random_walk_weights.txt
    - n_nodes: total number of nodes in the network
    
    returns:
    - node2rank: dict with key=left-out node, value=rank
    """
    node2rank = {}

    for left_out_node in left_out:
        netcore_LOO_file = os.path.join(netcore_LOO_dir, f"output_{left_out_node}", "random_walk_weights.txt")
        if not os.path.isfile(netcore_LOO_file):
            raise Exception(f"NetCore scores file not found for left-out node {left_out_node}: {netcore_LOO_file}")

        with open(netcore_LOO_file, 'r') as f:
            header = f.readline()
            if not header.startswith("node_index\t"):
                raise Exception(f"NetCore scores file problem for left-out node {left_out_node}, wrong format: {netcore_LOO_file}")

            rank = 1
            found = False
            for line in f:
                split_line = line.rstrip().split("\t")
                if len(split_line) != 4:
                    raise Exception(f"Invalid line in NetCore scores file for left-out node {left_out_node}: {line}")

                (idx, node, score, pvalue) = split_line
                if node == left_out_node:
                    node2rank[left_out_node] = rank
                    found = True
                    break
                rank += 1

            if not found:
                # If the left-out node is not found in the NetCore scores file, assign score=0 (rank==n_nodes)
                node2rank[left_out_node] = n_nodes

    return(node2rank)


def generate_random_ranks(n_left_out, network_size):
    """
    Generate ranks with a random classifier

    arguments:
    - n_left_out: number of left-out genes
    - network_size: total number of nodes in the network

    returns:
    - random_ranks: list of average ranks for the left-out genes across n iterations
    """

    # use the integers method of a Generator instance instead of numpy.random.randint
    # https://numpy.org/doc/2.1/reference/random/generated/numpy.random.randint.html
    rng = numpy.random.default_rng()

    # 1D numpy array for left-out gene ranks
    # size=len(left-out genes)
    random_ranks = numpy.zeros(n_left_out, dtype=numpy.uint64)

    n = 10000  # 10,000 iterations for convergence
    for i in range(n):
        # 1D numpy array for left-out gene ranks
        # size=len(left-out genes)
        tmp = rng.integers(low=1, high=network_size, endpoint = True, size=random_ranks.shape, dtype=numpy.uint64)
        tmp = numpy.sort(tmp)
        random_ranks += tmp

    for j in range(len(random_ranks)):
        avg_rank = random_ranks[j] / n
        random_ranks[j] = avg_rank
    
    return(random_ranks)


def ranks_to_curve(ranks, network_size):
    """
    Computes a curve and AUC of rank vs cumulative distributions of left-out nodes

    arguments:
    - ranks: list of ranks for the left-out nodes
    - network_size: total number of nodes in the network

    returns:
    - curve: list of cumulative distributions of left-out nodes
    - AUC: normalized area under the curve
    """
    curve = []
    
    for i in range(1, network_size + 1):
        curve.append(sum(rank <= i for rank in ranks))

    AUC = 0.0
    for i in range(1, len(curve)):
        AUC += (curve[i] + curve[i-1]) / 2

    max_possible = (len(curve) - 1) * len(ranks)  # maximum area: all left-out ranked at the top (rank=1)
    AUC_norm = AUC / max_possible

    return(curve, AUC_norm)


def plot_CDF(BFWalk_curve, BFWalk_AUC, random_curve, random_AUC, network_size, pheno="", out="CDF.png",
             multixrank_curve=None, multixrank_AUC=None, netcore_curve=None, netcore_AUC=None):
    x = range(network_size)
    matplotlib.pyplot.plot(x, BFWalk_curve, label="BFWalk (AUC={:.3f})".format(BFWalk_AUC), color="#D81B60")
    if multixrank_curve:
        matplotlib.pyplot.plot(x, multixrank_curve, label="MultiXrank (AUC={:.3f})".format(multixrank_AUC), color="#FFC107")
    if netcore_curve:
        matplotlib.pyplot.plot(x, netcore_curve, label="NetCore (AUC={:.3f})".format(netcore_AUC), color="#1E88E5")
    matplotlib.pyplot.plot(x, random_curve, label="random classifier (AUC={:.3f})".format(random_AUC), color="#004D40")

    matplotlib.pyplot.xlabel("rank x", fontsize=12)
    matplotlib.pyplot.ylabel("Number of left-out genes where rank <= x", fontsize=12)
    matplotlib.pyplot.title(pheno)
    matplotlib.pyplot.xticks(fontsize=11)
    matplotlib.pyplot.yticks(fontsize=11)
    matplotlib.pyplot.legend(loc='lower right', fontsize=12)

    matplotlib.pyplot.savefig(out, dpi=500)

    logger.info(f"CDF curve saved to {out}")


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


def main(network_file, BFWalk_ranks_file, multixrank_ranks_file=None, netcore_LOO_dir=None,
         cdf_path=None, pheno="", weighted=False, directed=False):
    
    logger.info(f"Parsing network {network_file}")
    (edge_list, node2idx, idx2node) = data_parser.parse_network(network_file, weighted, directed)

    # construct a networkx Graph from the edge list
    edge_list_no_weights = []
    for edge in edge_list:
        edge_list_no_weights.append((idx2node[edge[0]], idx2node[edge[1]]))
    interactome = networkx.from_edgelist(edge_list_no_weights)

    # calculate the mean and median degree of nodes in the network
    network_degrees = []
    for node in interactome.nodes():
        network_degrees.append(interactome.degree(node))
    logger.info(f"node degree mean: {round(statistics.mean(network_degrees))}, median: {round(statistics.median(network_degrees))}")

    logger.info("Parsing BFWalk ranks")
    BFWalk_node2rank = parse_ranks(BFWalk_ranks_file)
    if multixrank_ranks_file:
        logger.info("Parsing MultiXrank ranks")
        multixrank_node2rank = parse_ranks(multixrank_ranks_file)
        assert len(BFWalk_node2rank) == len(multixrank_node2rank), "BFWalk and MultiXrank ranks files have different number of left-out nodes"
    if netcore_LOO_dir:
        logger.info("Parsing NetCore scores")
        netcore_node2rank = netcore_scores_to_ranks(BFWalk_node2rank.keys(), netcore_LOO_dir, len(node2idx))
        assert len(BFWalk_node2rank) == len(netcore_node2rank), "BFWalk and NetCore ranks files have different number of left-out nodes"
    logger.info(f"Found {len(BFWalk_node2rank)} left-out nodes")

    # calculate the mean and median degree of the left-out genes (seeds)
    seeds_degrees = []
    for protein in BFWalk_node2rank:
        seeds_degrees.append(interactome.degree(protein))
    logger.info(f"left-out degree mean: {round(statistics.mean(seeds_degrees))}, median: {round(statistics.median(seeds_degrees))}")

    logger.info("Calculating CDF curves and AUCs")
    (BFWalk_curve, BFWalk_AUC) = ranks_to_curve(BFWalk_node2rank.values(), len(node2idx))
    logger.info(f"AUC: {BFWalk_AUC:.3f} (BFWalk)")

    optional_kwargs = {}  # dict for plotting
    if multixrank_ranks_file is not None:
        (multixrank_curve, multixrank_AUC) = ranks_to_curve(multixrank_node2rank.values(), len(node2idx))
        logger.info(f"AUC: {multixrank_AUC:.3f} (MultiXrank)")
        optional_kwargs.update({'multixrank_curve': multixrank_curve, 'multixrank_AUC': multixrank_AUC})

    if netcore_LOO_dir is not None:
        (netcore_curve, netcore_AUC) = ranks_to_curve(netcore_node2rank.values(), len(node2idx))
        logger.info(f"AUC: {netcore_AUC:.3f} (NetCore)")
        optional_kwargs.update({'netcore_curve': netcore_curve, 'netcore_AUC': netcore_AUC})

    random_ranks = generate_random_ranks(len(BFWalk_node2rank), len(node2idx))
    (random_curve, random_AUC) = ranks_to_curve(random_ranks, len(node2idx))

    # "Are BFWalk ranks higher than MultiXrank ranks?"
    p_value_BFWalk_vs_multixrank = scipy.stats.wilcoxon(list(multixrank_node2rank.values()), list(BFWalk_node2rank.values()), alternative="greater").pvalue
    logger.info(f"p-value (ranks BFWalk vs MultiXrank): {format_pvalue(p_value_BFWalk_vs_multixrank)}")

    # "Are BFWalk ranks higher than NetCore ranks?"
    p_value_BFWalk_vs_netcore = scipy.stats.wilcoxon(list(netcore_node2rank.values()), list(BFWalk_node2rank.values()), alternative="greater").pvalue
    logger.info(f"p-value (ranks BFWalk vs NetCore): {format_pvalue(p_value_BFWalk_vs_netcore)}")

    phenotype=""
    if pheno == "MMAF":
        phenotype = "MMAF"
    elif pheno == "HYPCARD":
        phenotype = "hypertrophic cardiomyopathy"
    elif pheno == "CKD":
        phenotype = "chronic kidney disease"
    elif pheno == "DYSCHROM":
        phenotype = "dyschromatopsia"

    cdf_path.parent.mkdir(parents=True, exist_ok=True)  # Path.parent of a bare filename returns Path("."), and mkdir on "."
    plot_CDF(BFWalk_curve=BFWalk_curve,
             BFWalk_AUC=BFWalk_AUC,
             random_curve=random_curve,
             random_AUC=random_AUC,
             network_size=len(interactome.nodes()),
             pheno=phenotype,
             out=cdf_path,
             **optional_kwargs)


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
        Validation of BFWalk ranks by plotting the cumulative distribution function (CDF) of left-out genes ranks.
        The CDF curve shows the number of left-out genes with rank <= x for each rank x.
        The area under the curve (AUC) is also calculated and shown in the legend.
        """)
    parser.add_argument('--network',
                        help="Path to the network SIF file",
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--BFWalk_ranks',
                        help="Path to the BFWalk ranks file (TSV with header, columns: NODE, RANK)",
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--multixrank_ranks',
                        help="Path to the MultiXrank ranks file (TSV with header, columns: NODE, RANK)",
                        type=pathlib.Path,
                        required=False,
                        default=None)
    parser.add_argument('--netcore_LOO_dir',
                        help="Path to the directory with NetCore scores files for left-out nodes, " \
                        "eg. for each left-out node scores for all nodes in the interactome are in: " \
                        "netcore-output/{PHENOTYPE}/output_{LEFT-OUT-NODE}/random_walk_weights.txt",
                        type=pathlib.Path,
                        required=False,
                        default=None)
    parser.add_argument('--cdf',
                        help="Path to the output figure for the CDF curve (default: CDF.png)",
                        type=pathlib.Path,
                        required=False,
                        default="CDF.png")
    parser.add_argument('--pheno',
                        help="Name of the phenotype (for plotting only)",
                        type=str,
                        required=False,
                        default=""),
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
             args.BFWalk_ranks,
             multixrank_ranks_file=args.multixrank_ranks,
             netcore_LOO_dir=args.netcore_LOO_dir,
             cdf_path=args.cdf,
             pheno=args.pheno,
             weighted=args.weighted,
             directed=args.directed)

    except Exception as e:
        # details on the issue should be in the exception name, print to stderr and die
        sys.stderr.write(f"ERROR in {script_name} : {repr(e)}\n")
        sys.exit(1)