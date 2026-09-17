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

import statistics
import numpy
import networkx
import matplotlib.pyplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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


def calculate_rank_difference(BFWalk_node2rank, other_method_node2rank, network):
    ranks_BFWalk = []
    ranks_other = []
    node_degrees = []

    for node in BFWalk_node2rank:
        ranks_BFWalk.append(BFWalk_node2rank[node])
        ranks_other.append(other_method_node2rank[node])
        node_degrees.append(network.degree(node))
    assert len(ranks_BFWalk) == len(ranks_other) == len(node_degrees)

    rank_diff = []
    rank_diff_abs = []
    negative_rank_degrees = []
    positive_rank_degrees = []

    for i in range(len(ranks_BFWalk)):
        diff = ranks_BFWalk[i] - ranks_other[i]
        rank_diff.append(diff)
        rank_diff_abs.append(abs(diff))

        if diff < 0:
            negative_rank_degrees.append(node_degrees[i])
        else:
            positive_rank_degrees.append(node_degrees[i])

    logger.info(f"absolute rank difference mean: {round(statistics.mean(rank_diff_abs))}, median: {round(statistics.median(rank_diff_abs))}")
    logger.info(f"negative rank degree mean: {round(statistics.mean(negative_rank_degrees))}, median: {round(statistics.median(negative_rank_degrees))}")
    logger.info(f"positive rank degree mean: {round(statistics.mean(positive_rank_degrees))}, median: {round(statistics.median(positive_rank_degrees))}")

    return(rank_diff, negative_rank_degrees, positive_rank_degrees, node_degrees)


def plot_rankVsDeg(rank_diff, negative_rank_degrees, positive_rank_degrees, node_degrees, interactome, other_method, out):
    fig, ax = matplotlib.pyplot.subplots(figsize=(7, 6))

    ax.scatter(x=rank_diff, y=node_degrees, s=1, alpha=0.8, zorder=3)
    ax.set_xlim(-(len(interactome)), len(interactome))
    ax.axvline(0, color='grey', linestyle='--', linewidth=1.5, zorder=2)

    ax.set_xlabel(f"Rank difference (BFWalk rank - {other_method} rank)", fontsize=12)
    ax.set_ylabel("Node degree", fontsize=12)

    ax.grid(True, linestyle='--', which='major',
                    color='grey', alpha=0.2)

    inset = inset_axes(ax, width="30%", height="30%", loc=2, borderpad=4)

    inset.boxplot([negative_rank_degrees, positive_rank_degrees],
                flierprops=dict(marker='o', markersize=2, markerfacecolor='black'))
    inset.set_xticklabels(["Rank diff < 0", "Rank diff > 0"], fontsize=8, rotation=45)
    inset.set_ylabel("Node degree", fontsize=8)
    inset.tick_params(axis='y', labelsize=8)

    inset.grid(True, linestyle="--", alpha=0.2)

    matplotlib.pyplot.savefig(out, dpi=500)


def main(network_file, phenotypes, BFWalk_out_dir, multixrank_out_dir=None, netcore_out_dir=None,
         rankVsDeg_dir="./", weighted=False, directed=False):

    if not (multixrank_out_dir or netcore_out_dir):
        logger.warning("Provide MultiXrank and/or NetCore output directories")
        raise Exception("No MultiXrank and/or NetCore output directories provided")

    rankVsDeg_dir.parent.mkdir(parents=True, exist_ok=True)  # Path.parent of a bare filename returns Path("."), and mkdir on "."

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

    # dicts to store node-rank pairs for all phenotypes combined
    BFWalk_node2rank = {}
    multixrank_node2rank = {}
    netcore_node2rank = {}

    for phenotype in phenotypes:
        logger.info(f"Phenotype: {phenotype}")
        logger.info("Parsing BFWalk ranks")
        # search all subdirectories of the BFWalk output directory for "ranks_LOO.tsv"
        BFWalk_ranks_file = os.path.join(BFWalk_out_dir, phenotype, "ranks_LOO.tsv")
        BFWalk_node2rank_pheno = parse_ranks(BFWalk_ranks_file)
        BFWalk_node2rank.update(BFWalk_node2rank_pheno)

        if multixrank_out_dir:
            logger.info("Parsing MultiXrank ranks")
            multixrank_ranks_file = os.path.join(multixrank_out_dir, phenotype, "ranks_LOO.tsv")
            multixrank_node2rank_pheno = parse_ranks(multixrank_ranks_file)
            multixrank_node2rank.update(multixrank_node2rank_pheno)
            logger.info(f"BFWalk vs MultiXrank ({phenotype}):")
            (rank_diff, negative_rank_degrees, positive_rank_degrees, node_degrees) = calculate_rank_difference(BFWalk_node2rank_pheno,
                                                                                                                multixrank_node2rank_pheno,
                                                                                                                interactome)

        if netcore_out_dir:
            logger.info("Parsing NetCore scores")
            netcore_LOO_dir = os.path.join(netcore_out_dir, phenotype)
            netcore_node2rank_pheno = netcore_scores_to_ranks(BFWalk_node2rank_pheno.keys(), netcore_LOO_dir, len(node2idx))
            netcore_node2rank.update(netcore_node2rank_pheno)
            logger.info(f"BFWalk vs NetCore ({phenotype}):")
            (rank_diff, negative_rank_degrees, positive_rank_degrees, node_degrees) = calculate_rank_difference(BFWalk_node2rank_pheno,
                                                                                                                netcore_node2rank_pheno,
                                                                                                                interactome)

    if multixrank_out_dir:
        assert len(BFWalk_node2rank) == len(multixrank_node2rank), "BFWalk and MultiXrank ranks files have different number of left-out nodes"
        logger.info(f"BFWalk vs MultiXrank (all phenotypes):")
        (rank_diff, negative_rank_degrees, positive_rank_degrees, node_degrees) = calculate_rank_difference(BFWalk_node2rank,
                                                                                                            multixrank_node2rank,
                                                                                                            interactome)
        rankVsDeg_path = os.path.join(rankVsDeg_dir, "all_rank_vs_deg_BFWalk_vs_RWR.png")
        plot_rankVsDeg(rank_diff, negative_rank_degrees, positive_rank_degrees, node_degrees, interactome, "RWR", rankVsDeg_path)

    if netcore_out_dir:
        assert len(BFWalk_node2rank) == len(netcore_node2rank), "BFWalk and NetCore ranks files have different number of left-out nodes"
        logger.info(f"BFWalk vs NetCore (all phenotypes):")
        (rank_diff, negative_rank_degrees, positive_rank_degrees, node_degrees) = calculate_rank_difference(BFWalk_node2rank,
                                                                                                            netcore_node2rank,
                                                                                                            interactome)
        rankVsDeg_path = os.path.join(rankVsDeg_dir, "all_rank_vs_deg_BFWalk_vs_NetCore.png")
        plot_rankVsDeg(rank_diff, negative_rank_degrees, positive_rank_degrees, node_degrees, interactome, "NetCore", rankVsDeg_path)

    logger.info(f"Found {len(BFWalk_node2rank)} left-out nodes")


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
        Validation of BFWalk ranks by plotting the left-out genes ranks vs their degrees.
        It assumes that for each method the output directory provided as input
        contain one subdirectory per phenotype and searches all these subdirectories for ranks files
        (or scores files for NetCore).
        """)
    parser.add_argument('--network',
                        help="Path to the network SIF file",
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--phenotypes',
                        nargs='+',
                        help="List of phenotypes",
                        required=True)
    parser.add_argument('--BFWalk_out_dir',
                        help="Path to the BFWalk output directory containing one subdirectory per phenotype (given by --phenotypes).",
                        type=pathlib.Path,
                        required=True,)
    parser.add_argument('--multixrank_out_dir',
                        help="Path to the MultiXrank output directory",
                        type=pathlib.Path,
                        required=False,
                        default=None)
    parser.add_argument('--netcore_out_dir',
                        help="Path to the NetCore output directory",
                        type=pathlib.Path,
                        required=False,
                        default=None)
    parser.add_argument('--rankVsDeg',
                        help="Path to the output directory for ranks vs degree plots (default: ./)",
                        type=pathlib.Path,
                        required=False,
                        default="./")
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
             args.phenotypes,
             args.BFWalk_out_dir,
             multixrank_out_dir=args.multixrank_out_dir,
             netcore_out_dir=args.netcore_out_dir,
             rankVsDeg_dir=args.rankVsDeg,
             weighted=args.weighted,
             directed=args.directed)

    except Exception as e:
        # details on the issue should be in the exception name, print to stderr and die
        sys.stderr.write(f"ERROR in {script_name} : {repr(e)}\n")
        sys.exit(1)