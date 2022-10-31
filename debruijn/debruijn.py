#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    #fastq_file = "eva71_two_reads.fq"
    with open(fastq_file, "r") as f1:
        for i in f1:
            yield next(f1).strip()
            next(f1)
            next(f1)



def cut_kmer(read, kmer_size):
    for i in range(len(read)):
        if i + kmer_size <= len(read):

            yield read[i:i + kmer_size]



def build_kmer_dict(fastq_file, kmer_size):
    kmer_dico = {}
    for i in read_fastq(fastq_file):
        for j in cut_kmer(i, kmer_size):
            if j in kmer_dico:
                kmer_dico[j] += 1
            else:
                kmer_dico[j] = kmer_dico.get(j, 0) + 1
    return kmer_dico



def build_graph(kmer_dico):
    kmer_graph = nx.DiGraph()
    for key, value in kmer_dico.items():
        kmer_graph.add_edge(key[:-1], key[1:], weight=value)
    return kmer_graph



def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for node in path_list:
        if not delete_entry_node and not delete_sink_node:  # if both False
            graph.remove_nodes_from(node[1:-1])  # remove first and last node
        elif delete_entry_node and not delete_sink_node:
            graph.remove_nodes_from(node[:-1])  # remove first node
        elif not delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(node[1:])  # remove last node
        else:
            graph.remove_nodes_from(node)  # if both True
    return graph


def std(data):
    return statistics.stdev(data)



def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    best_length_p = []
    best_weight_p = []
    for comp in range(len(path_list)):  # We go through the list of paths
        if weight_avg_list[comp] == max(weight_avg_list):
            best_weight_p.append(path_list[comp])  # Stock path with max weight

        for comp_bis in range(len(best_weight_p)):
            len_max = len(best_weight_p[comp_bis])
            if len_max > 0:
                best_length_p.append(best_weight_p[comp_bis])

    best_path = random.choice(best_length_p)
    print(random.choice(best_length_p))  # Without this line, test not passed

    path_list.pop(path_list.index(best_path))
    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph




def path_average_weight(graph, path):
    avg = 0
    len_path = len(path) - 1  # -1 avoid index of out range
    for comp in range(len_path):
        avg += graph[path[comp]][path[comp + 1]]["weight"]
    return avg / len_path

def solve_bubble(graph, ancestor_node, descendant_node):
    path_length = []
    weight_avg_list = []
    path_list = list(nx.all_simple_paths(graph, source=ancestor_node, target=descendant_node))
    for path in path_list:
        path_length.append(len(path))
        weight_avg_list.append(path_average_weight(graph, path))
    graph = select_best_path(graph, path_list, path_length, weight_avg_list)
    return graph



def simplify_bubbles(graph):
    ancestors = []
    descendants = []
    for node in graph.nodes:
        successor_list = list(graph.successors(node))  # collect a list of successors of node
        predecessor_list = list(graph.predecessors(node))  # same with predecessors
        if len(successor_list) > 1:  # if successor -> ancestor
            ancestors.append(node)
        if len(predecessor_list) > 1:  # if predecessor -> descendant
            descendants.append(node)

    for comp_anc in range(len(ancestors)):
        for comp_desc in range(len(descendants)):
            if ancestors[comp_anc] != descendants[comp_desc]:  # if bubble
                graph = solve_bubble(graph, ancestors[comp_anc], descendants[comp_desc])
    return graph



def solve_entry_tips(graph, starting_nodes):
    lst_path = []
    wei_path = []
    len_path = []
    for node in starting_nodes:
        for descendant in nx.descendants(graph,
                                         node):  # descendants(G, source) Returns all nodes reachable from source in G.
            predecessor = list(graph.predecessors(descendant))
            if len(predecessor) > 1:
                for path in nx.all_simple_paths(graph, node, descendant):
                    lst_path.append(path)
                    len_path.append(len(path))
                    wei_path.append(path_average_weight(graph, path))
    graph = select_best_path(graph, lst_path, len_path, wei_path, True, False)
    return (graph)


def solve_out_tips(kmer_graph, ending_nodes):
    """lst_path = []
    wei_path = []
    len_path = []
    for node in ending_nodes:
        for ancestor in nx.ancestors(kmer_graph,
                                     node):  # ancestors(G, source) Returns all nodes having a path to source in G.
            successor = list(kmer_graph.successors(ancestor))
            if len(successor) > 1:
                for path in nx.all_simple_paths(kmer_graph, node, ancestor):
                    lst_path.append(path)
                    len_path.append(len(path))
                    wei_path.append(path_average_weight(kmer_graph, path))
    graph = select_best_path(kmer_graph, lst_path, len_path, wei_path, False, True)
    return (kmer_graph)"""
  path_list = []
    path_length = []
    weight_avg_list = []
    for node in kmer_graph.nodes:
        successors_list = list(kmer_graph.successors(node))
        if len(successors_list) > 1:
            for end in ending_nodes:
                for path in nx.all_simple_paths(kmer_graph, node, end):
                    path_list.append(path)

    if len(path_list) != 0:
        for path in path_list:
            path_length.append(len(path))
            weight_avg_list.append(path_average_weight(kmer_graph, path))
        graph = select_best_path(kmer_graph, path_list, path_length, weight_avg_list,
                                 delete_entry_node=False, delete_sink_node=True)
    return kmer_graph




def get_starting_nodes(kmer_graph):
    start = []
    for node in kmer_graph.nodes:
        if not list(kmer_graph.predecessors(node)):
            start.append(node)
    return start


def get_sink_nodes(kmer_graph):
   """ sink = []
    for node in kmer_graph.nodes:
        successor = list(kmer_graph.successors(node))
        if not successor:
            sink.append(node)
        return sink"""
   sink_nodes = []
   for node in kmer_graph.nodes:
       successor = kmer_graph.successors(node)  # Returns an iterator over successor nodes of n
       if not list(successor):  # if empty list
           sink_nodes.append(node)
   return sink_nodes


def get_contigs(graph, starting_nodes, ending_nodes):
    contigs_list = []
    for input_node in starting_nodes:
        for output_node in ending_nodes:
            for path in nx.all_simple_paths(graph, source=input_node, target=output_node):
                contig = path[0]  # first node
                for comp in path[1:]:  # path without first element
                    contig += comp[1:]  # add nodes to contig
                contigs_list.append((contig, len(contig)))
    return contigs_list



def save_contigs(contigs_list, output_file):
    file = open(output_file, "w")
    index = 0
    for contig, size_contig in contigs_list:
        file.write(">contig_" + str(index) + " len=" + str(size_contig) + "\n")
        file.write(fill(contig, width=80) + "\n")
        index += 1



def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Lecture du fichier et construction du graphe:
    dic_kmer = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(dic_kmer)
    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    # Résolution des bulles:
    simplify_bubbles = simplify_bubbles(graph)
    # Resolution des points d'entrée et de sortie
    entry_tips = solve_entry_tips(simplify_bubbles, starting_nodes)
    out_tips = solve_out_tips(entry_tips, sink_nodes)
    # Ecriture des contigs:
    contigs = get_contigs(out_tips, starting_nodes, sink_nodes)
    save_contigs(contigs, args.output_file)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
