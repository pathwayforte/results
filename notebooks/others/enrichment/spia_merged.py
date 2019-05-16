# -*- coding: utf-8 -*-

"""Merge SPIA files from the three resources to build the SPIA merged data set"""

import logging
import os

import click
import pybel
from bio2bel_reactome import Manager as ReactomeManager
from pybel import from_pickle
from pybel.struct.mutation import collapse_all_variants, collapse_to_genes
from pybel_tools.analysis.spia import bel_to_spia_matrices
from tqdm import tqdm

from pathme.constants import KEGG, REACTOME, WIKIPATHWAYS
from pathme.constants import KEGG_BEL, REACTOME_BEL, WIKIPATHWAYS_BEL, SPIA_DIR
from pathme.export_utils import get_all_pickles, yield_all_children, spia_matrices_to_excel
from pathme.normalize_names import normalize_graph_names
from pathme.pybel_utils import flatten_complex_nodes
from pathway_forte.mappings import get_equivalent_mappings_dict

logger = logging.getLogger(__name__)


def get_equivalent_pathways(
        mappings,
        reactome_manager,
        kegg_dir,
        reactome_dir,
        wikipathways_dir
):
    """Get equivalent pathways and merge them"""
    graphs_to_merge = list()
    merged_pathways = list()

    for resource, pathway in mappings:

        if resource == KEGG:
            graph = get_kegg_graph(os.path.join(kegg_dir, f"{pathway}.pickle"))

        elif resource == REACTOME:
            graph = get_reactome_graph(reactome_manager, reactome_dir, f"{pathway}.pickle")

        elif resource == WIKIPATHWAYS:
            graph = get_wp_graph(os.path.join(wikipathways_dir, f"{pathway}.pickle"))

        else:
            raise ValueError(f"{resource} not found")

        graphs_to_merge.append(graph)
        merged_pathways.append((resource, pathway))

    merged_graph = pybel.union(graphs_to_merge)

    return merged_graph, merged_pathways


def get_kegg_graph(file):
    pathway_graph = from_pickle(file)
    normalize_graph_names(pathway_graph, KEGG)
    return pathway_graph


def get_wp_graph(file):
    pathway_graph = from_pickle(file)
    normalize_graph_names(pathway_graph, WIKIPATHWAYS)
    return pathway_graph


def get_reactome_graph(reactome_manager, reactome_dir, file):
    # Load BELGraph
    pathway_graph = from_pickle(os.path.join(reactome_dir, file))

    # Check if pathway has children to build the merge graph
    pathway_id = file.strip('.pickle')

    # Look up in Bio2BEL Reactome
    pathway = reactome_manager.get_pathway_by_id(pathway_id)

    # Log if it is not present
    if not pathway:
        logger.warning(f'{pathway_id} not found in database')

    # Check if there are children and merge them on the fly
    for child in yield_all_children(pathway):

        child_file_path = os.path.join(reactome_dir, f"{child.resource_id}.pickle")
        if not os.path.exists(child_file_path):
            logger.warning(f'{child.resource_id} pickle does not exist')
            continue

        # Load the pickle and union it
        child_graph = pybel.from_pickle(child_file_path)
        pathway_graph += child_graph

    # Normalize graph names
    normalize_graph_names(pathway_graph, REACTOME)

    return pathway_graph


def merge_spia_files(
        mapping_dict: dict,
        kegg_dir: str,
        reactome_dir: str,
        wikipathways_dir: str,
        output: str = SPIA_DIR

):
    """Merges PathMe pickles to generate the SPIA merged dataset.

    :param kegg_path: directory to KEGG pickles
    :param reactome_path: directory to Reactome pickles
    :param wikipathways_path: directory to WikiPathways pickles
    :param output: output directory
    """
    kegg_pickles, reactome_pickles, wp_pickles = get_all_pickles(kegg_dir, reactome_dir, wikipathways_dir)

    all_pickles = kegg_pickles + reactome_pickles + wp_pickles

    click.echo(f'A total of {len(all_pickles)} files will be exported')

    iterator = tqdm(all_pickles, desc='Exporting SPIA excel files')

    # Call Reactome manager and check that is populated
    reactome_manager = ReactomeManager()

    if not reactome_manager.is_populated():
        logger.warning('Reactome Manager is not populated')

    # Keep track of which pathways have been already merged to skip them
    pathways_already_merged = set()

    # Load each pickle and export it as excel file
    for file in iterator:
        if not file.endswith('.pickle'):
            continue

        if file in kegg_pickles:

            pathway_graph = get_kegg_graph(os.path.join(kegg_dir, file))

            pathway_id = file.strip("_unflatten.pickle")

            # Skip if the pathway has been already merged
            if (KEGG, pathway_id) in pathways_already_merged:
                continue

            if (KEGG, pathway_id) in mapping_dict:
                graph_equivalents, merged_pathways = get_equivalent_pathways(
                    mapping_dict[(KEGG, pathway_id)], reactome_manager, kegg_dir, reactome_dir, wikipathways_dir
                )

                # Merge graph with equivalents
                pathway_graph += graph_equivalents

                for pathway in merged_pathways:
                    pathway_id += f"|{pathway[1]}"
                    pathways_already_merged.add(pathway)

                print(pathway_id)

        elif file in reactome_pickles:
            pathway_graph = get_reactome_graph(reactome_manager, reactome_dir, file)

            pathway_id = file.strip(".pickle")

            # Skip if the pathway has been already merged
            if (REACTOME, pathway_id) in pathways_already_merged:
                continue

            if (REACTOME, pathway_id) in mapping_dict:
                graph_equivalents, merged_pathways = get_equivalent_pathways(
                    mapping_dict[(REACTOME, pathway_id)], reactome_manager, kegg_dir, reactome_dir, wikipathways_dir
                )

                # Merge graph with equivalents
                pathway_graph += graph_equivalents

                for pathway in merged_pathways:
                    pathway_id += f"|{pathway[1]}"
                    pathways_already_merged.add(pathway)

                print(pathway_id)


        elif file in wp_pickles:
            pathway_graph = get_wp_graph(os.path.join(wikipathways_dir, file))

            pathway_id = file.strip(".pickle")

            # Skip if the pathway has been already merged
            if (WIKIPATHWAYS, pathway_id) in pathways_already_merged:
                continue

            if (WIKIPATHWAYS, pathway_id) in mapping_dict:
                graph_equivalents, merged_pathways = get_equivalent_pathways(
                    mapping_dict[(WIKIPATHWAYS, pathway_id)], reactome_manager, kegg_dir, reactome_dir, wikipathways_dir
                )

                # Merge graph with equivalents
                pathway_graph += graph_equivalents

                for pathway in merged_pathways:
                    pathway_id += f"|{pathway[1]}"
                    pathways_already_merged.add(pathway)

                print(pathway_id)

        else:
            logger.warning(f'Unknown pickle file: {file}')
            continue

        # Explode complex nodes
        flatten_complex_nodes(pathway_graph)

        # Collapse nodes
        collapse_all_variants(pathway_graph)
        collapse_to_genes(pathway_graph)

        spia_matrices = bel_to_spia_matrices(pathway_graph)

        output_file = os.path.join(output, f"{pathway_id}.xlsx")

        if os.path.isfile(output_file):
            continue

        # Export excel file representing the connectivity matrix of the BEL Graph
        spia_matrices_to_excel(spia_matrices, output_file)


if __name__ == "__main__":
    # Mapping dictionary
    compath_equivalent_mappings = get_equivalent_mappings_dict()
    merge_spia_files(compath_equivalent_mappings, KEGG_BEL, REACTOME_BEL, WIKIPATHWAYS_BEL)
