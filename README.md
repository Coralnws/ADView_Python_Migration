# ADViewpy

## Introduction

**ADViewpy** project is a Python port of ADView system, originally developed by Professor Zipeng Liu. ADViewpy aims to provide visual comparison of phylogenetic trees, usually between one reference tree and a collection of hundreds of trees.

- **ADView System by Professor Zipeng Liu:** [GitHub](https://github.com/zipengliu) | [ADView Website](http://phylo-adview.site/)

The alpha version of this tool has been uploaded to PyPI.

## Features

- Visual comparison of phylogenetic trees
- Easy integration with existing Python workflows

## Installation

### Installation Steps

1. Install ADViewpy package from PyPI:

   ```
   pip install adviewpy
   ```

2. Import the ADViewpy package and initialize it in your Python script:

   ```
   import ADViewpy
   
   # Initialize the ADViewpy tool
   ADViewpy.init()
   ```

## Usage

#### `init(treefile)`

Initializes ADViewpy and creates the reference tree.

**Parameter:**

- `treefile` (str): Path to the reference tree file in Newick format.

**Return:**

- `ADViewpy` instance



#### `add_tree_collection(treefile,namefile)`

Add tree collection to ADViewpy

**Parameters:**

- `treefile` (str): Path to the tree collection file (each line contains a tree in newick format).
- `nanefile` (str): Path to the tree collection's name file, a name of tree each line corresponding to tree collection file.



#### `set_outgroup(outgroup_taxon[])`

Specify outgroup of phylogenetic tree, ADViewpy will re-root the trees according to outgroup provdied.

**Parameter:**

- `outgroup_taxon[]` (list): List of taxa names in outgroup

  

#### `reference_tree(exact_match_range[],support_value_range[])`

Visualization of reference tree, user can find branch with respect to attribute range. 

**Parameters:**

- `exact_match_range[]` (list): A list to represent node exact match percentage range, for example [80, 100], which represents to find nodes with exact match percentage between 80% and 100%.

- `support_value_range` (list): A list to represent node support value range, for example [50, 80], which represents to find nodes with support value between 50 and 80.

  

#### `AD(view,scale,context_level,max_ad,ad_interval,tree_id,tree_name,sort,escape_taxa_as_context_block,show_block_proportional,subtree_independent,show_tree_name,filter,differentiate_inexact_match)`

Visualization of individual and cluster aggregated dendrograms, to show topological relationships between the selected focal subtree.

**Common parameters:**

- `view` (str): `AD Individual` to show individual aggregated dendrograms, `AD Cluster` to show cluster aggregated dendrograms.
- `scale` (foat): Set the scale of aggregated dendrograms, default scale is 0.1.

**Parameters in individual view:**

- `context_level` (int): Control how many context blocks to show.
- `max_ad` (int): Set the maximum number of aggregated dendrograms to display.
- `ad_interval`(list): List to represent tree id range, for example [1, 20], which represents to show aggregated dendrograms with tree id between 1 and 20.
- `tree_id` (list): Each integer in list represents a tree ID. Only trees matching these IDs will be displayed.
- `tree_name`(list): Each str in list represents a tree name. Only trees matching these names will be displayed.
- `sort`(str): Specify how to sort the displayed aggregated dendrograms, can choose from `id`, `name`, and `rf distance`.
- `escape_taxa_as_context_block`(Boolean): Determine whether escape taxa should be considered as a context block.
- `show_block_proportional` (Boolean): Determine whether to use color coverage to indicate the percentage of highlighted taxa within subtree
- `subtree_independent` (Boolean):  Determine whether to merge overlapping subtrees into a single block.
- `show_tree_name` (Boolean): Determine whether to show tree id and tree name on top of each aggregated dendrogram.
- `filter` (str):  Specify whether the conditions are `include` or `exclude` for the following parameters: ad_interval, tree_id,tree_name. 

**Parameter in cluster view:**

- `differentiate_inexact_match` (Boolean): Choose whether to further distinguish exact and inexact matches in cluster

  

#### `pairwise_comparison(compare_tree)`

Pairwise compare two trees in detail. Allow to pairwise compare one tree in the collection with the reference tree or compare two trees in tree collection.

To compare a tree in the collection with the reference tree, user can choose a tree from AD()'s output view with mouse click.

**Parameter:**

`compare_tree` (int/str/list): Specify the trees for pairwise comparison. A tree can be specified using IDs or names. Providing an ID or name selects a single tree from the tree collection to compare with a reference tree. If a list is provided containing two trees, it indicates comparison between trees within the collection.



#### `tree_distribution()`

Showing distributions of trees in terms of their agreement and conflicts of leaf memberships for user-selected subtrees.



#### `tree_distance()`

Display the distances between each trees using a scatter plot.

**Return:**

- `Plotly` scatter showing distances between trees computed by the t-SNE dimensionality reduction technique.



#### `export_image(view,filename,*args)`

Exports generated view as an image

**Parameter:**

`view` (str): Selects which view to export. Choices include `Reference Tree`,`AD Individual`,`AD Cluster`,`Tree Distribution`,`Tree Distance`,`Pairwise Comparison`. 

