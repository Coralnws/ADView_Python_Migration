import copy

import dendropy
from Utils import *
from ipycanvas import Canvas
from IPython.display import display
import os
import myCanvas as myCanvas
import rtCanvas as rtCanvas
import tcCanvas as tcCanvas
import pairwiseCanvas as pairwiseCanvas
import importlib
importlib.reload(rtCanvas)
importlib.reload(tcCanvas)
import math
import numpy as np
import plotly.express as px
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE


class ADpy:
    # Testing
    output = []
    output2 = ""

    def __init__(self,treefile = None,type="newick"):
        # Reference Tree related
        self.rt = None  # Reference Tree
        self.rt_taxa = []
        self.rt_canvas = None  # Reference Tree's Canvas
        self.rt_view_support = False  # Whether to show internal node's support value
        self.outgroup = []
        self.rt_alter = False
        self.rt_file = treefile
        self.tree_schema = type
        self.tree_distance_matrix = None
        self.default_rt = None

        # AD/Tree collection related
        self.tc = []  # Tree Collection
        self.tc_taxa = []
        self.subtree_list = []
        self.ad_individual_canvas = None # Individual ADs' Canvas
        self.ad_cluster_canvas = None # Cluster ADs' Canvas
        self.individual_canvas_parameters = {}   # { 'parameter_name' : value }
        self.cluster_canvas_parameters = {}
        self.ad_parameter_alter = True # Record whether ad_canvas's parameter change

        # Paiwise canvas related
        self.pairwise_canvas = None


        self.default_subtree_attribute()
        # Reaf tree file and construct reference tree
        self.read_rt(treefile = treefile, type = type)


    #####  Public Functions ####
    def reference_tree(self,view_support=False):
        # Calculate height of canvas
        height = get_leaf_node_amount(self.rt) * RT_Y_INTERVAL + RT_Y_INTERVAL

        if not self.rt_canvas and not self.pairwise_canvas:
            self.default_subtree_attribute()
            self.default_rt = None
            self.rt_view_support = view_support
            self.create_pairwise_rt_canvas(alter_type=BOTH)

        if self.check_parameter_alter(view_support):
            self.create_pairwise_rt_canvas(alter_type=RT)

        return self.rt_canvas


    def set_outgroup(self,outgroup_taxon):
        # Reconstruct Tree
        self.read_rt(treefile=self.rt_file, type=self.tree_schema)

        # ["Uronema sp", "Monomastix opisthostigma", "Pyramimonas parkeae", "Nephroselmis pyriformis"]
        mrca = self.rt.mrca(taxon_labels=outgroup_taxon)
        self.rt.reroot_at_edge(mrca.edge)
        self.outgroup = outgroup_taxon

        if len(self.tc) > 1:
            for tree in self.tc:
                mrca = tree.mrca(taxon_labels=self.outgroup)
                tree.reroot_at_edge(mrca.edge)

            self.corresponding_branches()

            if self.tree_distance_matrix:
                self.generate_tree_distance_matrix()

        self.rt_alter = True

    def add_tree_collection(self,treefile=None,type="newick",namefile=None):
        self.tc = dendropy.TreeList.get_from_path(treefile, schema=type,taxon_namespace = self.rt.taxon_namespace)

        # Read tree collection name
        if namefile:
            file = open(namefile, "r")

        for index,tree in enumerate(self.tc):
            # Set trees' index, taxa_list, missing taxa, name and outgroup
            tree.id = index + 1
            tree.taxa_list = [leaf.taxon.label for leaf in tree.leaf_nodes()]
            tree.missing = set(self.rt.taxa_list) - set(tree.taxa_list)
            tree.pairwise_canvas = None

            if namefile:
                tree_name = file.readline().strip()
                tree.name = tree_name

            if len(self.outgroup) > 1:
                mrca = tree.mrca(taxon_labels=self.outgroup)
                tree.reroot_at_edge(mrca.edge)

            tree.rf_distance = dendropy.calculate.treecompare.symmetric_difference(self.rt, tree)

        self.generate_tree_distance_matrix()
        # self.rt_alter = True
        self.corresponding_branches()
        
    def tree_collection(self,sort_by=ID):
        tc_list = self.sort_tc(self.tc,sorted)
         
        for tree in tc_list:
            print("Tree "  + str(tree.id))
            print(" " * 2 + "Name: " + tree.name)
            print(" " * 2 + "Distance: " + str(tree.rf_distance))
            print("\n")

    def sort_tc(self,tc_list,sort_by):
        if sort_by == ID:
            return sorted(tc_list, key=lambda x: (x.id),reverse=False)
        elif sort_by == RF_DISTANCE:
            return sorted(tc_list, key=lambda x: (x.rf_distance, x.name),reverse=False)
        elif sort_by == NAME:
            return sorted(tc_list, key=lambda x: (x.name),reverse=False)
        
    # Same effect as click on the common ancestor of these nodes
    def select_subtree(self,nodes=None):
        nodes_list = []
        subtree_root = None

        if not self.rt_canvas:
            self.rt_not_exist_error()
            return

        # Approximate taxon name
        for node in nodes:
            for leaf_node in self.rt.leaf_node_iter():
                if node.lower() in leaf_node.taxon.label.lower():
                    if len(nodes) < 2:
                        self.rt_canvas.draw_subtree_block(leaf_node.parent_node)
                        return

                    nodes_list.append(leaf_node.taxon.label)

        subtree_root = self.rt.mrca(taxon_labels=nodes_list)

        if self.rt_canvas:
            self.rt_canvas.draw_subtree_block(subtree_root)

    # Same effect as click on the leaf_node
    def select_taxa(self, node=None):
        if not self.rt_canvas:
            self.rt_not_exist_error()
            return

        # Approximate taxon name
        for leaf_node in self.rt.leaf_node_iter():
            if node.lower() in leaf_node.taxon.label.lower():
                if self.rt_canvas:
                    self.rt_canvas.draw_subtree_block(leaf_node)

        # Exact taxon name
        # node = self.rt.find_node_with_taxon_label(node)
        # self.rt_canvas.draw_subtree_block(node)

    # Show AD, default: show individual AD
    def AD(self,view=AD_INDIVIDUAL,scale=1.0,max_ad=None,context_level=2,ad_interval=[],tree_id=None,
           tree_name=None,filter=INCLUDE,sort=RF_DISTANCE,espace_taxa_as_context_block=True,show_block_proportional=True,
           subtree_independent=False,parameter_from_individual_ad=True,differentiate_inexact_match=True,
           show_tree_name=False):

        first_ad = None
        last_ad = None
        if len(ad_interval) > 0:
            first_ad = ad_interval[0]
            last_ad = ad_interval[1]

        if tree_id and type(tree_id) is not list:
            new_list = []
            new_list.append(tree_id)
            tree_id = new_list
        if tree_name and type(tree_name) is not list:
            tree_name = [tree_name]


        # Check if condition is logical
        # 1. First_ad < last_ad
        if first_ad and last_ad:
            if last_ad < first_ad:
                self.parameter_error("ad_interval")
                return
        # 2.  First_ad > 0
        if first_ad and first_ad <= 0:
            self.parameter_error("ad_interval")
            return
        # 3. Last_ad > 0
        if last_ad and last_ad <= 0:
            self.parameter_error("ad_interval")
            return
        # 4. tree_id > 0
        if tree_id:
            for check_id in tree_id:
                if check_id <= 0:
                    self.parameter_error("tree_id")
            return
        # 5. context_level >= 1
        if context_level < 1:
            self.parameter_error("context_level")
            return
        # Parameter: width, height, context levels, show labels, show colors
        if scale != 1.0:
            ad_per_row = (CANVAS_MAX_WIDTH - (3 * DEFAULT_PADDING_BETWEEN_AD)) // (DEFAULT_AD_WIDTH * scale)
            if scale < 0.5:
                ad_per_row *= 3

        else:
            ad_per_row = DEFAULT_AD_PER_ROW

        if view == AD_INDIVIDUAL:
            # if not self.ad_individual_canvas or self.parameter_modified():
            if not self.ad_individual_canvas or self.ad_parameter_alter:
                if max_ad:
                    ad_row = math.ceil(max_ad / ad_per_row)
                else:
                    ad_row = math.ceil(len(self.tc) / ad_per_row)

                canvas_height = (ad_row * DEFAULT_AD_HEIGHT * scale) + (2 * DEFAULT_PADDING_BETWEEN_AD) +  ((ad_row -
                                                                                                             1) * DEFAULT_PADDING_BETWEEN_AD)
                self.ad_individual_canvas = tcCanvas.tcCanvas(layer = ad_row,adPy=self,view=view,
                                                              width=CANVAS_MAX_WIDTH,
                                                              height=canvas_height, scale=scale,max_ad=max_ad,
                                                              ad_per_row=ad_per_row, context_level=context_level,
                                                              first_ad=first_ad, last_ad=last_ad,tree_id=tree_id,
                                                              tree_name=tree_name, sort_by=sort,
                                                              espace_taxa_as_context_block=espace_taxa_as_context_block,
                                                              show_block_proportional=show_block_proportional,
                                                              subtree_independent=subtree_independent,show_tree_name=show_tree_name)
                display(self.ad_individual_canvas)
                # for layer in range(0, self.ad_individual_canvas.layer):
                for index, ad_tree in enumerate(self.ad_individual_canvas.ad_list):
                    self.ad_individual_canvas.draw_ad_tree(ad_tree)
                    self.ad_individual_canvas[ad_tree.located_layer].flush()

                # return self.ad_individual_canvas

        # Parameter: differentiate inexact match, differentiate sister-group relationships
        elif view == AD_CLUSTER:
            if not self.ad_cluster_canvas or self.ad_parameter_alter:
                self.ad_cluster_canvas = tcCanvas.tcCanvas(adPy=self,view=view,width=CANVAS_MAX_WIDTH,height=150,
                                                           scale=scale,context_level=context_level,
                                                           ad_per_row=ad_per_row,
                                                           show_block_proportional=show_block_proportional,
                                                           subtree_independent=subtree_independent,
                                                           parameter_from_individual_ad=parameter_from_individual_ad,
                                                           espace_taxa_as_context_block=espace_taxa_as_context_block,
                                                           differentiate_inexact_match=differentiate_inexact_match)

                return self.ad_cluster_canvas

    def tree_distance(self):
        # self.tree_distance_matrix = np.random.rand(69, 69)
        # self.tree_distance_matrix = (self.tree_distance_matrix + self.tree_distance_matrix.T) / 2
        # np.fill_diagonal(self.tree_distance_matrix, 0)

        # mds = MDS(n_components=2, dissimilarity='precomputed',random_state=42)
        # self.tree_point_coordinates = mds.fit_transform(self.tree_distance_matrix)
        # self.tree_point_coordinates -= self.tree_point_coordinates.min(axis=0)

        tsne = TSNE(n_components=2)
        self.tree_point_coordinates = tsne.fit_transform(self.tree_distance_matrix)
        # self.tree_point_coordinates -= self.tree_point_coordinates.min(axis=0)
        # self.tree_point_coordinates = self.tree_point_coordinates.astype(np.float_)

        self.x_coor = []
        self.y_coor = []
        for coordinate in self.tree_point_coordinates:
            self.x_coor.append(coordinate[0])
            self.y_coor.append(coordinate[1])

        id_name = [f"0 : {self.rt.name}"]
        for tc_tree in self.tc:
            id_name.append(f"{tc_tree.id} : {tc_tree.name}")

        color = ['Reference Tree']
        for i in range(len(self.tree_point_coordinates)-1):
            color.append('Tree Collection')

        fig = px.scatter(x=self.x_coor, y=self.y_coor, color=color,
                         title='Tree Distance')

        # Customize hover template
        fig.update_traces(hovertemplate='%{text}',text=id_name, hoverinfo='text')
        # Remove legend labels
        fig.update_layout(dragmode=False,showlegend=False,width=600, height=600,plot_bgcolor=TREE_NAME_BG,xaxis=dict(
            color=BLANK),yaxis=dict(color=BLANK))

        # fig = px.scatter(x=self.tree_point_coordinates[:, 0], y=self.tree_point_coordinates[:, 1], color=color,
        #                  title='Tree Distance')
        # fig.update_layout(showlegend=False)
        # fig.update_traces(text=id_name, hoverinfo='text', hovertemplate='%{text}')


        fig.show()

    def remove_taxa(self,tree,taxa_list=None):
        remove_taxa_list = []
        for i in range(10,len(self.rt.taxa_list)-1):
            taxa_list = (self.rt.taxa_list[i])

        # Iterate over the tree and remove the specified taxa
        for taxon_label in remove_taxa_list:
            self.rt.prune_taxa_with_labels(taxon_label)

        # Save the modified tree to a new file
        tree.write(path="modified_tree_file.tre", schema="newick")

    def pairwise_comparison(self,compare_tree=None):
        # adPy,layer, width, height,rt_canvas,tc_tree
        # Testing
        if not self.rt_canvas and not self.pairwise_canvas:
            self.create_pairwise_rt_canvas()

        if self.tc:
            compare_tree = self.tc[1]
        # 36

        self.pairwise_canvas.compare_tc_tree(compare_tree)

        return self.pairwise_canvas

    def tree_distribution_view(self):
        print(self.subtree_list)
        self.tree_distribution_view = tcCanvas.tcCanvas(adPy=self, view=GET_TREE_DISTRIBUTION)

        for subtree in self.subtree_list:
            print(f"Subtree {subtree.label}")
            index = 1
            print(subtree.topology_list)
            for topology,tree_list in subtree.topology_list.items():
                print(f"Topology {index} : {len(tree_list)}")
                print(topology)
                index += 1

        #### Internal Functions ####
    def read_rt(self,treefile,type):
        self.rt = dendropy.Tree.get(path=treefile, schema=type)

        # Filename as tree name
        rt_label = os.path.basename(treefile)
        self.rt.id = 0
        self.rt.name = rt_label
        self.rt.level = 0
        self.rt.pairwise_canvas = None

        # Get tree's taxa/leaf node
        self.rt.taxa_list = [leaf.taxon.label for leaf in self.rt.leaf_nodes()]

    def generate_tree_distance_matrix(self):
        dimension = len(self.tc) + 1
        self.tree_distance_matrix = np.zeros((dimension, dimension))

        for index,tree in enumerate(self.tc):
            self.tree_distance_matrix[0,index + 1] = tree.rf_distance

        for i in range(0,len(self.tc)):
            tree_compare = self.tc[i]
            # print("Compare with tree " + str(i))
            for j in range(i+1,len(self.tc)):
                # print(" " * 5 + "Tree " + str(j))
                distance = dendropy.calculate.treecompare.symmetric_difference(tree_compare,self.tc[j])
                self.tree_distance_matrix[i + 1, j + 1] = distance / 1.1
                # print(f"[{i+1},{j+1}] = " + str(distance))

        self.tree_distance_matrix = np.triu(self.tree_distance_matrix) + np.triu(self.tree_distance_matrix, k=1).T

    def corresponding_branches(self):
        if self.tc == None:
            print("Tree Collection Not Exist")
            return

        for node in self.rt.postorder_node_iter():
            node_taxa = [leaf.taxon.label for leaf in node.leaf_nodes()]
            node.corr = [node]
            node.exact_match = 0
            node.exact_match_tree = []

            for tree in self.tc:
                target_set = set(node_taxa) - tree.missing
                self.rt.missing = tree.missing - set(node_taxa)

                node.corr_similarity = 0
                node.corr.append(0)

                for tc_node in tree.postorder_node_iter():
                    if not hasattr(tc_node,'corr'):
                        tc_node.corr = None
                        tc_node.corr_similarity = 0

                    similarity = self.get_similarity(target_set=target_set,node=node,tc_node=tc_node)
                    if similarity > node.corr_similarity:
                        node.corr[tree.id] = tc_node
                        node.corr_similarity = similarity

                        if similarity > tc_node.corr_similarity:
                            tc_node.corr = node
                            tc_node.corr_similarity = similarity

                        # print(similarity)
                        if similarity == 1.0:
                            node.exact_match += 1
                            node.exact_match_tree.append(tree)


    def get_similarity(self,target_set,node,tc_node):
        if tc_node.is_leaf():
            tc_node.card = 0 if tc_node.taxon.label in self.rt.missing else 1
            tc_node.intersect = 1 if tc_node.taxon.label in target_set else 0
        else:
            child_nodes = tc_node.child_nodes()
            tc_node.card = 0
            tc_node.intersect = 0
            for child in child_nodes:
                tc_node.card += child.card
                tc_node.intersect += child.intersect

        union = len(target_set) + tc_node.card - tc_node.intersect
        return tc_node.intersect/union

    def is_descendant(self,parent_node,child_node):
        for node in parent_node.child_nodes():
            if self.is_descendant(node,child_node):
                return True

        if parent_node == child_node:
            return True

        return False

    def rt_not_exist_error(self):
        print("<Error> : Reference Tree not exist.")
        print("Please ensure that you have called the reference_tree() function before calling this function.")

    def parameter_error(self,parameter):
        print(f"<Error> : {parameter} given was incorrect.")
        print("Please ensure that the information provided is logical.")

    # Manage Subtree
    def default_subtree_attribute(self):
        self.subtree_list = []
        self.subtree_label_used = [1, 1, 1, 1, 1]
        self.subtree_color_used = [1, 1, 1, 1, 1]

    def create_pairwise_rt_canvas(self,alter_type=BOTH):
        height = get_leaf_node_amount(self.rt) * RT_Y_INTERVAL + RT_Y_INTERVAL * 2

        if self.rt_alter:
            self.default_rt = None
            self.rt_alter = False

        # Check paramater alter
        if alter_type == BOTH:
            # self.default_rt_value()
            self.rt_canvas = rtCanvas.rtCanvas(self, width=CANVAS_MAX_WIDTH, height=height ,
                                               view_support=self.rt_view_support, default_rt=self.default_rt)
            self.pairwise_canvas = pairwiseCanvas.pairwiseCanvas(self, width=CANVAS_MAX_WIDTH, height=height)
        elif alter_type == RT:
            self.rt_canvas = rtCanvas.rtCanvas(self, width=CANVAS_MAX_WIDTH, height=height,
                                               view_support=self.view_support, default_rt=self.default_rt)
        elif alter_type == PAIRWISE:
            self.pairwise_canvas = pairwiseCanvas.pairwiseCanvas(self, width=CANVAS_MAX_WIDTH, height=height)


    def check_parameter_alter(self,view_support):
        # Check whether parameter has alter and set parameter as self.attribute
        if view_support != self.rt_view_support:
            self.rt_view_support = view_support
            return True

    def select_subtree_from_tree(self,node_selected):
        if not hasattr(node_selected, 'selected') or node_selected.selected == False:
            # Ignore if 5 subtree had selected
            if len(self.subtree_list) >= 5:
                return None

            # Choose block color
            color_index = self.subtree_color_used.index(1)
            self.subtree_color_used[color_index] = 0
            color = SUBTREE_COLOR_LIST[color_index]

            # Record new subtree - choose subtree label
            subtree_label_index = self.subtree_label_used.index(1)
            self.subtree_label_used[subtree_label_index] = 0
            label = SUBTREE_LABEL_LIST[subtree_label_index]

            # Create new subtree (class from myTree.py)
            new_subtree = Subtree(label=label, belong_tree=self.rt, root=node_selected, color=color)

            # Mark the layer as occupied and record corresponding subtree's label
            node_selected.subtree = new_subtree
            node_selected.selected = True

            self.rt_canvas.draw_subtree_block(node_selected,new_subtree)
            self.subtree_list.append(new_subtree)
            self.pairwise_canvas.draw_subtree_block(node_selected,select_tree=RT,new_subtree=new_subtree)

            if self.pairwise_canvas.tc_tree:
                self.pairwise_canvas.draw_tc_subtree_block(new_subtree)
                self.pairwise_canvas.draw_escape_taxa()


            # return new_subtree

        else:
            self.rt_canvas.remove_subtree_block(node_selected.subtree)
            self.pairwise_canvas.remove_subtree_block(node_selected.subtree)
            self.subtree_list.remove(node_selected.subtree)

            # Release Color

            color_index = SUBTREE_COLOR_LIST.index(node_selected.subtree.color)
            self.subtree_color_used[color_index] = 1
            # Release Subtree Label
            label_index = SUBTREE_LABEL_LIST.index(node_selected.subtree.label)
            self.subtree_label_used[label_index] = 1

            # Remove subtree from node and mark node as not-selected
            node_selected.subtree = None
            node_selected.selected = False

            # if len(self.subtree_list) == 0:
            #     self.rt_canvas.reset_subtree_canvas()
            #     self.pairwise_canvas.reset_subtree_canvas()


class Subtree:
    # rt = None  # Belong to which reference tree
    # label = None  # A,B,C,D,E
    # rtCanvas_index = None
    # root = None # Subtree root
    # color = None
    # block = None
    # label_width = None
    # leaf_set = []
    # ad_block_list = []   # index = tc_index, value = [AD_Block list] : to record corresponding block in tree from tree collection

    def __init__(self, label, belong_tree, root, color):
        self.label = label
        self.belong_tree = belong_tree
        self.root = root
        self.color = color
        self.corresponding_tc_subtree = None

        self.label_width = None
        self.leaf_set = set()

        self.topology_list = {}
    def subtree_to_string(self,node):
        if len(node.child_nodes()) == 0:
            return node.taxon.label
        else:
            child_strings = [self.subtree_to_string(child) for child in node.child_nodes()]
            return "(" + ",".join(child_strings) + ")"

    def check_and_set_topology(self,root,tc_tree_id,type=NODE_DIFFERENCE):
        if type == TOPOLOGY_DIFFERENCE:
            subtree_string = self.subtree_to_string(root)
        elif type == NODE_DIFFERENCE:
            children = [child.taxon.label for child in root.leaf_nodes()]
            children = sorted(children)
            subtree_string = ",".join(children)

        if subtree_string not in self.topology_list:
            self.topology_list[f"{subtree_string}"] = []

        self.topology_list[f"{subtree_string}"].append(tc_tree_id)

    def set_rtLayer(self, rtCanvas_index):
        self.rtCanvas_index = rtCanvas_index

    def get_leaf_nodes(self):
        self.leaf_set = set()
        for leaf_node in self.root.leaf_nodes():
            self.leaf_set.add(leaf_node.taxon.label)


    def set_block(self,block):
        self.block = block
        self.block_size = block.get_size()

    def set_pairwise_block(self,pairwise_block):
        self.pairwise_block = pairwise_block
        self.block_size = pairwise_block.get_size()





# AD
class AD_Tree:
    # Internal node is AD_Block
    # Leaf node is AD_Block
    # Tree_width = self.width - padding
    padding = DEFAULT_PADDING_BETWEEN_BLOCK  # Padding between block

    def __init__(self, id, tc_tree, tc_canvas):
        self.id = id  # id for AD as well as tree id
        self.root = None  # AD_Node
        self.tc_tree = tc_tree
        self.block_list = []
        self.tc_canvas = tc_canvas
        self.nested_cnt = 0
        self.is_nested = False  # Is/Not nested tree

        self.width = DEFAULT_AD_WIDTH  # Default width
        self.height = DEFAULT_AD_HEIGHT  # Default height
        self.topL = None
        self.botR = None
        self.index_in_row = 0
        self.set_default()
        self.space_required = 0

        self.located_layer = 0
        self.tree_name_block = None
        self.missing_taxa_list = None
        self.missing_taxa_block = None
        self.missing_taxa_segment = []

    def set_default(self):
        self.x = 8
        self.y = 8
        self.padding = 5  # Padding between block


    def set_position_size(self, x, y):
        self.topL = myCanvas.Point(x, y)
        self.width = self.tc_canvas.scale * DEFAULT_AD_WIDTH
        self.height = self.tc_canvas.scale * DEFAULT_AD_HEIGHT
        self.botR = myCanvas.Point(x + self.width, y + self.height)

    def set_nested_tree_size(self, nested_block):
        self.topL = nested_block.topL
        self.botR = nested_block.botR
        self.width = self.botR.x - self.topL.x

    def plot_tree(self, node=None, level=0):
        if node is None:
            node = self.root

        if node.type == ROOT or node.type == INTERNAL:
            internaL_node = node.node_or_block
            print(' ' * (level * 4) + node.type, end=" ")
            print(internaL_node)

        elif node.type == LEAF:
            block = node.node_or_block

            if block.type == SUBTREE_BLOCK:
                print(' ' * (level * 4) + block.belong_subtree.label + ':' + str(block.subtree_taxa_count))
            elif block.type == INDIVIDUAL_BLOCK:
                if block.belong_subtree:
                    print(' ' * (level * 4) + block.belong_subtree.label + ': INDIVIDUAL LEAF')
                else:
                    print(' ' * (level * 4) + 'INDV BLANK BLOCK')

        for child in node.children:
            self.plot_tree(child, level + 1)

    def generate_block_list(self):
        self.block_list = []

        for node in self.traverse_postorder():
            if node.type == LEAF:
                self.block_list.append(node.node_or_block)

                if node.nested_tree:
                    node.nested_tree.generate_block_list()

    def insert_node(self, parent, child):
        parent.children.append(child)

        # if child.type == LEAF:
        #     self.block_list.append(child.node_or_block)

    def traverse_postorder(self, node=None, node_list=None):
        if node == None:
            node = self.root
            node_list = []

        node.children = sorted(node.children, key=lambda x: x.child_index, reverse=False)
        for child in node.children:
            node_list = self.traverse_postorder(child, node_list)

        node_list.append(node)
        return node_list

    def individual_subtree_block_list(self):
        indv_block_list = []
        subtree_block_list = []
        unnested_block_list = []
        unnested_block_taxa = 0

        for block in self.block_list:
            if block.type == INDIVIDUAL_BLOCK or block.type == INDIVIDUAL_LEAF_BLOCK:
                indv_block_list.append(block)
            elif block.type == SUBTREE_BLOCK:
                if block.nested_tree:
                    subtree_block_list.append(block)
                else:
                    unnested_block_list.append(block)
                    unnested_block_taxa += block.subtree_taxa_count

        for block in unnested_block_list:
            subtree_block_list.append(block)

        return indv_block_list, subtree_block_list, unnested_block_taxa

    def get_all_subtree_block_list(self,subtree_block_list):
        for block in self.block_list:
            if block.type == SUBTREE_BLOCK:
                subtree_block_list.append(block)
                if block.nested_tree:
                    block.nested_tree.get_all_subtree_block_list(subtree_block_list)

    def ad_taxa_total(self):
        indv_block_list, subtree_block_list = self.individual_subtree_block_list()
        taxa_cnt = 0

        for block in subtree_block_list:
            taxa_cnt += block.subtree_taxa_count
            if block.nested_tree:
                taxa_cnt += block.nested_tree.ad_taxa_total()

        return taxa_cnt

    def ad_to_string(self, canvas,node=None, differentiate_inexact_match=True):

        if node == None:
            node = self.root
            newick_str = ""

        if node.type == LEAF:
            return self.get_node_type(canvas,node,differentiate_inexact_match=differentiate_inexact_match)
        else:
            newick_str = ""
            sorted_children = sorted(node.children, key=lambda x: x.type_prior, reverse=False)
            for child in sorted_children:
                nested_tree_string = ""
                if child.nested_tree:
                    nested_tree_string = child.nested_tree.ad_to_string(canvas=canvas,
                                                                        differentiate_inexact_match=differentiate_inexact_match)
                    newick_str += self.ad_to_string(canvas=canvas,node=child,
                                                    differentiate_inexact_match=differentiate_inexact_match) + '[' + nested_tree_string + '],'
                else:
                    newick_str += self.ad_to_string(canvas=canvas,node=child,
                                                    differentiate_inexact_match=differentiate_inexact_match) + ','

            if node.type == ROOT:
                return '(' + newick_str[:-1] + ')' + self.get_node_type(canvas,node,differentiate_inexact_match=differentiate_inexact_match) + ';'
            else:
                return '(' + newick_str[:-1] + ')' + self.get_node_type(canvas,node,differentiate_inexact_match=differentiate_inexact_match)

    def get_node_type(self, canvas,node,differentiate_inexact_match=True):
        if node.type == LEAF:
            block = node.node_or_block
            if block.type == SUBTREE_BLOCK:
                subtree_label = ""
                if block.is_subtree_duplicate:
                    for index,subtree in enumerate(block.belong_subtree):
                        str = ""
                        if differentiate_inexact_match:
                            if block.exact_match[index]:
                                str = "[exact]"
                            else:
                                str = "[inexact]"
                        elif not block.exact_match[index]:
                            canvas.inexact_block.append(subtree.label)

                        subtree_label += str + subtree.label + '&'
                    return subtree_label[:-1]
                else:
                    str = ""
                    if differentiate_inexact_match:
                        if block.exact_match:
                            str = "[exact]"
                        else:
                            str = "[inexact]"
                    elif not block.exact_match:
                        canvas.inexact_block.append(block.belong_subtree.label)
                        
                    return str + block.belong_subtree.label
            elif block.type == INDIVIDUAL_LEAF_BLOCK:
                return block.belong_subtree.label + "_INDEPENDENT"
            elif block.type == INDIVIDUAL_BLOCK:
                return block.type

        elif node.type == INTERNAL:
            return node.type
        else:
            return ""

    def check_in_range(self,x,y):
        return x >= self.topL.x and x <= self.botR.x and y >= self.topL.y and y <= self.botR.y





class AD_Node:
    def __init__(self, node_or_block, x_level, type, child_index, type_prior):
        self.node_or_block = node_or_block  # Dendropy node if internal node, AD_Block if leaf
        self.children = []
        self.x_level = x_level
        self.type = type  # LEAF or INTERNAL
        self.pos = None
        self.branch_head = None
        self.branch_tail = None
        self.nested_tree = None  # Root of nested tree
        self.is_elided = False
        self.child_index = child_index
        self.type_prior = type_prior

    def insert_child(self, ad_node):
        self.children.append(ad_node)

    def construct_branch(self, x, y, scale):
        self.branch_head = myCanvas.Point(x - (AD_BRANCH_LENGTH * scale), y)
        self.branch_tail = myCanvas.Point(x, y)

    def construct_elided_branch(self, x, y, scale):
        self.branch_head = myCanvas.Point(x - (ELIDED_BRANCH_LENGTH * scale), y)
        self.branch_tail = myCanvas.Point(x, y)

    # Standardize the starting x-coordinate of all children's branches
    def unify_children_branches(self, x):
        for child in self.children:
            child.branch_head.x = x


class AD_Topology:
    def __init__(self, ad_string, sample_ad_tree):
        self.ad_string = ad_string
        self.sample_ad_tree = sample_ad_tree
        self.ad_tree_list = []
        self.tree_count = 1
        self.block_list = []
        self.generate_block_list()

    def add_ad(self, ad_tree):
        self.ad_tree_list.append(ad_tree)
        self.tree_count += 1

    def generate_block_list(self):
        self.sample_ad_tree.get_all_subtree_block_list(self.block_list)

    def set_block_inexact_match(self,subtree_label):
        for block in self.block_list:
            if block.is_subtree_duplicate:
                for index,subtree in enumerate(block.belong_subtree):
                    if subtree_label == subtree.label:
                        block.exact_match[index] = False
            else:
                if subtree_label == block.belong_subtree.label:
                    block.exact_match = False

def print_tree(tree):
    print(tree.as_ascii_plot())

def get_node_list(tree):
    return tree.leaf_node_iter()

def get_leaf_node_amount(tree):
    return sum(1 for node in tree.leaf_node_iter())

