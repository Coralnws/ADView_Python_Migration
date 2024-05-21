import copy
import re
import dendropy
from ipycanvas import Canvas,hold_canvas
from IPython.display import display
import os
import plotly.io as pio
from myUtils import *
from myCanvas import *
from rtCanvas import rtCanvas
from tcCanvas import tcCanvas
from treeDistributionView import TreeDistributionView
from pairwiseCanvas import pairwiseCanvas
import math
import numpy as np
import plotly.express as px
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

def init(treefile = None,type="newick"):
    new_adpy = ADViewpy(treefile, type)
    return new_adpy

class ADViewpy:
    output = []
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
        self.tree_distribution_tmp = None
        self.tc_canvas_tmp = None
        self.ad_parameter_alter = True # Record whether ad_canvas's parameter change

        self.tree_distribution_view = None
        self.tree_distance_fig = None

        self.ad_cluster_canvas_export = None
        self.ad_individual_canvas_export = None

        self.rt_exact_match_range = None
        self.rt_support_value_range = None

        self.pairwise_tree = None


        # Paiwise canvas related
        self.pairwise_canvas = None
        self.default_subtree_attribute()
        # Reaf tree file and construct reference tree
        self.read_rt(treefile = treefile, type = type)


    # ==========================  Public Functions ================================== #

    # reference_tree() : Function to show reference tree visualization
    # Parameter: view_support (Whether to write support value of internal node on tree)
    # Draw reference tree on canvas
    # User can select subtree through mouse click
    def reference_tree(self,view_support=False,show=True,exact_match_range=None,support_value_range=None):
        # Calculate height of canvas
        height = get_leaf_node_amount(self.rt) * RT_Y_INTERVAL + RT_Y_INTERVAL
        if exact_match_range and not self.check_attribute_range(exact_match_range):
            self.parameter_error("exact_match_range")
            return

        if support_value_range and not self.check_attribute_range(support_value_range):
            self.parameter_error("support_value_range")
            return

        if not self.rt_canvas and not self.pairwise_canvas:
            self.default_subtree_attribute()
            self.default_rt = None
            self.rt_view_support = view_support
            self.rt_exact_match_range = exact_match_range
            self.rt_support_value_range = support_value_range
            self.create_pairwise_rt_canvas(alter_type=BOTH)


        else:
            check_parameter = self.check_parameter_alter(view_support, exact_match_range, support_value_range)

            if check_parameter == REDRAW:
                self.create_pairwise_rt_canvas(alter_type=RT)

            if check_parameter == FILTER_NODE:
                self.rt_canvas.draw_filter_node(exact_match_range, support_value_range)

        if show:
            return self.rt_canvas

    # set_outgroup() : Function to set trees' outgroup
    # Parameter : outgroup_taxon (python list)
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
        self.ad_parameter_alter = True

    # add_tree_collection() : Function to set tree collection
    # Parameter: treefile(filepath to read trees file), namefile(filepath to read tree name file)
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
            tree.missing_node_list = []
            for taxa in tree.missing:
                self.generate_missing_node(tree,taxa)

            tree.missing_node_list = sorted(tree.missing_node_list, key=lambda x: x.label, reverse=False)

            tree.pairwise_canvas = None

            if namefile:
                tree_name = file.readline().strip()
                tree.name = tree_name
            else:
                tree.name = ""

            if len(self.outgroup) > 1:
                mrca = tree.mrca(taxon_labels=self.outgroup)
                tree.reroot_at_edge(mrca.edge)

            tree.rf_distance = dendropy.calculate.treecompare.symmetric_difference(self.rt, tree)

        self.generate_tree_distance_matrix()
        # self.rt_alter = True
        self.corresponding_branches()

    # tree_collection() : Function to get tree info from tree collection
    def tree_collection(self,sort_by=ID):
        tc_list = self.sort_tc(self.tc,sorted)

        for tree in tc_list:
            print("Tree "  + str(tree.id))
            print(" " * 2 + "Name: " + tree.name)
            print(" " * 2 + "Distance: " + str(tree.rf_distance))
            print("\n")


    # select_subtree() : Function to select subtree
    # Parameter: nodes(python list, contain taxa name as string)
    # Same effect as click on the most recent common ancestor of these taxa
    def select_subtree(self,nodes=None):
        nodes_list = []
        subtree_root = None

        if not self.rt_canvas and not self.pairwise_canvas:
            self.create_pairwise_rt_canvas()

        # Approximate taxon name
        for node in nodes:
            for leaf_node in self.rt.leaf_node_iter():
                if node.lower() in leaf_node.taxon.label.lower():
                    if len(nodes) < 2:
                        self.rt_canvas.draw_subtree_block(leaf_node.parent_node)
                        return

                    nodes_list.append(leaf_node.taxon.label)

        subtree_root = self.rt.mrca(taxon_labels=nodes_list)

        self.select_subtree_from_tree(subtree_root)
        # if self.rt_canvas:
        #     self.rt_canvas.draw_subtree_block(subtree_root)

    # select_taxa() : Function to choose a taxa(one and only one) as subtree
    # Parameter: node(string - taxa name)
    # Same effect as click on the leaf_node
    def select_taxa(self, node=None):
        if not self.rt_canvas:
            self.not_exist_error(tree="Reference Tree",pre_function="reference_tree()")
            return

        # Approximate taxon name
        for leaf_node in self.rt.leaf_node_iter():
            if node.lower() in leaf_node.taxon.label.lower():
                self.select_subtree_from_tree(leaf_node)

        # Exact taxon name
        # node = self.rt.find_node_with_taxon_label(node)
        # self.rt_canvas.draw_subtree_block(node)

    # AD() : Function to get Aggregated Dendrogram, Individual AD or cluster AD
    def AD(self,view=AD_INDIVIDUAL,scale=1.0,max_ad=None,context_level=2,ad_interval=None,tree_id=None,
           tree_name=None,filter=INCLUDE,sort=RF_DISTANCE,escape_taxa_as_context_block=True,show_block_proportional=True,
           subtree_independent=False,parameter_from_individual_ad=True,differentiate_inexact_match=True,
           show_tree_name=False,export=False):

        if len(self.subtree_list) == 0:
            print("<Error> : No Subtree Chosen")
            return

        if not self.tc:
            self.not_exist_error(tree="Tree Collection",pre_function="add_tree_collection()")

        first_ad = None
        last_ad = None
        if ad_interval:
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
        if tree_id is not None:
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
            if not self.ad_individual_canvas or self.ad_parameter_alter or export:
                if max_ad:
                    ad_row = math.ceil(max_ad / ad_per_row)
                else:
                    ad_row = math.ceil(len(self.tc) / ad_per_row)

                canvas_height = (ad_row * DEFAULT_AD_HEIGHT * scale) + (2 * DEFAULT_PADDING_BETWEEN_AD) +  ((ad_row -
                                                                                                             1) * DEFAULT_PADDING_BETWEEN_AD)
                if export:
                    self.ad_individual_canvas_export = tcCanvas(layer=ad_row, adPy=self, view=view,
                                                                width=CANVAS_MAX_WIDTH, filter=filter,
                                                                height=canvas_height, scale=scale, max_ad=max_ad,
                                                                ad_per_row=ad_per_row, context_level=context_level,
                                                                first_ad=first_ad, last_ad=last_ad, tree_id=tree_id,
                                                                tree_name=tree_name, sort_by=sort,
                                                                escape_taxa_as_context_block=escape_taxa_as_context_block,
                                                                show_block_proportional=show_block_proportional,
                                                                subtree_independent=subtree_independent,
                                                                show_tree_name=show_tree_name)

                    self.ad_individual_canvas_export.TREE_NAME_LAYER = -3
                    for index, ad_tree in enumerate(self.ad_individual_canvas_export.ad_list):
                        with hold_canvas(self.ad_individual_canvas_export):
                            self.ad_individual_canvas_export.paste_ad_tree_canvas(ad_tree)
                        self.ad_individual_canvas_export[ad_tree.located_layer].flush()

                else:
                    self.ad_individual_canvas = tcCanvas(layer=ad_row, adPy=self, view=view,
                                                         width=CANVAS_MAX_WIDTH, filter=filter,
                                                         height=canvas_height, scale=scale, max_ad=max_ad,
                                                         ad_per_row=ad_per_row, context_level=context_level,
                                                         first_ad=first_ad, last_ad=last_ad, tree_id=tree_id,
                                                         tree_name=tree_name, sort_by=sort,
                                                         escape_taxa_as_context_block=escape_taxa_as_context_block,
                                                         show_block_proportional=show_block_proportional,
                                                         subtree_independent=subtree_independent,
                                                         show_tree_name=show_tree_name)
                    display(self.ad_individual_canvas)
                    for index, ad_tree in enumerate(self.ad_individual_canvas.ad_list):
                        with hold_canvas(self.ad_individual_canvas):
                            self.ad_individual_canvas.paste_ad_tree_canvas(ad_tree)
                        # with hold_canvas(self.ad_individual_canvas):
                        #     self.ad_individual_canvas.draw_ad_tree(ad_tree)
                        self.ad_individual_canvas[ad_tree.located_layer].flush()
            else:
                display(self.ad_individual_canvas)


        elif view == AD_CLUSTER:
            if not self.ad_cluster_canvas or self.ad_parameter_alter or export:
                if export:
                    self.ad_cluster_canvas_export = tcCanvas(adPy=self,view=view,width=CANVAS_MAX_WIDTH,height=150,
                                                           scale=scale,context_level=context_level,
                                                           ad_per_row=ad_per_row,
                                                           show_block_proportional=show_block_proportional,
                                                           subtree_independent=subtree_independent,
                                                           parameter_from_individual_ad=False,
                                                           escape_taxa_as_context_block=escape_taxa_as_context_block,
                                                           differentiate_inexact_match=differentiate_inexact_match)
                else:
                    self.ad_cluster_canvas = tcCanvas(adPy=self,view=view,width=CANVAS_MAX_WIDTH,height=150,
                                                               scale=scale,context_level=context_level,
                                                               ad_per_row=ad_per_row,
                                                               show_block_proportional=show_block_proportional,
                                                               subtree_independent=subtree_independent,
                                                               parameter_from_individual_ad=parameter_from_individual_ad,
                                                               escape_taxa_as_context_block=escape_taxa_as_context_block,
                                                               differentiate_inexact_match=differentiate_inexact_match)

                return self.ad_cluster_canvas
            else:
                return self.ad_cluster_canvas

    # tree_distance() : Function to get all trees' distance scatter ()
    def tree_distance(self,export=False):
        tsne = TSNE(n_components=2)
        self.tree_point_coordinates = tsne.fit_transform(self.tree_distance_matrix)

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

        self.tree_distance_fig = px.scatter(x=self.x_coor, y=self.y_coor, color=color,
                         title='Tree Distance')

        # Customize hover template
        self.tree_distance_fig.update_traces(hovertemplate='%{text}',text=id_name, hoverinfo='text')
        # Remove legend labels
        self.tree_distance_fig.update_layout(dragmode=False,showlegend=False,width=600, height=600,plot_bgcolor=TREE_NAME_BG,xaxis=dict(
            color=BLANK),yaxis=dict(color=BLANK))


        if not export:
            self.tree_distance_fig.show()

    # pairwise_comparison() : Function to show trees' pairwise comparison
    # Compare reference tree with tree from tree collection / Compare two trees from tree collection
    # Parameter : compare_tree(integer or integer list - tree id ; string or string list - tree name)
    def pairwise_comparison(self,compare_tree=None):
        if not self.rt_canvas and not self.pairwise_canvas:
            self.create_pairwise_rt_canvas()

        compare_between_tc = False
        if compare_tree:
            if compare_tree == self.pairwise_tree:
                display(self.pairwise_canvas)
                return
            else:
                self.pairwise_tree = compare_tree

            if type(compare_tree) is list:
                tmp_list = []
                for tree in compare_tree:
                    if isinstance(tree, str):
                        tmp_tree = self.get_tree_by_name(tree)

                        if tmp_tree is None:
                            self.tree_not_exist()
                            return

                        tmp_list.append(tmp_tree)

                    elif isinstance(tree, int):

                        if tree > len(self.tc) or tree < 0:
                            self.tree_not_exist()
                            return

                        tmp_list.append(self.tc[tree - 1])
                    else:
                        self.parameter_error("Tree info")
                        return

                compare_tree = tmp_list
                compare_between_tc = True

                if compare_tree[0] == compare_tree[1]:
                    self.tree_is_same_error()
                    return

            else:
                if isinstance(compare_tree, str):
                    compare_tree = self.get_tree_by_name(compare_tree)
                    if compare_tree is None:
                        self.tree_not_exist()
                        return
                else:
                    if compare_tree > len(self.tc) or compare_tree < 0:
                        self.tree_not_exist()
                        return

                    compare_tree = self.tc[compare_tree - 1]

        else:
            if not self.ad_individual_canvas.tree_selected:
                self.no_tree_chosen_error()
                return

            compare_tree = self.ad_individual_canvas.tree_selected.tc_tree


        display(self.pairwise_canvas)
        self.pairwise_canvas.compare_tc_tree(compare_tree,compare_between_tc=compare_between_tc)

    def tree_distribution(self,test=False):
        self.rt_canvas.clear_subtree_compare_canvas()
        self.tc_canvas_tmp = tcCanvas(adPy=self, view=TREE_DISTRIBUTION, subtree_independent=False)

        self.tree_distribution_view = TreeDistributionView(adPy=self,
                                                                                tc_canvas_tmp=self.tc_canvas_tmp,
                                                                                test=test)

    def export_image(self,filename='new_image',filetype='png',view=RT,tree_id=None,tree_name=None,ad_interval=None,
                     context_level=2,
                     sort_by=ID,show_tree_name=None,scale=1.0,escape_taxa_as_context_block=True,show_block_proportional=True,
           subtree_independent=False,differentiate_inexact_match=True,max_ad=None):

        if not self.check_alphanumeric_underscore(filename):
            self.filename_error()
            return

        self.image_name = filename
        self.image_type = filetype
        self.export_canvas = None

        if view == RT:
            if not self.rt_canvas:
                self.not_exist_error("Reference Tree View", "reference_tree()")

            self.export_canvas = self.get_rt_canvas_image()

        elif view == AD_INDIVIDUAL:
            if not self.ad_individual_canvas or self.check_ad_parameter(tree_id=tree_id,tree_name=tree_name,\
                    ad_interval=ad_interval,
                     sort_by=sort_by,show_tree_name=show_tree_name,scale=scale,
                                       escape_taxa_as_context_block=escape_taxa_as_context_block,
                                       show_block_proportional=show_block_proportional,
                                       subtree_independent=subtree_independent,context_level=context_level,
                                                                        max_ad=max_ad) == CHANGED:

                self.AD(view=AD_INDIVIDUAL,tree_id=tree_id,tree_name=tree_name,ad_interval=ad_interval,
                     sort=sort_by,show_tree_name=show_tree_name,scale=scale,
                                       escape_taxa_as_context_block=escape_taxa_as_context_block,
                                       show_block_proportional=show_block_proportional,
                                       subtree_independent=subtree_independent,context_level=context_level,
                        max_ad=max_ad,export=True)
            else:
                self.ad_individual_canvas_export = self.ad_individual_canvas


            if not self.ad_individual_canvas_export:
                return

            self.export_canvas = self.get_ad_canvas_image()


        elif view == AD_CLUSTER:
            if not self.ad_cluster_canvas or self.check_cluster_parameter(scale=scale,context_level=context_level,
                                       escape_taxa_as_context_block=escape_taxa_as_context_block,
                                       show_block_proportional=show_block_proportional,
                                       subtree_independent=subtree_independent,
                                            differentiate_inexact_match=differentiate_inexact_match) == CHANGED:
                self.AD(view=AD_CLUSTER,scale=scale,context_level=context_level,
                                       escape_taxa_as_context_block=escape_taxa_as_context_block,
                                       show_block_proportional=show_block_proportional,
                                       subtree_independent=subtree_independent,
                                            differentiate_inexact_match=differentiate_inexact_match,export=True)

            else:
                self.ad_cluster_canvas_export = self.ad_cluster_canvas

            if not self.ad_individual_canvas_export:
                return

            self.export_canvas = self.get_cluster_canvas_image()

        elif view == TREE_DISTRIBUTION:
            if not self.tree_distribution_view or not self.tree_distribution_view.export_ready:
                self.not_exist_error("Tree Distribution View","tree_distribution()")
                return

            self.export_canvas = self.get_tree_distribution_image()

        elif view == PAIRWISE_COMPARISON:
            if not self.pairwise_canvas:
                self.not_exist_error("Pairwise Comparison View", "pairwise_comparison()")
                return

            self.export_canvas = self.get_pairwise_comparison_image()

        elif view == TREE_DISTANCE:
            if len(self.tc) == 0:
                print("Tree Collection Not Exist")
                return

            if not self.tree_distance_fig:
                self.tree_distance(export=True)

            colors = ['blue'] + ['red'] * (len(self.x_coor) - 1)

            # 绘制散点图
            plt.scatter(self.x_coor, self.y_coor, color=colors)
            plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                            labelbottom=False, labelleft=False)
            plt.title('Tree Distance')
            plt.savefig(f"{self.image_name}.{self.image_type}")

            plt.close()


        if not self.export_canvas and view != TREE_DISTANCE:
            print("<Error> : View type incorrect")
            return

        if self.export_canvas:
            self.export_canvas.observe(self.save_to_file, "image_data")


    ##-------------------------- Internal Functions -----------------------------
    #
    def sort_tc(self,tc_list,sort_by):
        if sort_by == ID:
            return sorted(tc_list, key=lambda x: (x.id),reverse=False)
        elif sort_by == RF_DISTANCE:
            return sorted(tc_list, key=lambda x: (x.rf_distance, x.name),reverse=False)
        elif sort_by == NAME:
            return sorted(tc_list, key=lambda x: (x.name),reverse=False)

    def get_pairwise_comparison_image(self):
        pairwise_canvas_image = Canvas(width=self.pairwise_canvas.width, height=self.pairwise_canvas.height,
                                 sync_image_data=True)
        pairwise_canvas_image.fill_style = BLANK
        pairwise_canvas_image.fill_rect(0, 0, self.pairwise_canvas.width, self.pairwise_canvas.height)


        if self.pairwise_canvas.compare_between_tc:
            start_index = 7
            end_index = 20
        else:
            start_index = 0
            end_index = 13

        with hold_canvas(pairwise_canvas_image):
            for i in range(start_index,end_index + 1):
                pairwise_canvas_image.draw_image(self.pairwise_canvas[i], 0, 0)

        pairwise_canvas_image.flush()
        return pairwise_canvas_image

    # Ready canvas for tree dirstribution
    def get_tree_distribution_image(self):
        view_tmp = self.tree_distribution_view
        total_height = (view_tmp.nodes_list_canvas.height + view_tmp.related_tc_tree_canvas.height +
                               view_tmp.cluster_canvas.height + DEFAULT_PADDING_BETWEEN_AD * 2)

        if total_height < view_tmp.rt_subtree_block_canvas.height:
            total_height = view_tmp.rt_subtree_block_canvas.height

        total_height += 50
        total_width = view_tmp.nodes_list_canvas.width + view_tmp.rt_subtree_block_canvas.width + 50

        tree_distribution_canvas = Canvas(width=total_width,height=total_height,sync_image_data=True)
        tree_distribution_canvas.fill_style = BLANK
        tree_distribution_canvas.fill_rect(0,0,total_width,tree_distribution_canvas.height)
        subtree_chosen = view_tmp.subtree_chosen

        if view_tmp.agree_rt:
            label_str = (f"This cluster(#trees={view_tmp.segment_button_clicked.description}) agrees with branch"
                         f" {subtree_chosen.label} in the reference tree.")
        else:
            label_str = (f"This cluster(#trees={view_tmp.segment_button_clicked.description}) disagrees with branch"
                         f" {subtree_chosen.label} in the reference tree.")

        with hold_canvas(tree_distribution_canvas):
            pointer_x = 20
            pointer_y = 20
            rt_subtree_block_canvas_x = pointer_x + view_tmp.nodes_list_canvas.width + 10
            tree_distribution_canvas.fill_style = BLACK
            tree_distribution_canvas.font = f'18px Times New Roman'
            tree_distribution_canvas.fill_text(label_str,pointer_x,pointer_y)

            pointer_y += 20
            tree_distribution_canvas.draw_image(view_tmp.related_tc_tree_canvas, pointer_x, pointer_y)
            tree_distribution_canvas.draw_image(view_tmp.rt_subtree_block_canvas, rt_subtree_block_canvas_x, pointer_y)
            pointer_y += view_tmp.related_tc_tree_canvas.height + 10

            tree_distribution_canvas.draw_image(view_tmp.cluster_canvas, pointer_x, pointer_y)
            pointer_y += view_tmp.cluster_canvas.height + 10

            tree_distribution_canvas.draw_image(view_tmp.nodes_list_canvas, pointer_x, pointer_y)
            pointer_y += view_tmp.nodes_list_canvas.height + 10

        tree_distribution_canvas.flush()

        return tree_distribution_canvas

    def check_ad_parameter(self,tree_id=None,tree_name=None,ad_interval=None,scale=1.0,
                     sort_by=None,show_tree_name=None,escape_taxa_as_context_block=True,show_block_proportional=True,
           subtree_independent=False,context_level=2,max_ad=None):
        canvas = self.ad_individual_canvas
        if not tree_id and not tree_name and not ad_interval and not sort_by and not show_tree_name and not max_ad:
            if (not canvas.target_tc_tree_index and not canvas.target_tc_tree_name and not canvas.first_ad and not
            canvas.last_ad and not canvas.sort_by and not canvas.show_tree_name):
                if (context_level == canvas.context_level and scale == canvas.scale and escape_taxa_as_context_block ==
                        canvas.escape_taxa_as_context_block and
                        show_block_proportional ==
                        canvas.show_block_proportional and subtree_independent == canvas.subtree_independent):

                    return NO_CHANGE
                else:
                    return CHANGED
            else:
                return CHANGED

        if tree_id and tree_id != canvas.target_tc_tree_index:
            return CHANGED

        if tree_name and tree_name != canvas.target_tc_tree_name:
            return CHANGED

        if ad_interval and ad_interval[0] != canvas.first_ad and ad_interval[1] != canvas.last_ad:
            return CHANGED

        if sort_by and sort_by != canvas.sort_by:
            return CHANGED

        if show_tree_name and show_tree_name != canvas.show_tree_name:
            return CHANGED

        if max_ad and max_ad != canvas.max_ad:
            return CHANGED

    def check_cluster_parameter(self, context_level=2,scale=1.0,escape_taxa_as_context_block=True,show_block_proportional=True,
           subtree_independent=False,differentiate_inexact_match=True):
        canvas = self.ad_cluster_canvas

        if (context_level == canvas.context_level and scale == canvas.scale and escape_taxa_as_context_block ==
                canvas.escape_taxa_as_context_block and show_block_proportional == canvas.show_block_proportional and
                subtree_independent == canvas.subtree_independent and differentiate_inexact_match == canvas.differentiate_inexact_match):

            return NO_CHANGE
        else:
            return CHANGED

    # Ready canvas for cluster ad
    def get_cluster_canvas_image(self):
        if not self.ad_cluster_canvas_export:
            return None
        cluster_export_canvas = Canvas(width=self.ad_cluster_canvas_export.width, height=
        self.ad_cluster_canvas_export.height, sync_image_data=True)
        cluster_export_canvas.fill_style = BLANK

        with hold_canvas(cluster_export_canvas):
            cluster_export_canvas.fill_rect(0, 0, cluster_export_canvas.width, cluster_export_canvas.height)
            cluster_export_canvas.draw_image(self.ad_cluster_canvas_export[-3], 0, 0)

        cluster_export_canvas.flush()
        return cluster_export_canvas

    # Ready canvas for tree collection ad
    def get_ad_canvas_image(self):
        if not self.ad_individual_canvas_export:
            return None
        ad_export_canvas = Canvas(width=self.ad_individual_canvas_export.width,height =
        self.ad_individual_canvas_export.height,sync_image_data=True)
        ad_export_canvas.fill_style = BLANK


        with hold_canvas(ad_export_canvas):
            ad_export_canvas.fill_rect(0, 0, ad_export_canvas.width, ad_export_canvas.height)
            ad_export_canvas.draw_image(self.ad_individual_canvas_export[-3],0,0)

        ad_export_canvas.flush()

        return ad_export_canvas

    #  Ready canvas for rt export
    def get_rt_canvas_image(self):
        rt_canvas_image = Canvas(width = self.rt_canvas.width,height = self.rt_canvas.height,sync_image_data=True)
        rt_canvas_image.fill_style = BLANK

        with hold_canvas(rt_canvas_image):
            rt_canvas_image.fill_rect(0, 0, self.rt_canvas.width, self.rt_canvas.height)

            for i in range(6):
                rt_canvas_image.draw_image(self.rt_canvas[i],0,0)

            rt_canvas_image.draw_image(self.rt_canvas[-1], 0, 0)


        rt_canvas_image.flush()

        return rt_canvas_image

    def save_to_file(self,*args, **kwargs):
        # self.export_canvas.to_file("save_file.png")
        self.output.append("in save_to_file")
        self.export_canvas.to_file(f"{self.image_name}.{self.image_type}")

    def print_tree_distribution(self):
        if not self.ad_individual_canvas and not self.ad_cluster_canvas and not self.tree_distribution_tmp:
            tmp = tcCanvas(adPy=self, view=TREE_DISTRIBUTION)

        for subtree in self.subtree_list:
            print(f"Subtree {subtree.label}")
            index = 1
            print(subtree.topology_list)
            for topology,tree_list in subtree.topology_list.items():
                print(f"Topology {index} : {len(tree_list)}")
                print(topology)
                index += 1

    def read_rt(self,treefile,type):
        self.rt = dendropy.Tree.get(path=treefile, schema=type)

        # Filename as tree name
        rt_label = os.path.basename(treefile)
        self.rt.id = 0
        self.rt.name = rt_label
        self.rt.level = 0
        self.rt.pairwise_canvas = None
        self.rt.missing = set()

        self.rt.taxa_list = []
        self.rt.internal_node = []

        self.rt.taxa_list = [leaf.taxon.label for leaf in self.rt.leaf_nodes()]



    def generate_missing_node(self,tree,taxa_name):
        new_taxon = dendropy.Taxon(label=taxa_name)
        new_node = dendropy.Node(taxon=new_taxon, label=taxa_name)
        new_node.is_missing = True

        tree.missing_node_list.append(new_node)

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
                missing = tree.missing - set(node_taxa)

                node.corr_similarity = 0
                node.corr.append(0)

                for tc_node in tree.postorder_node_iter():
                    if not hasattr(tc_node,'corr'):
                        tc_node.corr = None
                        tc_node.corr_similarity = 0

                    similarity = self.get_similarity(target_set=target_set,node=node,tc_node=tc_node,missing=missing)
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

                # If no corresponding branch, if is missing taxa in target tree
                if node.is_leaf() and node.corr[tree.id] == 0:
                    for tc_node in tree.missing_node_list:
                        if tc_node.label == node.taxon.label:
                            node.corr[tree.id] = tc_node
                            node.corr_similarity = 1.0
                            tc_node.corr = node
                            tc_node.corr_similarity = 1.0

                            break

    def get_similarity(self,target_set,node,tc_node,missing):
        if tc_node.is_leaf():
            tc_node.card = 0 if tc_node.taxon.label in missing else 1
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


    def not_exist_error(self,tree,pre_function):
        print(f"<Error> : {tree} Not Exist.")
        print(f"Please ensure that you have called the {pre_function} function before calling this function.")

    def parameter_error(self,parameter):
        print(f"<Error> : {parameter} given was incorrect.")
        print("Please ensure that the information provided is logical.")


    def no_tree_chosen_error(self):
        print(f"<Error> : No tree was chosen nor given.")
        print("Please choose a tree from AD() OR provide tree id/tree name.")

    def tree_not_exist(self):
        print(f"<Error> : Tree is not exist.")
        print("Please choose an existing tree.")

    def filename_error(self):
        print(f"<Error> : Please ensure that the filename contains only alphanumeric characters and underscores.")

    def tree_is_same_error(self):
        print(f"<Error> : Tree chosen is same.")
        print("Please choose different trees.")

    def check_alphanumeric_underscore(self,str):
        # 定义正则表达式模式，表示字符串只包含字母、数字和下划线
        pattern = r'^[a-zA-Z0-9_]+$'
        # 使用 re.match() 函数进行匹配
        if re.match(pattern, str):
            return True
        else:
            return False



    # Manage Subtree
    def default_subtree_attribute(self):
        self.subtree_list = []
        self.subtree_label_used = [1, 1, 1, 1, 1]
        self.subtree_color_used = [1, 1, 1, 1, 1]

    def create_pairwise_rt_canvas(self,alter_type=BOTH):
        height = get_leaf_node_amount(self.rt) * RT_Y_INTERVAL + RT_Y_INTERVAL * 5

        if self.rt_alter:
            self.default_rt = None
            self.rt_alter = False

        # Check paramater alter
        if alter_type == BOTH:
            # self.default_rt_value()
            self.rt_canvas = rtCanvas(self, width=CANVAS_MAX_WIDTH, height=height,
                                      view_support=self.rt_view_support, default_rt=self.default_rt,
                                      exact_match_range=self.rt_exact_match_range,
                                      support_value_range=self.rt_support_value_range)
            self.pairwise_canvas = pairwiseCanvas(self, width=CANVAS_MAX_WIDTH, height=height)
        elif alter_type == RT:
            self.rt_canvas = rtCanvas(self, width=CANVAS_MAX_WIDTH, height=height,
                                      view_support=self.rt_view_support, default_rt=self.default_rt,
                                      exact_match_range=self.rt_exact_match_range,
                                      support_value_range=self.rt_support_value_range)
        elif alter_type == PAIRWISE:
            self.pairwise_canvas = pairwiseCanvas(self, width=CANVAS_MAX_WIDTH, height=height)


    def check_parameter_alter(self,view_support):
        # Check whether parameter has alter and set parameter as self.attribute
        if view_support != self.rt_view_support:
            self.rt_view_support = view_support
            return True

    def check_attribute_range(self,interval_list):

        if len(interval_list) != 2:
            return False

        if not isinstance(interval_list[0], int) or not isinstance(interval_list[1], int):
            return False

        if interval_list[0] > interval_list[1]:
            return False

        return True

    def select_subtree_from_tree(self,node_selected):
        if not hasattr(node_selected, 'selected') or node_selected.selected == False:
            self.ad_parameter_alter = True
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
            self.pairwise_canvas.draw_subtree_block(node_selected,select_tree=RT,new_subtree=new_subtree,align=LEFT)

            if self.pairwise_canvas.tc_tree:
                self.pairwise_canvas.draw_tc_subtree_block(new_subtree,align=RIGHT)

                self.pairwise_canvas.draw_escape_taxa(align=RIGHT)
                if self.pairwise_canvas.compare_between_tc:
                    self.pairwise_canvas.draw_tc_subtree_block(new_subtree,align=LEFT)
                    self.pairwise_canvas.draw_escape_taxa(align=LEFT)

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

    def get_subtree(self,label):
        for subtree in self.subtree_list:
            if subtree.label == label:
                return subtree

    def get_tree_by_name(self,name):
        for tree in self.tc:
            if tree.name == name:
                return tree

        return None



def print_tree(tree):
    print(tree.as_ascii_plot())

def get_node_list(tree):
    return tree.leaf_node_iter()

def get_leaf_node_amount(tree):
    return sum(1 for node in tree.leaf_node_iter())


