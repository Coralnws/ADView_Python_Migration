import dendropy
from Utils import *
import os
import myCanvas as myCanvas
import importlib
importlib.reload(myCanvas)
import math

class myTree:
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

        # AD/Tree collection related
        self.tc = []  # Tree Collection
        self.tc_taxa = []
        self.ad_individual_canvas = None # Individual ADs' Canvas
        self.ad_cluster_canvas = None # Cluster ADs' Canvas
        self.ad_parameters = {}   # { 'parameter_name' : value }
        self.ad_parameter_alter = True # Record whether ad_canvas's parameter change

        # Reaf tree file and construct reference tree
        self.read_rt(treefile = treefile, type = type)


    #####  Pubic Functions ####
    def reference_tree(self,view_support=False):
        # Calculate height of canvas
        height = get_leaf_node_amount(self.rt) * RT_Y_INTERVAL + RT_Y_INTERVAL

        if not self.rt_canvas or self.rt_canvas.view_support != view_support or self.rt_alter:
            self.rt_canvas = myCanvas.rtCanvas(self,width = CANVAS_MAX_WIDTH,height = height,view_support = view_support,tc = self.tc)
            self.subtree_list = self.rt_canvas.subtree_list

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

        self.rt_alter = True

    def add_tree_collection(self,treefile=None,type="newick",namefile=None):
        self.tc = dendropy.TreeList.get_from_path(treefile, schema=type)

        # Read tree collection name
        if namefile:
            file = open(namefile, "r")

        for index,tree in enumerate(self.tc):
            # Set trees' index, taxa_list, missing taxa, name and outgroup
            tree.id = index
            tree.taxa_list = [leaf.taxon.label for leaf in tree.leaf_nodes()]
            tree.missing = set(self.rt.taxa_list) - set(tree.taxa_list)

            if namefile:
                tree_name = file.readline().strip()
                tree.name = tree_name

            if len(self.outgroup) > 1:
                mrca = tree.mrca(taxon_labels=self.outgroup)
                tree.reroot_at_edge(mrca.edge)

        self.corresponding_branches()

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
    def select_leaf_node(self, node=None):
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
    def AD(self,view=AD_INDIVIDUAL,scale=1.0,max_ad=None,context_level=2,first_ad=None,last_ad=None,tree_index=None,tree_name=None):  # view = "Individual" / "Cluster"
        # Check if condition is logical
        # 1. First_ad < last_ad
        if first_ad and last_ad:
            if last_ad < first_ad:
                self.parameter_error("first_ad/last_ad")
                return
        # 2.  First_ad > 0
        if first_ad and first_ad <= 0:
            self.parameter_error("first_ad")
            return
        # 3. Last_ad > 0
        if last_ad and last_ad <= 0:
            self.parameter_error("last_ad")
            return
        # 4. tree_index > 0
        if tree_index and type(tree_index) is not list and tree_index <= 0:
            self.parameter_error("tree_index")
            return
        # 5. context_level >= 1
        if context_level < 1:
            self.parameter_error("context_level")
            return

        self.subtree_list = self.rt_canvas.subtree_list
        # Parameter: width, height, context levels, show labels, show colors

        if view == AD_INDIVIDUAL:
            # if not self.ad_individual_canvas or self.parameter_modified():
            if not self.ad_individual_canvas or self.ad_parameter_alter:
                if scale != 1.0:
                    ad_per_row = (CANVAS_MAX_WIDTH - (3 * DEFAULT_AD_PADDING)) // (DEFAULT_AD_WIDTH * scale)
                else:
                    ad_per_row = DEFAULT_AD_PER_ROW

                if max_ad:
                    ad_row = math.ceil(max_ad / ad_per_row)
                else:
                    ad_row = math.ceil(len(self.tc) / ad_per_row)

                canvas_height = (ad_row * DEFAULT_AD_HEIGHT * scale) + (2 * DEFAULT_AD_PADDING) + ((ad_row-1) * DEFAULT_AD_PADDING)
                self.tc_canvas = myCanvas.tcCanvas(tree=self,width=CANVAS_MAX_WIDTH,height=canvas_height,scale=scale,max_ad=max_ad,\
                                                   ad_per_row=ad_per_row,context_level=context_level,first_ad=first_ad,last_ad=last_ad,\
                                                   tree_index=tree_index,tree_name=tree_name)

                return self.tc_canvas

        # Parameter: differentiate inexact match, differentiate sister-group relationships
        elif view == AD_CLUSTER:
            if not self.ad_cluster_canvas or self.parameter_modified():
                pass

    #### Internal Functions ####
    def read_rt(self,treefile,type):
        self.rt = dendropy.Tree.get(path=treefile, schema=type)

        # Filename as tree name
        rt_label = os.path.basename(treefile)
        self.rt.label = rt_label

        # Get tree's taxa/leaf node
        self.rt.taxa_list = [leaf.taxon.label for leaf in self.rt.leaf_nodes()]

    def get_tc(self):
        return self.tc

    def corresponding_branches(self):
        if self.tc == None:
            print("Tree Collection Not Exist")
            return

        for node in self.rt.postorder_node_iter():
            node_taxa = [leaf.taxon.label for leaf in node.leaf_nodes()]
            node.corr = []
            node.exact_match = 0
            node.exact_match_tree = []

            for tree in self.tc:
                target_set = set(node_taxa) - tree.missing
                self.rt.missing = tree.missing - set(node_taxa)

                node.corr_similarity = 0
                node.corr.append(0)

                for tc_node in tree.postorder_node_iter():
                    similarity = self.get_similarity(target_set=target_set,node=node,tc_node=tc_node)
                    if similarity > node.corr_similarity:
                        node.corr[tree.id] = tc_node
                        node.corr_similarity = similarity

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

        if node == child_node:
            return True

        return False

    def rt_not_exist_error(self):
        print("<Error> : Reference Tree not exist.")
        print("Please ensure that you have called the reference_tree() function before calling this function.")

    def parameter_error(self,parameter):
        print(f"<Error> : {parameter} given was incorrect.")
        print("Please ensure that the information provided is logical.")

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

    def __init__(self,label,rt,root,color,block):
        self.label = label
        self.rt = rt
        self.root = root
        self.color = color

        if block:
            self.block_size = root.block.get_size()
            self.block = block

        self.label_width = None
        self.leaf_set = set()


    def set_rtLayer(self,rtCanvas_index):
        self.rtCanvas_index = rtCanvas_index

    def get_leaf_nodes(self):
        self.leaf_set = set()
        for leaf_node in self.root.leaf_nodes():
            self.leaf_set.add(leaf_node.taxon.label)

# AD
class AD_Tree:
    # Internal node is AD_Block
    # Leaf node is AD_Block
    # Tree_width = self.width - padding
    padding = 5  # Padding between block

    def __init__(self,id,tc_tree,tc_canvas):
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
        self.set_default()
    def set_default(self):
        self.x = 8
        self.y = 8
        self.padding = 5  # Padding between block

    def set_position_size(self,x,y):
        self.topL = myCanvas.Point(x,y)
        self.width = self.tc_canvas.scale * DEFAULT_AD_WIDTH
        self.height = self.tc_canvas.scale * DEFAULT_AD_HEIGHT
        self.botR = myCanvas.Point(x+self.width,y+self.height)

    def set_nested_tree_size(self,nested_block):
        self.topL = nested_block.topL
        self.botR = nested_block.botR
        self.width = self.botR.x - self.topL.x

    def plot_tree(self, node=None, level=0):
        if node is None:
            node = self.root

        if node.type == ROOT or node.type == INTERNAL:
            internaL_node = node.node_or_block
            print(' ' * (level * 4) + node.type ,end=" ")
            print(internaL_node)

        elif node.type == LEAF:
            block = node.node_or_block

            if block.type == SUBTREE_BLOCK:
                print(' ' * (level * 4) + block.belong_subtree.label + ':' + str(block.taxa_count))
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

    def insert_node(self,parent,child):
        parent.children.append(child)

        # if child.type == LEAF:
        #     self.block_list.append(child.node_or_block)

    def traverse_postorder(self,node=None,node_list=None):
        if node == None:
            node = self.root
            node_list=[]

        node.children = sorted(node.children, key=lambda x: x.child_index,reverse=False)
        for child in node.children:
            node_list = self.traverse_postorder(child,node_list)

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
                    unnested_block_taxa += block.taxa_count

        for block in unnested_block_list:
            subtree_block_list.append(block)

        return indv_block_list,subtree_block_list,unnested_block_taxa

    # def get_subtree_block_list(self,order=BY_LEVEL):
    #     block_list = []
    #     if order == BY_LEVEL:
    #         for block in self.block_list:
    #             if block.type == SUBTREE_BLOCK or block.type == INDIVIDUAL_LEAF_BLOCK:
    #                 block_list.append(block)
    #
    #     return sorted(block_list, key=lambda x: x.x_level,reverse=False)

    # Calculate allocatable space for subtree block
    # If has missing taxa, just make changes in this function to reserve space for missing taxa block
    def allocatable_ad_space(self):
        indv_block_list,subtree_block_list = self.individual_subtree_block_list()
        indv_block_cnt = len(indv_block_list)
        subtree_block_cnt = len(subtree_block_list) + self.nested_cnt
        # print(subtree_block_cnt)
        total_block_cnt = len(self.block_list)
        scale = self.tc_canvas.scale

        if scale < 1.0:
            self.y = self.y * scale
            self.x = self.x * scale
            self.padding = self.padding * scale

        block_space = self.height - (2 * self.y) - (2 * self.padding) - (DEFAULT_INDV_BLOCK_HEIGHT * scale * indv_block_cnt) - \
                      ((total_block_cnt-1) * self.padding) - (BLOCK_MINIMUM_HEIGHT[subtree_block_cnt % MAX_SUBTREE] * scale * subtree_block_cnt)

        return block_space




    def ad_taxa_total(self):
        indv_block_list, subtree_block_list = self.individual_subtree_block_list()
        taxa_cnt = 0

        for block in subtree_block_list:
            taxa_cnt += block.taxa_count
            if block.nested_tree:
                taxa_cnt += block.nested_tree.ad_taxa_total()

        return taxa_cnt

    def ad_to_newick(self, ad_tree):
        pass

    def compare_topology(self):
        pass

class AD_Node:
    def __init__(self,node_or_block,x_level,type,child_index):
        self.node_or_block = node_or_block  # Dendropy node if internal node, AD_Block if leaf
        self.children = []
        self.x_level = x_level
        self.type = type   # LEAF or INTERNAL
        self.pos = None
        self.branch_head = None
        self.branch_tail = None
        self.nested_tree = None # Root of nested tree
        self.is_elided = False
        self.child_index = child_index

    def insert_child(self,ad_node):
        self.children.append(ad_node)

    def construct_branch(self,x,y,scale):
        self.branch_head = myCanvas.Point(x - (AD_BRANCH_LENGTH * scale),y)
        self.branch_tail = myCanvas.Point(x ,y)

    def construct_elided_branch(self,x,y,scale):
        self.branch_head = myCanvas.Point(x - (ELIDED_BRANCH_LENGTH * scale),y)
        self.branch_tail = myCanvas.Point(x ,y)

    # Standardize the starting x-coordinate of all children's branches
    def unify_children_branches(self,x):
        for child in self.children:
            child.branch_head.x = x

def print_tree(tree):
    print(tree.as_ascii_plot())

def get_node_list(tree):
    return tree.leaf_node_iter()

def get_leaf_node_amount(tree):
    return sum(1 for node in tree.leaf_node_iter())

'''   Test Code  '''



# for edge in tree.postorder_edge_iter():
#     print(edge)
# get_adjacent_edges()
# te = [i for i in self.tail_node.incident_edges() if i is not self]
'''
edges = mytree.rt.edges()
for edge in edges:
    tail_node = edge.tail_node
    head_node = edge.head_node
    if edge.is_internal():
        print("====Internal branch====")
        if tail_node and tail_node.is_leaf():
            print("Tail Node (Leaf):", tail_node.taxon.label)
        elif tail_node:
            print("Tail Node (Internal):", tail_node.label)
        if head_node and head_node.is_leaf():
            print("Head Node (Leaf):", head_node.taxon.label)
        elif head_node:
            print("Head Node (Internal):", head_node.label)


        # neighbor_edges = head_node.incident_edges()
        # print("   ---Neighbout edges---")
        # # 遍历相邻边并找到连接到叶节点的分支
        # for neighbor_edge in neighbor_edges:
        #     if head_node is not neighbor_edge.tail_node:
        #         continue
        #     edge_tail_node = neighbor_edge.tail_node
        #     edge_head_node = neighbor_edge.head_node
        #     if neighbor_edge.is_internal():
        #         print("is internal")
        #         if edge_tail_node and edge_tail_node.is_leaf():
        #             print("Tail Node (Leaf):", edge_tail_node.taxon.label)
        #         elif edge_tail_node:
        #             print("Tail Node (Internal):", edge_tail_node.label)
        #         if edge_head_node and edge_head_node.is_leaf():
        #             print("Head Node (Leaf):", edge_head_node.taxon.label)
        #         elif edge_head_node:
        #             print("Head Node (Internal):", edge_head_node.label)
        #     else:
        #         print("is leaf")
        #         if edge_tail_node.is_leaf():
        #             print("Tail Node (Leaf):", edge_tail_node.taxon.label)
        #         else:
        #             print("Tail Node (Internal):", edge_tail_node.label)
        #         if edge_head_node.is_leaf():
        #             print("Head Node (Leaf):", edge_head_node.taxon.label)
        #         else:
        #             print("Head Node (Internal):", edge_head_node.label)

    else:
        print("====Leaf====")
        if tail_node.is_leaf():
            print("Tail Node (Leaf):", tail_node.taxon.label)
        else:
            print("Tail Node (Internal):", tail_node.label)
        if head_node.is_leaf():
            print("Head Node (Leaf):", head_node.taxon.label)
        else:
            print("Head Node (Internal):", head_node.label)
'''



# for tree in mytree.tc:
#     taxa_names2 = set(taxon.label for taxon in tree.taxon_namespace)
#     print(taxa_names2)
#     missing_taxa = taxa_names1 - taxa_names2
#     print(len(missing_taxa))



"""
def filter_node_selected(node, x, y):
    if hasattr(node, 'x1'):
        if node.x1 <= x <= node.x2 and node.y1 <= y <= node.y2:
            print(node.taxon)
            return True
        return False


def add_coordinate(node):
    if node.is_leaf() and node.taxon.label == "Uronema sp":
        node.x1 = 10
        node.x2 = 20
        node.y1 = 10
        node.y2 = 20

mytree = myTree(treefile = "../../Data/69species/astral.FAA.trim50genes.final.tre")
print(get_leaf_node_amount(mytree.rt))
mytree.rt.find_node(lambda node: add_coordinate(node))
node_selected = mytree.rt.find_node(lambda node: filter_node_selected(node, 13, 13))
print(node_selected)
"""