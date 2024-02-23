import dendropy
from Utils import *
import os
import myCanvas as myCanvas
import importlib
importlib.reload(myCanvas)

class myTree:
    rt = None  # Reference Tree
    tc = None  # Tree Collection
    rt_canvas = None # Reference Tree's Canvas
    ad_individual_canvas = None # Individual ADs' Canvas
    ad_cluster_canvas = None # Cluster ADs' Canvas
    rt_view_support = False # Whether to show internal node's support value
    ad_parameters = {}   # { 'parameter_name' : value }

    def __init__(self,treefile = None,type="newick"):
        self.read_rt(treefile = treefile, type = type)

    def read_rt(self,treefile,type):
        self.rt = dendropy.Tree.get(path=treefile, schema=type)
        rt_label = os.path.basename(treefile)
        self.rt.label = rt_label

    def reference_tree(self,view_support=False):
        # Calculate height of canvas
        height = get_leaf_node_amount(self.rt) * Y_INTERVAL + Y_INTERVAL

        if not self.rt_canvas or self.rt_canvas.view_support != view_support:
            self.rt_canvas = myCanvas.rtCanvas(self.rt,width = CANVAS_MAX_WIDTH,height = height,view_support = view_support)

        return self.rt_canvas

    def set_outgroup(self,outgroup_taxon):
        # ["Uronema sp", "Monomastix opisthostigma", "Pyramimonas parkeae", "Nephroselmis pyriformis"]
        mrca = self.rt.mrca(taxon_labels=outgroup_taxon)
        self.rt.reroot_at_edge(mrca.edge, update_bipartitions=False)

    def add_tree_collection(self,treefile=None,type="newick",namefile=None):
        self.tc = dendropy.TreeList.get_from_path(treefile, schema=type)

        # Read tree collection names
        if namefile:
            file = open(namefile, "r")
            for tree in self.tc:
                tree_name = file.readline().strip()
                tree.label = tree_name

    # Same effect as click on the common ancestor of these nodes
    def select_subtree(self,nodes=None):
        nodes_list = []
        subtree_root = None

        # Approximate taxon name
        for node in nodes:
            for leaf_node in self.rt.leaf_node_iter():
                if node.lower() in leaf_node.taxon.label.lower():
                    if len(nodes) < 2:
                        self.rt_canvas.draw_subtree_block(leaf_node.parent_node)
                        return

                    nodes_list.append(leaf_node.taxon.label)

        subtree_root = self.rt.mrca(taxon_labels=nodes_list)

        self.rt_canvas.draw_subtree_block(subtree_root)

    # Same effect as click on the leaf_node
    def select_leaf_node(self, node=None):
        # Approximate taxon name
        for leaf_node in self.rt.leaf_node_iter():
            if node.lower() in leaf_node.taxon.label.lower():
                self.rt_canvas.draw_subtree_block(leaf_node)

        # Exact taxon name
        # node = self.rt.find_node_with_taxon_label(node)
        # self.rt_canvas.draw_subtree_block(node)

    # Show AD, default: show individual AD
    def AD(self,view=INDIVIDUAL,):  # view = "Individual" / "Cluster"
        # Parameter: width,height,context levels, show labels, show colors
        if view == INDIVIDUAL:
            if not self.ad_individual_canvas or self.parameter_modified():
                pass


        # Parameter: differentiate inexact match, differentiate sister-group relationships
        elif view == CLUSTER:
            if not self.ad_cluster_canvas or self.parameter_modified():
                pass

    def get_tc(self):
        return self.tc

    def corresponding_branches(self):
        if self.tc == None:
            print("Tree Collection Not Exist")
            return



class Subtree:
    rt = None  # Belong to which reference tree
    label = None  # A,B,C,D,E
    rtCanvas_index = None
    root = None # Subtree root
    color = None
    block = None
    label_width = None

    def __init__(self,label,rt,root,color,block):
        self.label = label
        self.rt = rt
        self.root = root
        self.color = color
        self.block_size = root.block.get_size()
        self.block = block

    def set_rtLayer(self,rtCanvas_index):
        self.rtCanvas_index = rtCanvas_index




def print_tree(tree):
    print(tree.as_ascii_plot())

def get_node_list(tree):
    return tree.leaf_node_iter()

def get_leaf_node_amount(tree):
    return sum(1 for node in tree.leaf_node_iter())

'''   Test Code  '''

'''
mytree = myTree(treefile = "Data/69species/astral.FAA.trim50genes.final.tre")
mytree.add_tree_collection(treefile = "Data/69species/tc.tre",namefile="Data/69species/names.txt")
tree = mytree.rt

# for edge in tree.postorder_edge_iter():
#     print(edge)
# get_adjacent_edges()
# te = [i for i in self.tail_node.incident_edges() if i is not self]
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