from ipycanvas import Canvas,MultiCanvas
from myUtils import *
import math

# Use Multicanvas which has multiple canvas layer
class MyCanvas(MultiCanvas):
    tree = None
    line_dashes = [[2, 2], [5, 5], [8, 8], [10, 10]]

    def __init__(self,layer,width, height):
        super().__init__(layer,width = width, height = height, sync_image_data=True)
        self.layer = layer

    def draw_dots(self,x,y,layer_index = -1,canvas=None):
        if not canvas:
            canvas = self[layer_index]
        canvas.fill_style = RED
        canvas.fill_circle(x, y, 1)
        canvas.fill_style = BLACK

    def draw_solid_line(self,head,tail,color,layer_index = -1,canvas=None):
        if not canvas:
            canvas = self[layer_index]
        canvas.stroke_style = color
        canvas.begin_path()
        canvas.move_to(head.x, head.y)
        canvas.line_to(tail.x,tail.y)
        canvas.stroke()

    def draw_dotted_line(self,head,tail,color,line_dash_index=0,layer_index = -1,canvas=None):
        if not canvas:
            canvas = self[layer_index]
        canvas.stroke_style = color
        canvas.set_line_dash(self.line_dashes[line_dash_index])
        canvas.begin_path()
        canvas.move_to(head.x, head.y)
        canvas.line_to(tail.x,tail.y)
        canvas.stroke()

        canvas.set_line_dash([])

    # Draw rectangle on specific layer
    def draw_rec(self,x,y,width,height,color,layer_index = -1,canvas=None):
        if not canvas:
            canvas = self[layer_index]
        canvas.fill_style = color
        canvas.fill_rect(x, y, width,height)

    def draw_frame(self,x,y,width,height,color,layer_index = -1,canvas=None):
        if not canvas:
            canvas = self[layer_index]
        canvas.stroke_style = color
        canvas.stroke_rect(x, y, width, height)

    def draw_dotted_frame(self,x,y,width,height,color,line_dash_index=0,layer_index = -1,canvas=None):
        if not canvas:
            canvas = self[layer_index]
        canvas.stroke_style = color
        canvas.set_line_dash(self.line_dashes[line_dash_index])
        canvas.stroke_rect(x, y, width, height)

        canvas.set_line_dash([])

    def draw_elided_branch(self,head,tail,color,layer_index = -1,canvas=None):
        self.draw_dotted_line(head,tail,color,line_dash_index=0,layer_index=layer_index,canvas=canvas)

        x = head.x + (tail.x - head.x) / 2
        y = tail.y - 5
        start_point = Point(x + 1,y)
        end_point = Point(x - 2,y + 10)

        self.draw_solid_line(head=start_point,tail=end_point,color=color,layer_index=layer_index,canvas=canvas)

    def crop_and_paste_canvas(self, paste_canvas,layer_index,copy_list,paste_list):
        x = copy_list[0]
        y = copy_list[1]
        width = copy_list[2]
        height = copy_list[3]

        crop = self[layer_index].get_image_data(x,y,width,height)

        x = paste_list[0]
        y = paste_list[1]

        paste_canvas.put_image_data(crop, x,y)
        # copy_canvas.observe(crop_and_paste_canvas, "image_data")

    def get_subtree(self,subtree_list,label):
        for subtree in subtree_list:
            if subtree.label == label:
                return subtree

        return None

class Point():
    def __init__(self,x,y):
        self.x = x
        self.y = y

class Block():
    def __init__(self,topL,botR):
        self.topL = topL
        self.botR = botR
        self.branch_x = None
        self.segment_list = []

        if topL and botR:
            self.calculate_width_height()
            self.topR = Point(self.botR.x, self.topL.y)
            self.botL = Point(self.topL.x, self.botR.y)

    def check_in_range(self,point):
        return point.x >= self.topL.x and point.x <= self.botR.x and point.y >= self.topL.y and point.y <= self.botR.y

    def print_block(self):
        print("topL = (",self.topL.x,",",self.topL.y,")")
        print("botR = (",self.botR.x,",",self.botR.y,")")

    def get_size(self):
        return self.botR.y - self.topL.y

    def calculate_width_height(self):
        self.width = self.botR.x - self.topL.x
        self.height = self.botR.y - self.topL.y

    def check_nested_block(self,block,select_tree=None):
        if not select_tree or select_tree == RT:
            return self.check_in_range(block.topL) and block.botR.y <= self.botR.y
        elif select_tree == TC:
            return self.check_in_range(block.botR) and block.botR.y <= self.botR.y

# Block for ad_tree
class AD_Block(Block):
    # id = None  # Block id in ad_tree

    def __init__(self, width,height,subtree,ad_tree,type,topL=None, botR=None):
        super().__init__(topL, botR)
        self.root = None  # To record this block represent which branch. If represent a subtree,root is subtree root
        self.taxa_list = set()
        self.subtree_taxa_count = 0   # subtree_taxa in block
        self.duplicate_subtree_taxa_count = []  # taxa_count of each subtree in same block
        self.exact_match = None
        self.color = None
        self.belong_subtree = subtree
        self.belong_ad_tree = ad_tree
        self.nested_tree = None
        self.is_elided = False
        self.is_subtree_duplicate = False

        self.type = type
        self.width = width
        self.height = height


class Missing_Taxa_Segment(Block):
    def __init__(self, subtree, ad_tree, missing_taxa,topL=None, botR=None):
        super().__init__(topL, botR)
        self.ad_tree = ad_tree
        self.subtree = subtree # Missing taxa from which subtree
        self.missing_taxa_list = missing_taxa # Include which taxa

class Subtree_Block_Segment(Block):
    def __init__(self, subtree,belong_block,topL=None, botR=None):
        super().__init__(topL, botR)
        self.subtree = subtree
        self.belong_block = belong_block
        self.belong_block.segment_list.append(self)


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

    def subtree_to_string(self, node):
        if len(node.child_nodes()) == 0:
            return node.taxon.label
        else:
            child_strings = [self.subtree_to_string(child) for child in node.child_nodes()]
            return "(" + ",".join(child_strings) + ")"

    def check_and_set_topology(self, root, tc_tree_id, type=NODE_DIFFERENCE):
        if type == TOPOLOGY_DIFFERENCE:
            subtree_string = self.subtree_to_string(root)

        elif type == NODE_DIFFERENCE:
            children = [child.taxon.label for child in root.leaf_nodes()]
            nodes_list = root.leaf_nodes()
            children = sorted(children)
            subtree_string = ",".join(children)

        if tc_tree_id == 0:
            return

        if subtree_string not in self.topology_list:
            self.topology_list[f"{subtree_string}"] = [nodes_list, []]

        self.topology_list[f"{subtree_string}"][1].append(tc_tree_id)

    def set_rtLayer(self, rtCanvas_index):
        self.rtCanvas_index = rtCanvas_index

    def get_leaf_nodes(self):
        self.leaf_set = set()
        for leaf_node in self.root.leaf_nodes():
            self.leaf_set.add(leaf_node.taxon.label)

    def set_block(self, block):
        self.block = block
        self.block_size = block.get_size()

    def set_pairwise_block(self, pairwise_block):
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
        self.topL = Point(x, y)
        self.width = self.tc_canvas.scale * DEFAULT_AD_WIDTH
        self.height = self.tc_canvas.scale * DEFAULT_AD_HEIGHT
        self.botR = Point(x + self.width, y + self.height)

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
                if type(block.belong_subtree) is list:
                    str = ""
                    for subtree in block.belong_subtree:
                        str += subtree.label + "&"
                    print(' ' * (level * 4) + str[:-1] + ':' )
                else:
                    print(' ' * (level * 4) + block.belong_subtree.label + ':')
            elif block.type == INDIVIDUAL_BLOCK or block.type == INDIVIDUAL_LEAF_BLOCK:
                if block.belong_subtree:
                    if type(block.belong_subtree) is list:
                        for subtree in block.belong_subtree:
                            str += subtree.label + "&"

                        print(' ' * (level * 4) + str[:-1] + ': INDIVIDUAL LEAF')
                    else:
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

    def get_all_subtree_block_list(self, subtree_block_list):
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

    def ad_to_string(self, canvas, node=None, differentiate_inexact_match=True, is_rt=False):
        if node == None:
            node = self.root
            newick_str = ""

        if node.type == LEAF:
            return self.get_node_type(canvas, node, differentiate_inexact_match=differentiate_inexact_match)
        else:
            newick_str = ""
            sorted_children = sorted(node.children, key=lambda x: x.type_prior, reverse=False)
            for child in sorted_children:
                nested_tree_string = ""
                if child.nested_tree:
                    nested_tree_string = child.nested_tree.ad_to_string(canvas=canvas,
                                                                        differentiate_inexact_match=differentiate_inexact_match)
                    newick_str += self.ad_to_string(canvas=canvas, node=child,
                                                    differentiate_inexact_match=differentiate_inexact_match) + '[' + nested_tree_string + '],'
                else:
                    newick_str += self.ad_to_string(canvas=canvas, node=child,
                                                    differentiate_inexact_match=differentiate_inexact_match) + ','

            if node.type == ROOT:
                return '(' + newick_str[:-1] + ')' + self.get_node_type(canvas, node,
                                                                        differentiate_inexact_match=differentiate_inexact_match) + ';'
            else:
                return '(' + newick_str[:-1] + ')' + self.get_node_type(canvas, node,
                                                                        differentiate_inexact_match=differentiate_inexact_match)

    def get_node_type(self, canvas, node, differentiate_inexact_match=True, is_rt=False):
        if node.type == LEAF:
            block = node.node_or_block
            if block.type == SUBTREE_BLOCK:
                subtree_label = ""
                if block.is_subtree_duplicate:
                    for index, subtree in enumerate(block.belong_subtree):
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

    def check_in_range(self, x, y):
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
        self.branch_head = Point(x - (AD_BRANCH_LENGTH * scale), y)
        self.branch_tail = Point(x, y)

    def construct_elided_branch(self, x, y, scale):
        self.branch_head = Point(x - (ELIDED_BRANCH_LENGTH * scale), y)
        self.branch_tail = Point(x, y)

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
        self.agree_rt = False

    def add_ad(self, ad_tree):
        self.ad_tree_list.append(ad_tree)
        self.tree_count += 1

    def generate_block_list(self):
        self.sample_ad_tree.get_all_subtree_block_list(self.block_list)

    def set_block_inexact_match(self, subtree_label):
        for block in self.block_list:
            if block.is_subtree_duplicate:
                for index, subtree in enumerate(block.belong_subtree):
                    if subtree_label == subtree.label:
                        block.exact_match[index] = False
            else:
                if subtree_label == block.belong_subtree.label:
                    block.exact_match = False


def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))

def rgb_to_hex(rgb_color):
    return '#{:02x}{:02x}{:02x}'.format(*rgb_color)

def lighten_color(hex_color, amount=10):
    rgb_color = hex_to_rgb(hex_color)
    lightened_rgb = tuple(min(255, c + amount) for c in rgb_color)
    return rgb_to_hex(lightened_rgb)
