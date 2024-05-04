from ipycanvas import Canvas,MultiCanvas,hold_canvas
from .myCanvas import *
from .myUtils import *
import math
import copy
import time

# Canvas for Reference Tree
class pairwiseCanvas(MyCanvas):
    # 0-5 = rt tree layer
    # 6 = empty layer
    # 7-13 = tc_tree RIGHT layer + escape taxa layer
    # 14-20 = tc_tree LEFT layer + escape taxa layer
    # 14 = HOVER_NODE_LAYER
    RT_LAYER = 5   # layer to draw reference tree
    TC_LEFT_LAYER = 20
    TC_RIGHT_LAYER = 13  # layer to draw tc_tree
    EMPTY_LAYER = 6
    HOVER_NODE_LAYER = -1  # layer to draw hover block
    SUBTREE_COMPARE_LAYER =  -2 # layer to compare subtree with
    RIGHT_ESCAPE_TAXA_LAYER = 12
    LEFT_ESCAPE_TAXA_LAYER = 19
    output = []

    def __init__(self,adPy, width, height, tc_tree=None):
        super().__init__(22 , width=width, height=height)
        # layer = subtree_count

        self.adPy = adPy
        self.rt = adPy.rt
        # tc_tree = copy.deepcopy(adPy.rt)
        self.tc_tree = tc_tree
        self.max_level = None

        self.left_escape_taxa_list = {}  # {'node':[subtree.color]}
        self.right_escape_taxa_list = {}
        self.left_escape_taxa_subtree = []
        self.right_escape_taxa_subtree = []

        self.rt_layer_block_list = {}
        self.tc_left_layer_block_list = {}
        self.tc_right_layer_block_list = {}

        self.rt_sorted_layer_list = []
        self.tc_left_sorted_layer_list = []
        self.tc_right_sorted_layer_list = []

        self.rt_layer_occupied = [0, 0, 0, 0, 0]
        self.tc_left_layer_occupied = [0, 0, 0, 0, 0]
        self.tc_right_layer_occupied = [0, 0, 0, 0, 0]

        self.subtree_label_used = self.adPy.subtree_label_used
        self.subtree_color_used = self.adPy.subtree_color_used

        self.tc_left_subtree_list = []
        self.tc_right_subtree_list = []

        self.left_tree_width = 0
        self.right_tree_width = self.width
        self.left_tree_height = 0
        self.right_tree_height = 0

        self.node_hover = None

        self.last_draw_time = time.time()
        self.last_mouse_position = None

        self.left_missing_taxa_list = []
        self.right_missing_taxa_list = []

        self.rt_section_list = []
        self.left_tree_section_list = []
        self.right_tree_section_list = []
        self.inner_section = True
        self.last_draw_time = time.time()
        self.last_mouse_position = None
        self.setup_section_list([self.rt_section_list,self.left_tree_section_list,self.right_tree_section_list])

        if self.rt.level > 0:
            self.max_level = self.rt.level
        else:
            self.max_level = self.get_max_level(self.rt)

        self.node_similarity_checked = False
        self.compare_between_tc = False

        self.max_level += 5
        self.default_value(align=LEFT)
        with hold_canvas(self):
            self.draw_tree(self[self.RT_LAYER], type=RT,align=LEFT, node=self.rt.seed_node, level=0)

        self.on_mouse_down(self.mouse_clicked)

        self.on_mouse_move(self.mouse_hover)

    def setup_section_list(self,tree_section_list):
        for section_list in tree_section_list:
            section_cnt = math.ceil(self.height/ MIN_SECTION_HEIGHT)

            for i in range(section_cnt):
                section_list.append([])
                if self.inner_section:
                    for j in range(3):
                        section_list[i].append([])

    def default_value(self,align):
        if align == LEFT:
            self.left_tree_width = 0
            self.left_tree_height = 0
            self.x = 30
            self.left_escape_taxa_list = {}  # {'node':[subtree.color]}
            self.left_escape_taxa_subtree = []
            self.tc_left_layer_block_list = {}
            self.tc_left_sorted_layer_list = []
            self.tc_left_layer_occupied = [0, 0, 0, 0, 0]
            self.tc_left_subtree_list = []
            self.left_tree_section_list = []
            self.setup_section_list([self.left_tree_section_list])
        elif align == RIGHT:
            self.x = self.width - 200
            self.right_escape_taxa_list = {}
            self.right_escape_taxa_subtree = []
            self.tc_right_layer_block_list = {}
            self.tc_right_sorted_layer_list = []
            self.tc_right_layer_occupied = [0, 0, 0, 0, 0]
            self.tc_right_subtree_list = []
            self.right_tree_width = self.width
            self.right_tree_height = 0
            self.right_tree_section_list = []
            self.setup_section_list([self.right_tree_section_list])

        # Common value
        self.y = 30

    def get_max_level(self,tree):
        max_level = 0
        for node in tree.leaf_nodes():
            if node.level() > max_level:
                max_level = node.level()

        return max_level

    def draw_tree(self, draw_canvas, type,align, node, level=0):   # align = left or right
        # Default : Tree root as reference tree's root

        node.is_missing = False
        for child in node.child_node_iter():
            self.draw_tree(draw_canvas=draw_canvas,type=type, align=align, node=child, level=level + 1)

        # Set canvas style
        draw_canvas.fill_style = BLACK
        draw_canvas.font = LEAF_FONT


        if node.is_leaf():
            if self.compare_between_tc:
                node.pairwise_similarity = 1.0
            # Calculate node's x,y coordinates
            x = PAIRWISE_X_INTERVAL * self.max_level + PAIRWISE_X_INTERVAL

            if align == RIGHT:
                x = self.width - PAIRWISE_RIGHT_PADDING - x

            node.pos_in_pairwise = Point(x, self.y)  # Node position

            self.y += RT_Y_INTERVAL

            # Mouse click range of current leaf node
            if align == LEFT:
                range_topL_in_pairwise = Point(PAIRWISE_X_INTERVAL * level + PAIRWISE_X_INTERVAL, node.pos_in_pairwise.y - NODE_BLOCK_HEIGHT)
                range_botR_in_pairwise = Point(node.pos_in_pairwise.x + PAIRWISE_X_INTERVAL + len(str(
                    node.taxon.label)) * 11,
                                   range_topL_in_pairwise.y + NODE_BLOCK_HEIGHT + 3)
                node.mouse_range_in_pairwise = Block(range_topL_in_pairwise, range_botR_in_pairwise)
                node.pairwise_block = Block(range_topL_in_pairwise, range_botR_in_pairwise)

                max_x = self.width / 2 - 55

                if range_botR_in_pairwise.x > self.left_tree_width:
                    self.left_tree_width = max_x
                if range_botR_in_pairwise.y > self.left_tree_height:
                    self.left_tree_height = range_botR_in_pairwise.y

            elif align == RIGHT:
                range_topL_in_pairwise = Point(node.pos_in_pairwise.x - PAIRWISE_X_INTERVAL - (len(str(node.taxon.label))) * 12, node.pos_in_pairwise.y -
                                                                       NODE_BLOCK_HEIGHT)
                x = self.width - PAIRWISE_RIGHT_PADDING - (PAIRWISE_X_INTERVAL * level + PAIRWISE_X_INTERVAL)
                range_botR_in_pairwise = Point(x,range_topL_in_pairwise.y + NODE_BLOCK_HEIGHT + 3)
                node.mouse_range_in_pairwise = Block(range_topL_in_pairwise, range_botR_in_pairwise)
                node.pairwise_block = Block(range_topL_in_pairwise, range_botR_in_pairwise)

                if range_topL_in_pairwise.x < self.right_tree_width:
                    self.right_tree_width = range_topL_in_pairwise.x
                if range_botR_in_pairwise.y > self.right_tree_height:
                    self.right_tree_height = range_botR_in_pairwise.y

            self.draw_leaf_node(draw_canvas, node,align)

            # Testing
            # self.draw_dots(node.pos_in_pairwise.x + X_INTERVAL - 5 , node.pos_in_pairwise.y - NODE_BLOCK_HEIGH  T)
            # self.draw_dots(node.pos_in_pairwise.x + X_INTERVAL + len(str(node.taxon.label)) * 9, node.pos_in_pairwise.y + 2)
            # self.draw_rec(node.pos_in_pairwise.x + X_INTERVAL - 2,node.pos_in_pairwise.y - NODE_BLOCK_HEIGHT,X_INTERVAL + len(str(node.taxon.label)) * 9,NODE_BLOCK_HEIGHT+2)
            # self.draw_frame(node.pairwise_block.topL.x,node.pairwise_block.topL.y,node.pairwise_block.width,node.pairwise_block.height,RED)

        else:
            if self.compare_between_tc and not self.node_similarity_checked:
                if node.corr:
                    compare_node = node.corr.corr[self.tc_tree[1].id]
                    self.check_similarity(node,compare_node)
                else:
                    node.pairwise_similarity = 1.0


            child_nodes = node.child_nodes()
            # To determine the midpoint of vertical line connecting multiple children
            first_child = child_nodes[0] if child_nodes else None
            last_child = child_nodes[-1] if child_nodes else None

            x = PAIRWISE_X_INTERVAL * level + PAIRWISE_X_INTERVAL

            if align == RIGHT:
                x = self.width - PAIRWISE_RIGHT_PADDING - x

            tmp = first_child.pos_in_pairwise.y + last_child.pos_in_pairwise.y
            node.pos_in_pairwise = Point(x, tmp / 2)

            line_tail = PAIRWISE_X_INTERVAL * (level + 1) + PAIRWISE_X_INTERVAL
            # Mouse click range of current internal node
            if align == LEFT:
                range_topL_in_pairwise = Point(node.pos_in_pairwise.x + 3, first_child.pos_in_pairwise.y)

                range_botR_in_pairwise = Point(range_topL_in_pairwise.x + PAIRWISE_X_INTERVAL, last_child.pos_in_pairwise.y - PAIRWISE_X_INTERVAL)
                node.mouse_range_in_pairwise = Block(range_topL_in_pairwise, range_botR_in_pairwise)

            elif align == RIGHT:
                range_topL_in_pairwise = Point(node.pos_in_pairwise.x - 3 - PAIRWISE_X_INTERVAL, first_child.pos_in_pairwise.y)
                range_botR_in_pairwise = Point(node.pos_in_pairwise.x - 3, last_child.pos_in_pairwise.y - PAIRWISE_X_INTERVAL)
                node.mouse_range_in_pairwise = Block(range_topL_in_pairwise, range_botR_in_pairwise)
                line_tail = self.width - PAIRWISE_RIGHT_PADDING - line_tail

            # Drawing internal node's horizontal branch
            draw_canvas.begin_path()
            draw_canvas.move_to(node.pos_in_pairwise.x, node.pos_in_pairwise.y - 5)
            draw_canvas.line_to(line_tail, node.pos_in_pairwise.y - 5)
            draw_canvas.stroke()

            # Draw vertical branch
            draw_canvas.begin_path()
            draw_canvas.move_to(line_tail, first_child.pos_in_pairwise.y - 5)
            draw_canvas.line_to(line_tail, last_child.pos_in_pairwise.y - 5)
            draw_canvas.stroke()

            # Draw childs' horizontal line
            draw_canvas.begin_path()
            draw_canvas.move_to(line_tail,first_child.pos_in_pairwise.y - 5)
            draw_canvas.line_to(first_child.pos_in_pairwise.x , first_child.pos_in_pairwise.y - 5)
            draw_canvas.stroke()

            draw_canvas.begin_path()
            draw_canvas.move_to(line_tail, last_child.pos_in_pairwise.y - 5)
            draw_canvas.line_to(last_child.pos_in_pairwise.x,last_child.pos_in_pairwise.y - 5)
            draw_canvas.stroke()

            # Testing
            # self.draw_dots(node.pos_in_pairwise.x + X_INTERVAL - 4, node.pos_in_pairwise.y - X_INTERVAL)
            # self.draw_dots(node.pos_in_pairwise.x + X_INTERVAL + 4, node.pos_in_pairwise.y )
            # self.draw_rec(node.range_topL_in_pairwise.x,node.range_topL_in_pairwise.y,8,node.range_botR_in_pairwise.y - node.range_topL_in_pairwise.y,BEIGE)
            # self.draw_frame(node.pairwise_block.topL.x, node.pairwise_block.topL.y, node.pairwise_block.width,node.pairwise_block.height, RED)


        # To calculate and record subtree's block size
        if node.parent_node:
            if hasattr(node.parent_node, 'pairwise_block'):
                if node.pairwise_block.topL.y < node.parent_node.pairwise_block.topL.y:
                    node.parent_node.pairwise_block.topL = node.pairwise_block.topL

                if node.pairwise_block.botR.y > node.parent_node.pairwise_block.botR.y:
                    node.parent_node.pairwise_block.botR = node.pairwise_block.botR
            else:
                node.parent_node.pairwise_block = Block(node.pairwise_block.topL, node.pairwise_block.botR)

        node.selected = False
        if align == LEFT:
            if type == RT:
                self.insert_node_section_list(self.rt_section_list,node)
            else:
                self.insert_node_section_list(self.left_tree_section_list, node)
        else:
            self.insert_node_section_list(self.right_tree_section_list, node)


        # Testing
        # self.draw_dots(node.pos_in_pairwise.x,node.pos_in_pairwise.y - X_INTERVAL)
        # self.draw_dots(node.pos_in_pairwise.x + X_INTERVAL, node.pos_in_pairwise.y)
        # self.draw_dots(node.pos_in_pairwise.x + X_INTER0L.y,8,node.range_botR_in_pairwise.y - node.range_topL_in_pairwise.y,BEIGE)

    def insert_node_section_list(self, section_list,node):
        top_section_index = math.floor(node.mouse_range_in_pairwise.topL.y / MIN_SECTION_HEIGHT)
        bottom_section_index = math.floor(node.mouse_range_in_pairwise.botR.y / MIN_SECTION_HEIGHT)

        if not self.inner_section:
            if top_section_index == bottom_section_index:
                section_list[top_section_index].append(node)
            else:
                for i in range(top_section_index, bottom_section_index + 1):
                    section_list[i].append(node)

            return

        if top_section_index == bottom_section_index:
            section_midpoint = top_section_index * MIN_SECTION_HEIGHT + MIN_SECTION_HEIGHT / 2
            if node.mouse_range_in_pairwise.botR.y <= section_midpoint:
                section_list[top_section_index][1].append(node)
            elif node.mouse_range_in_pairwise.topL.y >= section_midpoint:
                section_list[top_section_index][2].append(node)
            else:
                section_list[top_section_index][0].append(node)
        else:
            for i in range(top_section_index, bottom_section_index + 1):
                section_midpoint = i * MIN_SECTION_HEIGHT + MIN_SECTION_HEIGHT / 2
                if i == 0:
                    if node.mouse_range_in_pairwise.topL.y <= section_midpoint:
                        section_list[i][1].append(node)
                    else:
                        section_list[i][2].append(node)
                elif i == bottom_section_index:
                    if node.mouse_range_in_pairwise.topL.y <= section_midpoint:
                        section_list[i][1].append(node)
                    else:
                        section_list[i][2].append(node)
                else:
                    section_list[i][0].append(node)
    def draw_leaf_node(self,canvas,node,align):
        canvas.text_align = align
        canvas.font = LEAF_FONT
        canvas.fill_style = BLACK
        if align == LEFT:
            canvas.fill_text(node.taxon.label, node.pos_in_pairwise.x + PAIRWISE_X_INTERVAL, node.pos_in_pairwise.y)
        elif align == RIGHT:
            canvas.fill_text(node.taxon.label, node.pos_in_pairwise.x - PAIRWISE_X_INTERVAL, node.pos_in_pairwise.y)
    def filter_node_selected(self,node, x , y):
        if node.mouse_range_in_pairwise.topL.x <= x <= node.mouse_range_in_pairwise.botR.x and node.mouse_range_in_pairwise.topL.y <= y <= node.mouse_range_in_pairwise.botR.y:
            return True
        return False

    def filter_node_from_section_list(self,section_list,x,y):

        node_selected = None
        inner_section_check = [1,1,1]

        section_index = math.floor(y / MIN_SECTION_HEIGHT)


        inner_section_check[0] = 0
        for node in section_list[section_index][0]:
            if self.filter_node_selected(node, x, y):
                return node

        section_midpoint = section_index * MIN_SECTION_HEIGHT + MIN_SECTION_HEIGHT / 2
        if y <= section_midpoint:
            subsection_list = section_list[section_index][1]
            inner_section_check[1] = 0
        else:
            subsection_list = section_list[section_index][2]
            inner_section_check[2] = 0

        for node in subsection_list:
            if self.filter_node_selected(node, x, y):
                return node

        subsection_index = inner_section_check.index(1)
        subsection_list = section_list[section_index][subsection_index]
        for node in subsection_list:
            if self.filter_node_selected(node, x, y):
                return node

        return node_selected

    def mouse_clicked(self, x, y):
        if self.compare_between_tc:
            return

        # node_selected = self.rt.find_node(lambda node: self.filter_node_selected(node, x, y))
        if self.node_hover and self.filter_node_selected(self.node_hover ,x, y):
            node_selected = self.node_hover
        else:
            node_selected = self.filter_node_from_section_list(self.rt_section_list,x, y)

        if not node_selected:
            return

        # Testing
        # self.draw_frame(node_selected.pairwise_block.topL.x,node_selected.pairwise_block.topL.y,10,10,BLACK)

        self.adPy.select_subtree_from_tree(node_selected)

    def mouse_hover(self, x, y):
        type = None
        align = None
        node_selected = None
        subtree_root = False

        current_time = time.time()
        if current_time - self.last_draw_time > 0.2:
            self.last_draw_time = current_time

            # Filter node
            if x <= self.left_tree_width:
                align = LEFT

                if self.compare_between_tc:
                    type = TC
                    node_selected = self.filter_node_from_section_list(self.left_tree_section_list, x, y)


                    if node_selected in self.left_escape_taxa_list:
                        subtree_root = True
                    for subtree in self.adPy.subtree_list:
                        if node_selected == subtree.root.corr[self.tc_tree[0].id]:
                            subtree_root = True

                else:
                    type = RT
                    node_selected = self.filter_node_from_section_list(self.rt_section_list, x, y)

            else:
                align = RIGHT
                node_selected = self.filter_node_from_section_list(self.right_tree_section_list, x, y)
                type = TC

                if node_selected in self.right_escape_taxa_list:
                    subtree_root = True

                if self.compare_between_tc:
                    tc_tree = self.tc_tree[1]
                else:
                    tc_tree = self.tc_tree

                for subtree in self.adPy.subtree_list:
                    if node_selected == subtree.root.corr[tc_tree.id]:
                        subtree_root = True
        else:
            return

        if not node_selected:
            self.node_hover = node_selected
            self[self.HOVER_NODE_LAYER].clear()
            self[self.HOVER_NODE_LAYER].flush()

        if node_selected != self.node_hover:
            self[self.HOVER_NODE_LAYER].clear()
            self[self.HOVER_NODE_LAYER].flush()

        if node_selected and node_selected != self.node_hover:
            self.node_hover = node_selected
            self.draw_hover_block(type=type,align=align,node=self.node_hover)
            if not node_selected.is_missing:
                self.draw_node_details(type,node_selected,subtree_root)

    def draw_node_details(self,type,node,subtree_root):
        # support , length , exact match
        label_length = 0

        if node.is_leaf():
            support = "Support : " + '-'
        else:
            support = "Support : "  + node.label


        length = node.edge.length if node.edge.length != None else 0
        length = "Length : " + str(length)

        if len(length) * 8 > label_length:
            label_length = len(length) * 8

        if type == TC and len(self.adPy.tc) > 0 and subtree_root:
            exact_match =  str(node.exact_match_percentage * 100)
        else:
            if type == TC:
                exact_match = str(node.corr_similarity * 100) + "%"
            else:
                exact_match = format(node.exact_match / len(self.adPy.tc) * 100, ".2f")

        exact_match = "Exact Match :  " + exact_match + "%"

        if len(exact_match) * 8 > label_length:
            label_length = len(exact_match) * 8

        similarity = ""
        if type == TC:
            if self.compare_between_tc:
                similarity = "Similarity :  " + str(node.pairwise_similarity * 100) + "%"
            else:
                similarity = "Similarity :  " + str(node.corr_similarity * 100) + "%"

        if len(similarity) * 10 > label_length:
            label_length = len(similarity) * 10

        label_x = self.left_tree_width + ((self.right_tree_width - self.left_tree_width - 50) / 2)
        label_y = node.pos_in_pairwise.y
        if node.pos_in_pairwise.y + 40 > self.height:
            label_y = self.height - 40

        if subtree_root:
            self.draw_rec(label_x-10 , label_y - 15, label_length + 10,65,BLACK,self.HOVER_NODE_LAYER)
        else:
            self.draw_rec(label_x-10, label_y - 15, label_length + 10, 50,BLACK,self.HOVER_NODE_LAYER)

        self[self.HOVER_NODE_LAYER].begin_path()
        self[self.HOVER_NODE_LAYER].fill_style = BLANK
        self[self.HOVER_NODE_LAYER].font = LEAF_FONT

        self[self.HOVER_NODE_LAYER].fill_text(support, label_x  , label_y)
        self[self.HOVER_NODE_LAYER].fill_text(length, label_x , label_y+ 15)
        if type == RT:
            self[self.HOVER_NODE_LAYER].fill_text(exact_match, label_x , label_y + 30)
        elif type == TC:
            if subtree_root:
                self[self.HOVER_NODE_LAYER].fill_text(exact_match, label_x,
                                                      label_y + 45)
            self[self.HOVER_NODE_LAYER].fill_text(similarity, label_x,label_y + 30)

    def draw_hover_block(self,type,align,node,corr=False):
        if align == LEFT:
            if not node.is_leaf():
                self.draw_rec(node.pos_in_pairwise.x + 2, node.pos_in_pairwise.y - RT_X_INTERVAL + 2 , RT_X_INTERVAL,
                              RT_X_INTERVAL - 2, BLACK, layer_index=self.HOVER_NODE_LAYER)
            else:
                if node.is_missing:
                    radius = 4
                    self[self.HOVER_NODE_LAYER].fill_style = BLACK
                    self[self.HOVER_NODE_LAYER].fill_arc(self.left_tree.missing_taxa_block.botR.x -
                                                         RT_X_INTERVAL - 2,
                                                         node.pos_in_pairwise.y + RT_X_INTERVAL,
                                                         radius, 0, 360)
                else:
                    width = node.pos_in_pairwise.x - node.mouse_range_in_pairwise.topL.x
                    self.draw_rec(node.mouse_range_in_pairwise.topL.x, node.pos_in_pairwise.y - RT_X_INTERVAL,width,
                                  RT_X_INTERVAL - 2, BLACK, layer_index = self.HOVER_NODE_LAYER)
            if corr:
                return
            else:
                if self.compare_between_tc:
                    node = node.corr.corr[self.tc_tree[1].id]
                    self.draw_hover_block(type=TC,align=RIGHT,node=node,corr=True)
                else:
                    self.draw_hover_block(type=TC, align=RIGHT, node=node.corr[self.tc_tree.id], corr=True)
        elif align == RIGHT:
            if not node.is_leaf():
                self.draw_rec(node.pos_in_pairwise.x - RT_X_INTERVAL - 2, node.pos_in_pairwise.y - RT_X_INTERVAL + 2 ,
                              RT_X_INTERVAL,
                              RT_X_INTERVAL - 2, BLACK, layer_index=self.HOVER_NODE_LAYER)
            else:
                if node.is_missing:
                    radius = 4
                    self[self.HOVER_NODE_LAYER].fill_style = BLACK
                    self[self.HOVER_NODE_LAYER].fill_arc(node.pos_in_pairwise.x - RT_X_INTERVAL,
                                                         node.pos_in_pairwise.y + RT_X_INTERVAL,
                                                              radius, 0, 360)
                else:
                    width = node.mouse_range_in_pairwise.botR.x - node.pos_in_pairwise.x
                    self.draw_rec(node.pos_in_pairwise.x, node.pos_in_pairwise.y - RT_X_INTERVAL, width,
                                  RT_X_INTERVAL - 2, BLACK, layer_index = self.HOVER_NODE_LAYER)
            if corr:
                return
            else:
                if self.compare_between_tc:
                    node = node.corr.corr[self.tc_tree[0].id]

                    self.draw_hover_block(type=TC,align=LEFT,node=node,corr=True)
                else:
                    self.draw_hover_block(type=RT,align=LEFT,node=node.corr,corr=True)

        self[self.HOVER_NODE_LAYER].flush()


    def write_block_label(self,subtree,layer_index,align=LEFT):
        node_amt = len(subtree.root.leaf_nodes())
        label_str = f"{subtree.label} : {node_amt}"
        subtree.label_width = (len(label_str) * 6)

        if align == LEFT:
            label_x = subtree.pairwise_block.botR.x - subtree.label_width
        elif align == RIGHT:
            label_x = subtree.pairwise_block.topL.x + 5

        label_y = subtree.pairwise_block.topL.y + 12
        self[layer_index].fill_style = BLACK
        self[layer_index].fill_text(label_str, label_x, label_y)

        # # Testing
        # self[layer_index].fill_text("layer: " + str(layer_index), subtree.pairwise_block.botR.x + 5, subtree.pairwise_block.topL.y + 5)

    def generate_subtree_block(self,node_selected,align=LEFT):
        if align == LEFT:
            if node_selected.is_leaf():
                rec_x = node_selected.pairwise_block.topL.x + RT_X_INTERVAL - 5
            else:
                rec_x = node_selected.mouse_range_in_pairwise.topL.x

            rec_y = node_selected.pairwise_block.topL.y
            rec_width = self.left_tree_width + RT_X_INTERVAL - rec_x
            rec_height = node_selected.pairwise_block.botR.y - rec_y

        elif align == RIGHT:
            rec_x = self.left_tree_width + TREE_PADDING
            rec_y = node_selected.pairwise_block.topL.y

            if node_selected.is_missing:
                rec_x = node_selected.pairwise_block.topL.x - 5
                rec_width = node_selected.pairwise_block.botR.x - rec_x - RT_X_INTERVAL
            elif node_selected.is_leaf():
                rec_width = node_selected.mouse_range_in_pairwise.botR.x - rec_x - RT_X_INTERVAL
            else:
                rec_width = node_selected.mouse_range_in_pairwise.botR.x - rec_x

            rec_height = node_selected.pairwise_block.botR.y - rec_y

        new_block = Block(Point(rec_x, rec_y), Point(rec_x + rec_width, rec_y + rec_height))

        return new_block

    def get_node_from_tree(self,tree,taxa_name):

        if taxa_name in tree.missing:
            for node in tree.missing_node_list:
                if node.label == taxa_name:
                    return node
        else:
            node = tree.find_node_with_taxon_label(taxa_name)

        return node

    def remove_escape_taxa(self,subtree_label,escape_taxa_list):
        escape_taxa_delete = []
        for escape_taxa,subtree_list in escape_taxa_list.items():

            for subtree in subtree_list:
                if subtree_label == subtree.label:
                    subtree_list.remove(subtree)
                    if len(subtree_list) == 0:
                        escape_taxa_delete.append(escape_taxa)
        for taxa in escape_taxa_delete:
            del escape_taxa_list[taxa]

        return


    def draw_escape_taxa(self,align):
        if align == LEFT:
            tree = self.left_tree
            layer_index = self.LEFT_ESCAPE_TAXA_LAYER
            escape_taxa_list = self.left_escape_taxa_list
        else:
            tree = self.right_tree
            layer_index = self.RIGHT_ESCAPE_TAXA_LAYER
            escape_taxa_list = self.right_escape_taxa_list

        self[layer_index].clear()
        for escape_taxa,subtree_list in escape_taxa_list.items():
            node = self.get_node_from_tree(tree,escape_taxa)
            block = node.escape_taxa_block
            if len(subtree_list) == 1:
                subtree = subtree_list[0]
                self.draw_rec(block.topL.x, block.topL.y, block.width, block.height,subtree.color,layer_index=layer_index)
            else:
                segment_width = block.width / len(subtree_list)
                for index,subtree in enumerate(subtree_list):
                    self.draw_rec(block.topL.x + segment_width * index, block.topL.y, segment_width, block.height,subtree.color,layer_index=layer_index)

    def record_escape_taxa(self,node_selected,subtree,align):
        if align == LEFT:
            taxa_list = self.left_escape_taxa_list
            taxa_subtree_list = self.left_escape_taxa_subtree
        else:
            taxa_list = self.right_escape_taxa_list
            taxa_subtree_list = self.right_escape_taxa_subtree

        node_selected.escape_taxa_block = self.generate_subtree_block(node_selected,align=align)
        taxa_subtree_list.append(subtree.label)
        if node_selected.taxon.label not in taxa_list:
            taxa_list[node_selected.taxon.label] = []

        taxa_list[node_selected.taxon.label].append(subtree)

    def record_missing_taxa(self,node_selected,subtree,align):
        if align == LEFT:
            taxa_list = self.left_missing_taxa_list
        else:
            taxa_list = self.right_missing_taxa_list

        if hasattr(node_selected,'subtree'):
            if type(node_selected.subtree) is list:
                node_selected.subtree.append(subtree)
            else:
                subtree_list = [node_selected.subtree,subtree]
                node_selected.subtree = subtree_list
        else:
            node_selected.subtree = subtree

        node_selected.selected = True
        taxa_list.append(node_selected)

    def draw_missing_taxa(self,align):
        if align == LEFT:
            taxa_list = self.left_missing_taxa_list
        else:
            taxa_list = self.right_missing_taxa_list

        for taxa in taxa_list:
            self.draw_subtree_block(taxa,select_tree=TC,new_subtree=taxa.subtree,align=align)

    def draw_subtree_block(self,node_selected,select_tree=None,new_subtree=None,subtree_from_rt=None,align=None):
        new_block = self.generate_subtree_block(node_selected, align=align)
        new_subtree.set_pairwise_block(new_block)

        if select_tree == RT:
            first_layer = 0
            last_layer = 5
        elif select_tree == TC:
            if align == RIGHT:
                first_layer = 7
                last_layer = 12
            else:
                first_layer = 14
                last_layer = 19

        if select_tree == RT:
            self.draw_subtree_list = sorted(self.adPy.subtree_list, key=lambda x: x.block_size,reverse=True)
        else:
            if align == LEFT:
                self.draw_subtree_list = sorted(self.tc_left_subtree_list, key=lambda x: x.block_size,reverse=True)
            else:
                self.draw_subtree_list = sorted(self.tc_right_subtree_list, key=lambda x: x.block_size, reverse=True)

        subtree_cnt = len(self.draw_subtree_list)
        tmp_sorted_layer_list = []
        for i in range(first_layer,first_layer + subtree_cnt):
            if select_tree == RT:
                subtree = self.draw_subtree_list[i]
            else:
                if align == RIGHT:
                    subtree = self.draw_subtree_list[i - 7]
                else:
                    subtree = self.draw_subtree_list[i - 14]

            subtree_block = subtree.pairwise_block
            tmp_sorted_layer_list.append(subtree.label)

            self[i].clear()
            if subtree.root.is_missing:
                subtree_block = subtree.root.pairwise_block


            self.draw_rec(subtree_block.topL.x, subtree_block.topL.y, subtree_block.width, subtree_block.height,subtree.color,
                      layer_index=i)

            if subtree.root.is_missing == False:
                self.write_block_label(subtree, i,align=align)

            self[i].label = subtree.label
            self[i].flush()

        if select_tree == RT:
            self.rt_sorted_layer_list = tmp_sorted_layer_list
        elif select_tree == TC:
            if align == LEFT:
                self.tc_left_sorted_layer_list = tmp_sorted_layer_list
            else:
                self.tc_right_sorted_layer_list = tmp_sorted_layer_list

    def remove_subtree_block(self,subtree):
        # Clear corresponding rt layer
        try:
            tc_clear_layer = None
            clear_layer = self.rt_sorted_layer_list.index(subtree.label)

            self[clear_layer].clear()

            if self.tc_tree:
                if self.compare_between_tc:
                    if subtree.label in self.tc_left_sorted_layer_list:
                        tc_clear_layer = self.tc_left_sorted_layer_list.index(subtree.label)
                        self[tc_clear_layer + 14].clear()

                if subtree.label in self.tc_right_sorted_layer_list:
                    tc_clear_layer = self.tc_right_sorted_layer_list.index(subtree.label)
                    self[tc_clear_layer + 7].clear()

            if self.compare_between_tc:
                left_subtree = self.get_subtree(self.tc_left_subtree_list, subtree.label)
                self.tc_left_subtree_list.remove(left_subtree)
                self.remove_escape_taxa(subtree.label,self.left_escape_taxa_list)

            right_subtree = self.get_subtree(self.tc_right_subtree_list, subtree.label)
            if right_subtree :
                self.tc_right_subtree_list.remove(right_subtree)
            self.remove_escape_taxa(subtree.label,self.right_escape_taxa_list)

            if self.compare_between_tc:
                self.draw_escape_taxa(LEFT)

            self.draw_escape_taxa(RIGHT)

        except Exception as err:
            self[-1].fill_text("Loading", self.width - 200, 20)


    def setup_tc_subtree_list(self,tree,align):
        for index,subtree in enumerate(self.adPy.subtree_list):
            self.draw_tc_subtree_block(subtree,align)

        self.draw_escape_taxa(align)
    def check_block_exact_match(self,tc_subtree,rt_subtree):
        if tc_subtree.leaf_set.issubset(rt_subtree.leaf_set):
            return True
        else:
            return False


    def draw_tc_subtree_block(self,subtree,align):
        if align == LEFT:
            subtree_list = self.tc_left_subtree_list
            tree = self.tc_tree[0]
        else:
            if self.compare_between_tc:
                tree = self.tc_tree[1]
            else:
                tree = self.tc_tree
            subtree_list = self.tc_right_subtree_list

        corresponding_subtree = subtree.root.corr[tree.id]

        # If taxa in corresponding subtree is all missing taxa in target tree

        new_subtree = Subtree(label=subtree.label, belong_tree=tree, root=corresponding_subtree,
                                   color=subtree.color)
        subtree_list.append(new_subtree)
        subtree.corresponding_tc_subtree = new_subtree

        if corresponding_subtree == 0:
            for leaf_node in subtree.root.leaf_nodes():
                tc_corr = leaf_node.corr[tree.id]

                self.record_escape_taxa(tc_corr, new_subtree, align)
            return


        new_subtree.get_leaf_nodes()
        subtree.get_leaf_nodes()

        exact_match = False
        if new_subtree.leaf_set == subtree.leaf_set:
            exact_match = True

        new_subtree.root.exact_match_percentage = len(subtree.leaf_set.intersection(new_subtree.leaf_set)) / len(
            subtree.leaf_set)

        if not exact_match:
            tmp_result = self.check_block_exact_match(new_subtree,subtree)
            if not tmp_result:
                # All leaf nodes as escape taxa
                for leaf_node in subtree.root.leaf_nodes():
                    tc_corr = leaf_node.corr[tree.id]
                    tc_corr.exact_match_percentage = new_subtree.root.exact_match_percentage
                    self.record_escape_taxa(tc_corr,new_subtree,align)
                return
            else:
                # Draw leaf nodes which not in new_subtree
                escape_taxa = subtree.leaf_set.difference(new_subtree.leaf_set)
                for leaf_node in subtree.root.leaf_nodes():
                    if leaf_node.taxon.label not in escape_taxa:
                        continue

                    tc_corr = leaf_node.corr[tree.id]
                    tc_corr.exact_match_percentage = new_subtree.root.exact_match_percentage
                    self.record_escape_taxa(tc_corr, new_subtree,align)

        self.check_duplicate_subtree(corresponding_subtree, subtree, new_subtree)
        for leaf_node in corresponding_subtree.leaf_nodes():
            leaf_node.selected = True


        self.draw_subtree_block(corresponding_subtree, select_tree=TC, new_subtree=new_subtree,
                                subtree_from_rt=subtree,align=align)

    def check_duplicate_subtree(self,corresponding_subtree,subtree,new_subtree):
        if not hasattr(corresponding_subtree, 'subtree') or not corresponding_subtree.subtree:
            corresponding_subtree.subtree = new_subtree

        elif hasattr(corresponding_subtree, 'subtree'):
            if type(corresponding_subtree.subtree) is not list and corresponding_subtree.subtree != subtree:
                corresponding_subtree.duplicate_subtree = True
                subtree_list = []
                subtree_list.append(corresponding_subtree.subtree)
                subtree_list.append(new_subtree)
                corresponding_subtree.subtree = subtree_list
            else:
                corresponding_subtree.duplicate_subtree = True
                corresponding_subtree.subtree.append(new_subtree)

        corresponding_subtree.selected = True

    def compare_tc_tree(self,compare_tree=None,compare_between_tc=False):
        self[self.TC_LEFT_LAYER].clear()
        self[self.TC_RIGHT_LAYER].clear()
        self.tc_tree = compare_tree
        self.compare_between_tc = compare_between_tc
        self.last_draw_time = time.time()
        self.reset_subtree_canvas()

        if compare_between_tc:
            self[self.EMPTY_LAYER].fill_style = BLANK
            self[self.EMPTY_LAYER].fill_rect(0,0,self.width,self.height)

            level_tmp = self.get_max_level(self.tc_tree[0])
            if level_tmp > self.max_level:
                self.max_level = level_tmp

            self.left_tree = self.tc_tree[0]
            self.default_value(align=LEFT)
            with (hold_canvas(self)):
                self.draw_tree(self[self.TC_LEFT_LAYER], type=TC,align=LEFT, node=self.tc_tree[0].seed_node, level=0)
                self.node_similarity_checked = True

                # Draw tree missing taxa
                if len(self.left_tree.missing) > 0:
                    self.generate_missing_taxa_block_in_tree(self.left_tree, LEFT)
                    self.draw_missing_taxa_block(align=LEFT)

            # Draw tc tree (align right)
            level_tmp = self.get_max_level(self.tc_tree[1])
            if level_tmp > self.max_level:
                self.max_level = level_tmp

            self.right_tree = self.tc_tree[1]
            self.default_value(align=RIGHT)
            with hold_canvas(self):
                self.draw_tree(self[self.TC_RIGHT_LAYER], type=TC,align=RIGHT, node=self.tc_tree[1].seed_node, level=0)

                # Draw tree missing taxa
                if len(self.right_tree.missing) > 0:
                    self.generate_missing_taxa_block_in_tree(self.right_tree, RIGHT)
                    self.draw_missing_taxa_block(align=RIGHT)

            if len(self.adPy.subtree_list) > 0:
                self.setup_tc_subtree_list(self.tc_tree[0],LEFT)
                self.setup_tc_subtree_list(self.tc_tree[1],RIGHT)
        else:
            self[self.EMPTY_LAYER].clear()
            self.left_tree = self.rt
            self.right_tree = self.tc_tree

            level_tmp = self.get_max_level(self.tc_tree)
            if level_tmp > self.max_level:
                self.max_level = level_tmp

            self.default_value(align=RIGHT)

            with (hold_canvas(self)):
                self.draw_tree(self[self.TC_RIGHT_LAYER], type = TC,align=RIGHT, node=self.tc_tree.seed_node, level=0)

                if len(self.right_tree.missing) > 0:
                    self.generate_missing_taxa_block_in_tree(self.right_tree,RIGHT)
                    self.draw_missing_taxa_block(align=RIGHT)

                    # self.draw_frame(self.right_tree_width, y, width, height, BLACK, layer_index=self.TC_RIGHT_LAYER)

            if len(self.adPy.subtree_list) > 0:
                self.setup_tc_subtree_list(self.tc_tree,RIGHT)

    def generate_missing_taxa_block_in_tree(self,tree,align):
        if align == RIGHT:
            x = self.width / 2 + PAIRWISE_X_INTERVAL * 2
            end_x = x + MISSING_TAXA_BLOCK_WIDTH
            y = self.right_tree_height + 20
            layer_index = self.TC_RIGHT_LAYER
        else:
            x = self.x
            end_x = x + MISSING_TAXA_BLOCK_WIDTH
            y = self.left_tree_height + 20

            layer_index = self.TC_LEFT_LAYER

        height = (len(tree.missing) + 1) * RT_Y_INTERVAL
        width = end_x - x

        tree.missing_taxa_block = Block(topL=Point(x, y), botR=Point(x + width, y + height))

        self.draw_dotted_frame(x, y, MISSING_TAXA_BLOCK_WIDTH, height, color=BLACK, line_dash_index=1,
                               layer_index=layer_index)

    def draw_missing_taxa_block(self,align):
        if align == RIGHT:
            block = self.right_tree.missing_taxa_block
            tree = self.right_tree
            layer_index = self.TC_RIGHT_LAYER
            label_x = block.topL.x + PAIRWISE_X_INTERVAL * 4
        else:
            block = self.left_tree.missing_taxa_block
            tree = self.left_tree
            layer_index = self.TC_LEFT_LAYER
            label_x = block.botR.x - PAIRWISE_X_INTERVAL * 4

        label_y = block.topL.y + PAIRWISE_X_INTERVAL

        for node in tree.missing_node_list:
            if align == LEFT:
                node.pos_in_pairwise = Point(block.topL.x + PAIRWISE_X_INTERVAL, label_y + 6)
            else:
                node.pos_in_pairwise = Point(label_x - 7, label_y + 6)

            if align == LEFT:
                end_x = label_x + 6
            else:
                end_x = block.botR.x - PAIRWISE_X_INTERVAL

            range_botR_in_pairwise = Point(end_x, label_y + RT_Y_INTERVAL + 3)

            node.mouse_range_in_pairwise = Block(node.pos_in_pairwise, range_botR_in_pairwise)
            node.pairwise_block = Block(node.pos_in_pairwise, range_botR_in_pairwise)

            label_y += RT_Y_INTERVAL

            self[layer_index].fill_style = BLACK
            self[layer_index].font = LEAF_FONT

            if align == RIGHT:
                self[layer_index].text_align = LEFT
            else:
                self[layer_index].text_align = RIGHT

            self[layer_index].fill_text(node.label, label_x, label_y)

            if align == RIGHT:
                self.insert_node_section_list(self.right_tree_section_list, node)
            else:
                self.insert_node_section_list(self.left_tree_section_list, node)


    def reset_subtree_canvas(self):
        for i in range(0,22):
            if i == 5 or i == 13 or i == 20:
                continue
            self[i].clear()

    def check_similarity(self,node,compare_node):

        first_leaf_nodes = {leaf.taxon.label for leaf in node.leaf_nodes()}
        second_leaf_nodes = {leaf.taxon.label for leaf in compare_node.leaf_nodes()}

        intersection_nodes = first_leaf_nodes.intersection(second_leaf_nodes)

        node.pairwise_similarity = len(intersection_nodes)/len(first_leaf_nodes)
        compare_node.pairwise_similarity = len(intersection_nodes) / len(second_leaf_nodes)



