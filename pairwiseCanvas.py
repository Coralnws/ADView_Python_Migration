from ipycanvas import Canvas,MultiCanvas
import ADpy as ADpy
from myCanvas import *
import tcCanvas as tcCanvas
from Utils import *
import math
import copy

# Canvas for Reference Tree
class pairwiseCanvas(MyCanvas):
    # 0-5 = rt tree layer
    # 6-11 = tc_tree layer
    # 12 = SUBTREE_COMPARE_LAYER
    # 13 = HOVER_NODE_LAYER
    RT_LAYER = 5   # layer to draw reference tree
    TC_LAYER = 12  # layer to draw tc_tree
    HOVER_NODE_LAYER = -1  # layer to draw hover block
    SUBTREE_COMPARE_LAYER =  -2 # layer to compare subtree with
    ESCAPE_TAXA_LAYER = 11


    def __init__(self,adPy, width, height, tc_tree=None):
        super().__init__(15 , width=width, height=height)
        # layer = subtree_count

        self.adPy = adPy
        self.rt = adPy.rt
        # tc_tree = copy.deepcopy(adPy.rt)
        self.tc_tree = tc_tree
        self.max_level = None

        self.escape_taxa_list = {}  # {'node':[subtree.color]}
        self.escape_taxa_subtree = []

        self.rt_layer_block_list = {}
        self.tc_layer_block_list = {}

        self.rt_sorted_layer_list = []
        self.tc_sorted_layer_list = []

        self.rt_layer_occupied = [0, 0, 0, 0, 0]
        self.tc_layer_occupied = [0, 0, 0, 0, 0]

        self.subtree_label_used = self.adPy.subtree_label_used
        self.subtree_color_used = self.adPy.subtree_color_used

        self.tc_subtree_list = []

        self.left_tree_width = 0
        self.right_tree_width = self.width

        if self.rt.level > 0:
            self.max_level = self.rt.level
        else:
            self.max_level = self.get_max_level(self.rt)

        self.max_level += 1

        # Draw reference tree (align left)
        self.default_value(align=LEFT)
        self.draw_tree(self[self.RT_LAYER],align=LEFT,node=self.rt.seed_node,level=0)
        # self.setup_subtree_list()

        if tc_tree:
            # Draw tc tree (align right)
            self.default_value(align=RIGHT)
            self.draw_tree(self[self.TC_LAYER],align=RIGHT,node=self.tc_tree.seed_node, level=0)
            # self.setup_subtree_list()

        self.on_mouse_down(self.mouse_clicked)
        self.on_mouse_move(self.mouse_hover)

    def default_value(self,align):
        if align == LEFT:
            self.x = 30
        elif align == RIGHT:
            self.x = self.width - 200

        # Common value
        self.y = 30

    def get_max_level(self,tree):
        max_level = 0
        for node in tree.leaf_nodes():
            if node.level() > max_level:
                max_level = node.level()

        return max_level

    def draw_tree(self, draw_canvas, align, node, level=0):   # align = left or right
        # Default : Tree root as reference tree's root

        for child in node.child_node_iter():
            self.draw_tree(draw_canvas=draw_canvas, align=align, node=child, level=level + 1)

        # Set canvas style
        draw_canvas.fill_style = BLACK
        draw_canvas.font = LEAF_FONT

        if node.is_leaf():
            # Calculate node's x,y coordinates
            x = PAIRWISE_X_INTERVAL * self.max_level + PAIRWISE_X_INTERVAL

            if align == RIGHT:
                x = self.width - PAIRWISE_RIGHT_PADDING - x

            node.pos_in_pairwise = Point(x, self.y)  # Node pos_in_pairwiseition

            self.y += RT_Y_INTERVAL

            # Mouse click range of current leaf node
            if align == LEFT:
                range_topL_in_pairwise = Point(PAIRWISE_X_INTERVAL * level + PAIRWISE_X_INTERVAL, node.pos_in_pairwise.y - NODE_BLOCK_HEIGHT)
                range_botR_in_pairwise = Point(node.pos_in_pairwise.x + PAIRWISE_X_INTERVAL + len(str(node.taxon.label)) * 12,
                                   range_topL_in_pairwise.y + NODE_BLOCK_HEIGHT + 3)
                node.mouse_range_in_pairwise = Block(range_topL_in_pairwise, range_botR_in_pairwise)
                node.pairwise_block = Block(range_topL_in_pairwise, range_botR_in_pairwise)

                max_x = node.pos_in_pairwise.x + PAIRWISE_X_INTERVAL + len(str(node.taxon.label)) * 10
                if range_botR_in_pairwise.x > self.left_tree_width:
                    self.left_tree_width = max_x

            elif align == RIGHT:
                range_topL_in_pairwise = Point(node.pos_in_pairwise.x - PAIRWISE_X_INTERVAL - (len(str(node.taxon.label))) * 12, node.pos_in_pairwise.y -
                                                                       NODE_BLOCK_HEIGHT)
                x = self.width - PAIRWISE_RIGHT_PADDING - (PAIRWISE_X_INTERVAL * level + PAIRWISE_X_INTERVAL)
                range_botR_in_pairwise = Point(x,range_topL_in_pairwise.y + NODE_BLOCK_HEIGHT + 3)
                node.mouse_range_in_pairwise = Block(range_topL_in_pairwise, range_botR_in_pairwise)
                node.pairwise_block = Block(range_topL_in_pairwise, range_botR_in_pairwise)

                if range_topL_in_pairwise.x < self.right_tree_width:
                    self.right_tree_width = range_topL_in_pairwise.x

            self.draw_leaf_node(draw_canvas, node,align)

            # Testing
            # self.draw_dots(node.pos_in_pairwise.x + X_INTERVAL - 5 , node.pos_in_pairwise.y - NODE_BLOCK_HEIGH  T)
            # self.draw_dots(node.pos_in_pairwise.x + X_INTERVAL + len(str(node.taxon.label)) * 9, node.pos_in_pairwise.y + 2)
            # self.draw_rec(node.pos_in_pairwise.x + X_INTERVAL - 2,node.pos_in_pairwise.y - NODE_BLOCK_HEIGHT,X_INTERVAL + len(str(node.taxon.label)) * 9,NODE_BLOCK_HEIGHT+2)
            # self.draw_frame(node.pairwise_block.topL.x,node.pairwise_block.topL.y,node.pairwise_block.width,node.pairwise_block.height,RED)

        else:
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


        # Testing
        # self.draw_dots(node.pos_in_pairwise.x,node.pos_in_pairwise.y - X_INTERVAL)
        # self.draw_dots(node.pos_in_pairwise.x + X_INTERVAL, node.pos_in_pairwise.y)
        # self.draw_dots(node.pos_in_pairwise.x + X_INTER0L.y,8,node.range_botR_in_pairwise.y - node.range_topL_in_pairwise.y,BEIGE)
    def draw_leaf_node(self,canvas,node,align):
        canvas.text_align = align
        canvas.font = LEAF_FONT
        canvas.fill_style = BLACK
        if align == LEFT:
            canvas.fill_text(node.taxon.label, node.pos_in_pairwise.x + PAIRWISE_X_INTERVAL, node.pos_in_pairwise.y)
        elif align == RIGHT:
            canvas.fill_text(node.taxon.label, node.pos_in_pairwise.x - PAIRWISE_X_INTERVAL, node.pos_in_pairwise.y)
    def filter_node_selected(self,node, x , y):
        # if node.mouse_range_in_pairwise.check_inRange(Point(x,y)):
        if node.mouse_range_in_pairwise.topL.x <= x <= node.mouse_range_in_pairwise.botR.x and node.mouse_range_in_pairwise.topL.y <= y <= node.mouse_range_in_pairwise.botR.y:
            return True
        return False
    def mouse_clicked(self, x, y):
        if x < self.left_tree_width:
            node_selected = self.rt.find_node(lambda node: self.filter_node_selected(node, x,y))
        # elif x >= self.right_tree_width:
        #     node_selected = self.tc_tree.find_node(lambda node: self.filter_node_selected(node, x, y))

        if not node_selected:
            return

        # Testing
        # self.draw_frame(node_selected.pairwise_block.topL.x,node_selected.pairwise_block.topL.y,10,10,BLACK)

        self.adPy.select_subtree_from_tree(node_selected)
        # self.draw_subtree_block(node_selected,type=RT)
    def mouse_hover(self, x, y):
        type = None

        if x <= self.left_tree_width:
            node_selected = self.rt.find_node(lambda node: self.filter_node_selected(node, x, y))
            type = RT
        elif x >= self.right_tree_width:
            node_selected = self.tc_tree.find_node(lambda node: self.filter_node_selected(node, x, y))
            type = TC

        if not node_selected:
            self.node_hover = node_selected
            self[self.HOVER_NODE_LAYER].clear()
            self[self.HOVER_NODE_LAYER].flush()

        if node_selected != self.node_hover:
            self[self.HOVER_NODE_LAYER].clear()
            self[self.HOVER_NODE_LAYER].flush()

        if node_selected and node_selected != self.node_hover:
            self.node_hover = node_selected
            self.draw_hover_block(type,node=self.node_hover)
            self.draw_node_details(type,node_selected)
    def draw_node_details(self,type,node):
        # support , length , exact match
        if node.is_leaf():
            support = '-'
        else:
            support = node.label

        length = node.edge.length if node.edge.length != None else 0

        if type == TC and len(self.adPy.tc) > 0:
            exact_match = format(node.exact_match / len(self.ad_Py.tc) * 100, ".2f")
        else:
            exact_match = node.corr_similarity

        self[self.HOVER_NODE_LAYER].begin_path()
        self[self.HOVER_NODE_LAYER].fill_style = BLACK
        self[self.HOVER_NODE_LAYER].font = LEAF_FONT

        label_x = self.left_tree_width + ((self.right_tree_width - self.left_tree_width - 50) / 2)
        label_y = node.pos_in_pairwise.y
        if node.pos_in_pairwise.y + 40 > self.height:
            label_y = self.height - 40

        self[self.HOVER_NODE_LAYER].fill_text("Support : " + support, label_x  , label_y)
        self[self.HOVER_NODE_LAYER].fill_text("Length :  " + str(length), label_x , label_y+ 15)
        self[self.HOVER_NODE_LAYER].fill_text("Exact Match :  " + str(exact_match) + "%", label_x , label_y + 30)
    def draw_hover_block(self,type,node,corr=False):
        if type == RT:
            if not node.is_leaf():
                self.draw_rec(node.pos_in_pairwise.x + 2, node.pos_in_pairwise.y - RT_X_INTERVAL + 2 , RT_X_INTERVAL,
                              RT_X_INTERVAL - 2, BLACK, layer_index=self.HOVER_NODE_LAYER)
            else:
                width = node.pos_in_pairwise.x - node.mouse_range_in_pairwise.topL.x
                self.draw_rec(node.mouse_range_in_pairwise.topL.x, node.pos_in_pairwise.y - RT_X_INTERVAL,width,
                              RT_X_INTERVAL - 2, BLACK, layer_index = self.HOVER_NODE_LAYER)
            if corr:
                return
            else:
                self.draw_hover_block(type=TC,node=node.corr[self.tc_tree.id],corr=True)

        else:
            if not node.is_leaf():
                self.draw_rec(node.pos_in_pairwise.x - RT_X_INTERVAL - 2, node.pos_in_pairwise.y - RT_X_INTERVAL + 2 ,
                              RT_X_INTERVAL,
                              RT_X_INTERVAL - 2, BLACK, layer_index=self.HOVER_NODE_LAYER)
            else:
                width = node.mouse_range_in_pairwise.botR.x - node.pos_in_pairwise.x
                self.draw_rec(node.pos_in_pairwise.x, node.pos_in_pairwise.y - RT_X_INTERVAL, width,
                              RT_X_INTERVAL - 2, BLACK, layer_index = self.HOVER_NODE_LAYER)
            if corr:
                return
            else:
                self.draw_hover_block(type=RT,node=node.corr,corr=True)

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

    def generate_subtree_block(self,node_selected,select_tree=RT):
        if select_tree == RT:
            if node_selected.is_leaf():
                rec_x = node_selected.pairwise_block.topL.x + RT_X_INTERVAL - 5
            else:
                rec_x = node_selected.mouse_range_in_pairwise.topL.x

            rec_y = node_selected.pairwise_block.topL.y
            rec_width = self.left_tree_width + RT_X_INTERVAL - rec_x
            rec_height = node_selected.pairwise_block.botR.y - rec_y

        elif select_tree == TC:
            rec_x = self.left_tree_width + TREE_PADDING
            rec_y = node_selected.pairwise_block.topL.y
            if node_selected.is_leaf():
                rec_width = node_selected.mouse_range_in_pairwise.botR.x - rec_x - RT_X_INTERVAL
            else:
                rec_width = node_selected.mouse_range_in_pairwise.botR.x - rec_x

            rec_height = node_selected.pairwise_block.botR.y - rec_y

        new_block = Block(Point(rec_x, rec_y), Point(rec_x + rec_width, rec_y + rec_height))

        return new_block

    # def generate_mix_block(self, block,subtree_list,new_subtree):
    #
    #     if len(block.segment_list) == 0:  # Need to set current block as segment
    #         new_segment = Subtree_Block_Segment(block,subtree_list[0],topL=block.topL)
    #         block.segment_list.append(new_segment)
    #
    #     height = block.height / len(subtree_list)
    #     y = block.topL.y
    #
    #     for segment in block.segment_list:
    #         segment.height = height
    #         y += height
    #
    #     new_segment = Subtree_Block_Segment(block,new_subtree,topL=Point(block.topL.x,y),botR=block.botR)

    def remove_escape_taxa(self,subtree_label):
        for escape_taxa,subtree_list in self.escape_taxa_list.items():
            for subtree in subtree_list:
                if subtree_label == subtree.label:
                    subtree_list.remove(subtree)
                    if len(subtree_list) == 0:
                        del self.escape_taxa_list[escape_taxa]

    def draw_escape_taxa(self):
        self[self.ESCAPE_TAXA_LAYER].clear()

        for escape_taxa,subtree_list in self.escape_taxa_list.items():
            block = escape_taxa.escape_taxa_block
            if len(subtree_list) == 1:
                subtree = subtree_list[0]
                self.draw_rec(block.topL.x, block.topL.y, block.width, block.height,subtree.color,layer_index=self.ESCAPE_TAXA_LAYER)
            else:
                segment_width = block.width / len(subtree_list)
                for index,subtree in enumerate(subtree_list):
                    self.draw_rec(block.topL.x + segment_width * index, block.topL.y, segment_width, block.height,subtree.color,layer_index=self.ESCAPE_TAXA_LAYER)

    def record_escape_taxa(self,node_selected,subtree):
        node_selected.escape_taxa_block = self.generate_subtree_block(node_selected,select_tree=TC)
        self.escape_taxa_subtree.append(subtree.label)
        if node_selected not in self.escape_taxa_list:
            self.escape_taxa_list[node_selected] = []

        self.escape_taxa_list[node_selected].append(subtree)

    def draw_subtree_block(self,node_selected,select_tree=None,new_subtree=None,subtree_from_rt=None):
        # new_block = self.generate_subtree_block(node_selected,select_tree=select_tree)
        #
        # new_subtree.set_pairwise_block(new_block)
        #
        # # Get the available layer_index and draw block on the canvas
        # layer_index = self.get_layer_index(new_subtree,select_tree = select_tree)
        # self[layer_index].clear()
        #
        # self.draw_rec(new_block.topL.x, new_block.topL.y, new_block.width, new_block.height, new_subtree.color,
        #               layer_index=layer_index)
        #
        # if select_tree == RT:
        #     self.write_block_label(new_subtree, layer_index,align=LEFT)
        #     self.rt_layer_occupied[layer_index] = 1
        #
        # elif select_tree == TC:
        #     self.write_block_label(new_subtree, layer_index,align=RIGHT)
        #     self.tc_layer_occupied[layer_index-6] = 1
        #
        # self[layer_index].flush()
        # self[layer_index].label = new_subtree.label

        new_block = self.generate_subtree_block(node_selected, select_tree=select_tree)
        new_subtree.set_pairwise_block(new_block)

        if select_tree == RT:
            first_layer = 0
            last_layer = 5
        elif select_tree == TC:
            first_layer = 6
            last_layer = 11

        if select_tree == RT:
            self.draw_subtree_list = sorted(self.adPy.subtree_list, key=lambda x: x.block_size,reverse=True)
        else:
            self.draw_subtree_list = sorted(self.tc_subtree_list, key=lambda x: x.block_size,reverse=True)

        subtree_cnt = len(self.draw_subtree_list)
        tmp_sorted_layer_list = []
        for i in range(first_layer,first_layer + subtree_cnt):
            if select_tree == RT:
                subtree = self.draw_subtree_list[i]
            else:
                subtree = self.draw_subtree_list[i - 6]

            subtree_block = subtree.pairwise_block
            tmp_sorted_layer_list.append(subtree.label)

            self[i].clear()
            self.draw_rec(subtree_block.topL.x, subtree_block.topL.y, subtree_block.width, subtree_block.height,subtree.color,
                          layer_index=i)

            # if select_tree == TC and node_selected.corr_similarity < 1.0:
            #     self.adPy.output.append(f"similarity = {node_selected.corr_similarity}")
            #     for leaf_node in subtree_from_rt.root.leaf_nodes():
            #         if leaf_node.corr[self.tc_tree.id].selected is False:
            #             escape_taxa_block = self.generate_subtree_block(leaf_node.corr[self.tc_tree.id])
            #             escape_taxa_block.topL.x = subtree_block.topL.x - LABEL_MAX_WIDTH
            #             escape_taxa_block.width = subtree_block.width - LABEL_MAX_WIDTH
            #             self.draw_rec(escape_taxa_block.topL.x, escape_taxa_block.topL.y, escape_taxa_block.width,
            #                           escape_taxa_block.height,subtree.color,layer_index=i)

            if select_tree == RT:
                self.write_block_label(subtree, i,align=LEFT)

            elif select_tree == TC:
                self.write_block_label(subtree, i,align=RIGHT)

            self.adPy.output.append(f"after write label")

            self[i].label = subtree.label
            self[i].flush()

        if select_tree == RT:
            self.rt_sorted_layer_list = tmp_sorted_layer_list
        elif select_tree == TC:
            self.tc_sorted_layer_list = tmp_sorted_layer_list

    def remove_subtree_block(self,subtree):
        # Clear corresponding rt layer
        try:
            tc_clear_layer = None
            clear_layer = self.rt_sorted_layer_list.index(subtree.label)

            # Testing
            # tmp_canvas = Canvas(width=self.width,height=self.height)
            # tmp_canvas.draw_image(self[clear_layer],0,0)
            # tmp_canvas2 = Canvas(width=self.width, height=self.height)
            # tmp_canvas2.draw_image(self[1], 0, 0)
            # self.adPy.output.append(tmp_canvas)
            # self.adPy.output.append(tmp_canvas2)

            self[clear_layer].clear()

            if self.tc_tree:
                tc_clear_layer = self.tc_sorted_layer_list.index(subtree.label)
                self[tc_clear_layer + 6].clear()

            tc_subtree = subtree.corresponding_tc_subtree

            if tc_subtree.label in self.escape_taxa_subtree:
                self.remove_escape_taxa(subtree.label)
                self.escape_taxa_subtree.remove(tc_subtree.label)

            if hasattr(tc_subtree.root, 'duplicate_subtree') and tc_subtree.root.duplicate_subtree:
                self.adPy.output.append("Has duplicate subtree")
                self.adPy.output.append(tc_subtree.root.subtree)
                self.adPy.output.append(tc_subtree)
                tmp_subtree_list = tc_subtree.root.subtree
                tc_subtree.root.subtree.remove(tc_subtree)
                self.adPy.output.append("Before check redraw")
                for tree in tc_subtree.root.subtree:
                    if tree.pairwise_block.topL.x > tc_subtree.pairwise_block.topL.x:
                        tree.pairwise_block.topL.x -= LABEL_MAX_WIDTH
                        tree.pairwise_block.width += LABEL_MAX_WIDTH

                        # redraw_layer = self.tc_sorted_layer_list.index(tree.label)
                        # self.adPy.output.append(f"Redraw mix block {tree.label}in layer{redraw_layer+6}")
                        # self[redraw_layer+6].clear()
                        # self.draw_rec(tree.pairwise_block.topL.x, tree.pairwise_block.topL.y,
                        #               tree.pairwise_block.width, tree.pairwise_block.height,
                        #               tree.color,
                        #               layer_index=redraw_layer+6)

                self.adPy.output.append("After check redraw")
                if len(tmp_subtree_list) == 1:
                    tc_subtree.root.subtree = tmp_subtree_list[0]
                    tc_subtree.root.duplicate_subtree = False
                self.adPy.output.append("After change list to object")
            else:
                tc_subtree.root.subtree = None
                tc_subtree.root.selected = None

            self.tc_subtree_list.remove(tc_subtree)

            self.draw_escape_taxa()

        except Exception as err:
            self[-1].fill_text("Error in remove subtree", self.width - 200, 20)

    def adjust_block_size(self,new_subtree,select_tree):
        new_block = new_subtree.pairwise_block
        # 1. To check smallest nested block
        if select_tree == RT:
            tmp_subtree_list = sorted(self.adPy.subtree_list, key=lambda x: (x.block_size),
                                                        reverse=False)
        elif select_tree == TC:
            tmp_subtree_list = sorted(self.tc_subtree_list, key=lambda x: (x.block_size),
                                      reverse=False)

        for subtree in tmp_subtree_list:
            if subtree.label == new_subtree.label:
                continue

            if subtree.pairwise_block.check_nested_block(new_block):
                if select_tree == RT:
                    new_block.botR.x = subtree.pairwise_block.botR.x
                    new_block.calculate_width_height()
                else:
                    new_block.topL.x = subtree.pairwise_block.topL.x
                    new_block.calculate_width_height()

                if subtree.pairwise_block.topL.y == new_block.topL.y:
                    if select_tree == RT:
                        new_block.botR.x -= LABEL_MAX_WIDTH
                        new_block.width -= LABEL_MAX_WIDTH
                        new_block.calculate_width_height()
                    else:
                        new_block.topL.x += LABEL_MAX_WIDTH
                        new_block.width -= LABEL_MAX_WIDTH

                break

        # for index, tree in enumerate(tmp_subtree_list):
        #     # If current is nested beneath new block
        #     if new_block.check_nested_block(tree.pairwise_block, select_tree=select_tree):
        #         if tree.pairwise_block.topL.y == new_block.topL.y:
        #             if select_tree == RT:
        #                 if tree.pairwise_block.botR.x < new_block.botR.x:
        #                     checked = True
        #                 else:
        #                     break
        #             else:
        #                 if tree.pairwise_block.topL.x > new_block.botR.x:
        #                     checked = True
        #                 else:
        #                     break
        # if checked:
        #     return new_block

        # tmp_subtree_list = sorted(tmp_subtree_list, key=lambda x: x.block_size, reverse=False)
        # # tmp_subtree_list = sorted(tmp_subtree_list, key=lambda x: x.block_size, reverse=False)
        # self.adPy.output.append(tmp_subtree_list)
        # for index,tree in enumerate(tmp_subtree_list):
        #     # Check if new block nested beneath this subtree's block
        #     self.adPy.output.append(f"Check subtree: {tree.label}")
        #     if tree.pairwise_block.check_nested_block(new_block,select_tree = select_tree):
        #         nested_index = index
        #         self.adPy.output.append(f"new_block nested in current subtree: {tree.label}")
        #
        #         if select_tree == RT:
        #             new_block.botR.x = tree.pairwise_block.botR.x
        #             new_block.calculate_width_height()
        #         else:
        #             self.adPy.output.append("check new block x ")
        #             new_block.topL.x = tree.pairwise_block.topL.x
        #             new_block.calculate_width_height()
        #
        #         # Check if y_coor is same → block's label overlapping
        #         if tree.pairwise_block.topL.y == new_block.topL.y:
        #             self.adPy.output.append(f"new_block overlap with current subtree : {tree.label}")
        #             if select_tree == RT:
        #                 new_block.botR.x -= LABEL_MAX_WIDTH
        #                 new_block.width -= LABEL_MAX_WIDTH
        #             else:
        #                 new_block.topL.x += LABEL_MAX_WIDTH
        #                 new_block.width -= LABEL_MAX_WIDTH
        #
        #         break
        #
        # if nested_index is not None:
        #     return new_block
        #
        # self.adPy.output.append(f"Middle of adjust block size, x = {new_block.topL.x}")
        # tmp_subtree_list = sorted(tmp_subtree_list, key=lambda x: x.block_size, reverse=True)
        # interval = 0
        # nested_block = None
        # # self.adPy.output.append("loop2")
        # for index,subtree in enumerate(tmp_subtree_list):
        #     if new_block.check_nested_block(subtree.pairwise_block,select_tree = select_tree):
        #         self.adPy.output.append("current_tree nested in new_block")
        #         if subtree.pairwise_block.topL.y == new_block.topL.y:
        #             interval += 1
        #             if select_tree == RT:
        #                 subtree.pairwise_block.botR.x = new_block.botR.x - interval * LABEL_MAX_WIDTH
        #                 subtree.pairwise_block.calculate_width_height()
        #             else:
        #                 subtree.pairwise_block.topL.x = new_block.topL.x + interval * LABEL_MAX_WIDTH
        #                 subtree.pairwise_block.calculate_width_height()
        #
        #             nested_block = subtree.pairwise_block
        #             self.adPy.output.append(f"nested_block = {subtree.label}")
        #         else:
        #
        #             if nested_block and nested_block.check_nested_block(subtree.pairwise_block,select_tree = select_tree):
        #
        #                 if select_tree == RT:
        #                     subtree.pairwise_block.botR.x = nested_block.botR.x
        #                     subtree.pairwise_block.calculate_width_height()
        #                 else:
        #                     subtree.pairwise_block.topL.x = nested_block.topL.x
        #                     subtree.pairwise_block.calculate_width_height()
        #             else:
        #                 if select_tree == RT:
        #                     subtree.pairwise_block.botR.x = new_block.botR.x
        #                     subtree.pairwise_block.calculate_width_height()
        #                 else:
        #                     subtree.pairwise_block.topL.x = new_block.topL.x
        #                     subtree.pairwise_block.calculate_width_height()
        #
        #         redraw_layer = tmp_sorted_layer_list.index(subtree.label)
        #         if select_tree == TC:
        #             redraw_layer += 6
        #
        #         self[redraw_layer].clear()
        #         self.draw_rec(subtree.pairwise_block.topL.x, subtree.pairwise_block.topL.y, subtree.pairwise_block.width, subtree.pairwise_block.height,
        #                       subtree.color,
        #                       layer_index=redraw_layer)
        #
        #         if interval < 3:
        #             align = LEFT
        #             if select_tree == TC:
        #                 align = RIGHT
        #
        #             self.write_block_label(subtree, redraw_layer,align=align)
        #         # self[redraw_layer].flush()
        # self.adPy.output.append(f"End of adjust block size, x = {new_block.topL.x}")

        # return new_block


    # Sort multicanvas layer and return available layer's index
    # Layer of larger subtree should below the smaller subtree
    def get_layer_index(self,subtree,select_tree):
        try:
            if select_tree == RT:
                self.rt_layer_block_list[subtree.label] = subtree.block_size  # { 'label' : block-size }
                tmp_layer_block_list = self.rt_layer_block_list
                tmp_layer_occupied = self.rt_layer_occupied
                self.rt_sorted_layer_list = sorted(self.rt_layer_block_list, key=lambda x: self.rt_layer_block_list[x],
                                              reverse=True)
                tmp_sorted_layer_list = self.rt_sorted_layer_list
            else:
                self.tc_layer_block_list[subtree.label] = subtree.block_size  # { 'label' : block-size }
                tmp_layer_block_list = self.tc_layer_block_list
                tmp_layer_occupied = self.tc_layer_occupied
                self.tc_sorted_layer_list = sorted(self.tc_layer_block_list, key=lambda x: self.tc_layer_block_list[x],
                                              reverse=True)
                tmp_sorted_layer_list = self.tc_sorted_layer_list

                # If no subtree selected
            if tmp_layer_occupied.count(1) == 0:
                if select_tree == RT:
                    return 0
                else:
                    return 6

            if all(subtree.block_size <= value for value in tmp_layer_block_list.values()):
                index = tmp_layer_occupied.index(0)
                if select_tree == TC:
                    index += 6

                return index

            # If need to sort multicanvas layer - shift layers in list to the right
            canvas_tmp = Canvas(width = self.width,height = self.height)
            next_index = tmp_layer_occupied.index(0)

            if select_tree == TC:
                next_index += 6

            subtree_index = tmp_sorted_layer_list.index(subtree.label)
            if select_tree == TC:
                subtree_index += 6
            for i in range(next_index, subtree_index , -1):
                canvas_tmp.clear()

                canvas_tmp.draw_image(self[i-1],0,0)
                self[i].clear()
                self[i].draw_image(canvas_tmp,0,0)

            if select_tree == TC:
                tmp_layer_occupied[next_index-6] = 1

            index = tmp_sorted_layer_list.index(subtree.label)
            if select_tree == TC:
                index += 6

            return index

        except Exception as err:
            self[-1].fill_text("Error in get layer index", self.width - 200, 20)

    # One subtree's block was removed, shift layers in list to the left
    def rearrange_canvas_layer(self,clear_layer_index,select_tree):
        self.adPy.output.append("select_tree = " + select_tree)

        if select_tree == RT:
            first_layer = 0
            last_layer = 5
        elif select_tree == TC:
            first_layer = 6
            last_layer = 11

        self.adPy.output.append("clear_layer_index = " + str(clear_layer_index))
        self.adPy.output.append("last_layer= " + str(last_layer))

        canvas_tmp = Canvas(width=self.width, height=self.height)

        # Move list to left
        for i in range(clear_layer_index,last_layer):
            canvas_tmp.clear()
            if i < last_layer-1:
                canvas_tmp.draw_image(self[i+1], 0, 0)
            self[i].clear()
            self[i].draw_image(canvas_tmp, 0, 0)
    def setup_tc_subtree_list(self):
        for index,subtree in enumerate(self.adPy.subtree_list):
            self.draw_tc_subtree_block(subtree)

        self.draw_escape_taxa()

    def check_block_exact_match(self,tc_subtree,rt_subtree):
        if tc_subtree.leaf_set.issubset(rt_subtree.leaf_set):
            return True
        else:
            return False


    def draw_tc_subtree_block(self,subtree):
        '''
        如果符合条件的话主要subtree就和平时一样画，如果发现有escape taxa的话另外处理，把escape taxa记录起来
        remove的时候如果有escape taxa就去检查，把subtree颜色从escape taxa里面删除（如果没有了颜色就删除escape taxa）
        同时duplicate subtree的mix block也记录起来 （换成像是adjust_block_size，可是只是缩短完全一致的block）
        '''
        self.adPy.output.append("check in draw_tc_block")
        corresponding_subtree = subtree.root.corr[self.tc_tree.id]

        new_subtree = ADpy.Subtree(label=subtree.label, belong_tree=self.tc_tree, root=corresponding_subtree,
                                   color=subtree.color)
        subtree.corresponding_tc_subtree = new_subtree

        new_subtree.get_leaf_nodes()
        subtree.get_leaf_nodes()
        exact_match = False
        self.adPy.output.append(new_subtree.leaf_set)
        self.adPy.output.append(subtree.leaf_set)
        if new_subtree.leaf_set == subtree.leaf_set:
            exact_match = True

        if not exact_match:
            self.adPy.output.append(f"Not exact match, subtree {subtree.label}")
            tmp_result = self.check_block_exact_match(new_subtree,subtree)
            if not tmp_result:
                # All leaf nodes as escape taxa
                for leaf_node in subtree.root.leaf_nodes():
                    tc_corr = leaf_node.corr[self.tc_tree.id]
                    self.adPy.output.append(f"Leaf node: {tc_corr.taxon.label}")
                    self.record_escape_taxa(tc_corr,new_subtree)
                return
            else:
                # Draw leaf nodes which not in new_subtree
                escape_taxa = subtree.leaf_set.difference(new_subtree.leaf_set)
                for leaf_node in subtree.root.leaf_nodes():
                    if leaf_node.taxon.label not in escape_taxa:
                        continue

                    tc_corr = leaf_node.corr[self.tc_tree.id]
                    self.record_escape_taxa(tc_corr, new_subtree)

        self.check_duplicate_subtree(corresponding_subtree, subtree, new_subtree)

        for leaf_node in corresponding_subtree.leaf_nodes():
            leaf_node.selected = True

        self.tc_subtree_list.append(new_subtree)
        self.draw_subtree_block(corresponding_subtree, select_tree=TC, new_subtree=new_subtree,subtree_from_rt=subtree)

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

    def compare_tc_tree(self,compare_tree=None):
        self.tc_tree = compare_tree
        if self.tc_tree:
            # Draw tc tree (align right)
            self.default_value(align=RIGHT)
            self.draw_tree(self[self.TC_LAYER],align=RIGHT,node=self.tc_tree.seed_node, level=0)
            self.adPy.output.append("before call setup subtree list")
            if len(self.adPy.subtree_list) > 0:
                self.setup_tc_subtree_list()

    def reset_subtree_canvas(self):
        for i in range(0,11):
            if i == 5:
                continue
            self[i].clear()
