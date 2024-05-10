from ipycanvas import Canvas, hold_canvas
from .myCanvas import *
from .myUtils import *
import math
import time

# Canvas for Reference Tree
class rtCanvas(MyCanvas):
    # Initial x,y coordinates
    x = 30
    y = 30
    SUBTREE_COMPARE_LAYER = -3
    NODE_HOVER_LAYER = -2
    RT_LAYER = -1

    def __init__(self,adPy,width,height,view_support,default_rt=None):
        super().__init__( 8, width = width, height = height)
        # Layer 7 - reference tree
        # Layer 6 - hover block
        self.adPy = adPy
        self.layer = 8
        self.rt = adPy.rt
        self.view_support = view_support
        self.tree_width = 1
        self.dupletree = None


        self.layer_block_list = {}
        self.layer_occupied = [0, 0, 0, 0, 0]
        self.subtree_label_used = self.adPy.subtree_label_used
        self.subtree_color_used = self.adPy.subtree_color_used
        self.sorted_layer_list = []

        self.node_hover = None

        self.rt.level = 0

        self.section_list = []
        self.inner_section = True
        self.last_draw_time = time.time()
        self.last_mouse_position = None

        self.setup_section_list()

        if not default_rt:
            self.default_rt = Canvas(width=self.width,height=self.height)
            with hold_canvas(self):
                self.draw_rt(self.default_rt,node=self.rt.seed_node)

            self.adPy.default_rt = self.default_rt


        self[self.RT_LAYER].draw_image(self.default_rt, 0, 0)

        self.on_mouse_down(self.mouse_clicked)
        self.on_mouse_move(self.mouse_hover)

    def setup_section_list(self):
        section_cnt = math.ceil(self.height/ MIN_SECTION_HEIGHT)

        for i in range(section_cnt):
            self.section_list.append([])
            if self.inner_section:
                for j in range(3):
                    self.section_list[i].append([])
            # section_top_border = i * MIN_SECTION_HEIGHT
            # new_dupletree = DupleTree(self,Point(0,section_top_border),Point(self.width/2,section_top_border +
            #                                                                  MIN_SECTION_HEIGHT))
            # self.section_list.append(new_dupletree)

    # Draw reference tree on rtcanvas
    def draw_rt(self,draw_canvas,node=None,level=0):

        if level == 0:
            draw_canvas.fill_style = BLACK
            draw_canvas.font = LEAF_FONT
            draw_canvas.fill_text(self.rt.name, self.x, self.y)
            self.y += RT_Y_INTERVAL + 5

        # Default : Tree root as reference tree's root
        if node is None:
            node = self.rt.seed_node


        node.is_missing = False
        for child in node.child_node_iter():
            self.draw_rt(draw_canvas=draw_canvas,node=child, level=level + 1)

        # Set canvas style
        draw_canvas.fill_style = BLACK
        draw_canvas.font = LEAF_FONT

        if level > self.rt.level:
            self.rt.level = level

        if node.is_leaf():
            # Calculate node's x,y coordinates
            if self.view_support:
                x = self.x * level
            else:
                x = RT_X_INTERVAL * level + RT_X_INTERVAL

            node.pos = Point(x,self.y)  # Node position

            self.y += RT_Y_INTERVAL
            node.level = level

            # Mouse click range of current leaf node
            range_topL = Point(node.pos.x, node.pos.y - NODE_BLOCK_HEIGHT)
            range_botR = Point(range_topL.x + RT_X_INTERVAL + len(str(node.taxon.label)) * 12,
                               range_topL.y + NODE_BLOCK_HEIGHT + 3)
            node.mouse_range = Block(range_topL,range_botR)
            node.block = Block(range_topL,range_botR)

            if range_botR.x > self.tree_width:
                self.tree_width = range_botR.x

            self.draw_leaf_node(draw_canvas,node)

            # Testing
            # self.draw_dots(node.pos.x + X_INTERVAL - 5 , node.pos.y - NODE_BLOCK_HEIGH  T)
            # self.draw_dots(node.pos.x + X_INTERVAL + len(str(node.taxon.label)) * 9, node.pos.y + 2)
            # self.draw_rec(node.pos.x + X_INTERVAL - 2,node.pos.y - NODE_BLOCK_HEIGHT,X_INTERVAL + len(str(node.taxon.label)) * 9,NODE_BLOCK_HEIGHT+2)

        else:
            child_nodes = node.child_nodes()

            # To determine the midpoint of vertical line connecting multiple children
            first_child = child_nodes[0] if child_nodes else None
            last_child = child_nodes[self.RT_LAYER] if child_nodes else None

            if self.view_support:
                x = self.x * level
            else:
                x = RT_X_INTERVAL * level + RT_X_INTERVAL

            tmp = first_child.pos.y + last_child.pos.y
            node.pos = Point(x, tmp / 2)

            # Mouse click range of current internal node
            range_topL = Point(node.pos.x + 3, first_child.pos.y)
            if self.view_support:
                range_botR = Point(self.x * (level+1), last_child.pos.y - RT_X_INTERVAL)
            else:
                range_botR = Point(range_topL.x + RT_X_INTERVAL, last_child.pos.y - RT_X_INTERVAL)
            node.mouse_range = Block(range_topL, range_botR)

            # Draw vertical branch
            draw_canvas.begin_path()
            draw_canvas.move_to(first_child.pos.x, first_child.pos.y - 5)
            draw_canvas.line_to(last_child.pos.x, last_child.pos.y - 5)
            draw_canvas.stroke()

            # Drawing horizontal branch
            draw_canvas.begin_path()
            draw_canvas.move_to(node.pos.x , node.pos.y - 5)
            if self.view_support:
                draw_canvas.line_to(self.x * (level+1), node.pos.y - 5)
            else:
                draw_canvas.line_to(RT_X_INTERVAL * (level + 1) + RT_X_INTERVAL, node.pos.y - 5)
            draw_canvas.stroke()

            if self.view_support:
                draw_canvas.fill_style = 'blue'
                if node.label:
                    draw_canvas.fill_text(node.label, node.pos.x + 3, node.pos.y - 8)
                draw_canvas.fill_style = BLACK

            # Testing
            # self.draw_dots(node.pos.x + X_INTERVAL - 4, node.pos.y - X_INTERVAL)
            # self.draw_dots(node.pos.x + X_INTERVAL + 4, node.pos.y )
            # self.draw_rec(node.range_topL.x,node.range_topL.y,8,node.range_botR.y - node.range_topL.y,BEIGE)

        # To calculate and record subtree's block size
        if node.parent_node:
            if hasattr(node.parent_node, 'block'):
                if node.block.topL.y < node.parent_node.block.topL.y:
                    node.parent_node.block.topL = node.block.topL

                if node.block.botR.y > node.parent_node.block.botR.y:
                    node.parent_node.block.botR = node.block.botR
            else:
                node.parent_node.block = Block(node.block.topL,node.block.botR)


        node.selected = False
        self.insert_node_section_list(node)
        # self.dupletree.insert(node=node,point=node.mouse_range.topL)
        # self.dupletree.insert(node=node, point=node.mouse_range.botR)
        # self.dupletree.insert(node=node, point=Point(node.mouse_range.botR.x, node.mouse_range.topL.y))
        # self.dupletree.insert(node=node, point= Point(node.mouse_range.topL.x, node.mouse_range.botR.y))
        # self.quads.insert(quads.Point(block.topL.x, block.topL.y, data=node))
        # self.quads.insert(quads.Point(block.topL.x, block.botR.y, data=node))
        # self.quads.insert(quads.Point(block.botR.x, block.topL.y, data=node))
        # self.quads.insert(quads.Point(block.botR.x, block.botR.y, data=node))

        # Testing
        # self.draw_dots(node.pos.x,node.pos.y - X_INTERVAL)
        # self.draw_dots(node.pos.x + X_INTERVAL, node.pos.y)
        # self.draw_dots(node.pos.x + X_INTER0L.y,8,node.range_botR.y - node.range_topL.y,BEIGE)

    def insert_node_section_list(self,node):
        top_section_index = math.floor(node.mouse_range.topL.y / MIN_SECTION_HEIGHT)
        bottom_section_index = math.floor(node.mouse_range.botR.y / MIN_SECTION_HEIGHT)

        if not self.inner_section:
            if top_section_index == bottom_section_index:
                self.section_list[top_section_index].append(node)
            else:
                for i in range(top_section_index, bottom_section_index + 1):
                    self.section_list[i].append(node)

            return

        if top_section_index == bottom_section_index:
            section_midpoint = top_section_index * MIN_SECTION_HEIGHT + MIN_SECTION_HEIGHT / 2
            if node.mouse_range.botR.y <= section_midpoint:
                self.section_list[top_section_index][1].append(node)
            elif node.mouse_range.topL.y >= section_midpoint:
                self.section_list[top_section_index][2].append(node)
            else:
                self.section_list[top_section_index][0].append(node)
        else:
            for i in range(top_section_index, bottom_section_index + 1):
                section_midpoint = i * MIN_SECTION_HEIGHT + MIN_SECTION_HEIGHT / 2
                if i == 0:
                    if node.mouse_range.topL.y <= section_midpoint:
                        self.section_list[i][1].append(node)
                    else:
                        self.section_list[i][2].append(node)
                elif i == bottom_section_index:
                    if node.mouse_range.topL.y <= section_midpoint:
                        self.section_list[i][1].append(node)
                    else:
                        self.section_list[i][2].append(node)
                else:
                    self.section_list[i][0].append(node)

    # Drawing leaf node on canvas
    def draw_leaf_node(self,canvas,node):
        canvas.font = LEAF_FONT
        canvas.fill_style = BLACK
        canvas.begin_path()
        canvas.move_to(node.pos.x, node.pos.y - 5)
        canvas.line_to(node.pos.x + RT_X_INTERVAL, node.pos.y - 5)
        canvas.stroke()
        canvas.fill_text(node.taxon.label, node.pos.x + RT_X_INTERVAL * 2, node.pos.y)

    # Filter function
    def filter_node_selected(self,node, x , y):
        if node.mouse_range.topL.x <= x <= node.mouse_range.botR.x and node.mouse_range.topL.y <= y <= node.mouse_range.botR.y:
            return True
        return False
    def filter_node_from_section_list(self,x,y):
        node_selected = None
        if x > self.tree_width:
            return None

        inner_section_check = [1,1,1]

        section_index = math.floor(y / MIN_SECTION_HEIGHT)
        inner_section_check[0] = 0
        for node in self.section_list[section_index][0]:
            if self.filter_node_selected(node, x, y):
                return node

        section_midpoint = section_index * MIN_SECTION_HEIGHT + MIN_SECTION_HEIGHT / 2
        if y <= section_midpoint:
            subsection_list = self.section_list[section_index][1]
            inner_section_check[1] = 0
        else:
            subsection_list = self.section_list[section_index][2]
            inner_section_check[2] = 0

        for node in subsection_list:
            if self.filter_node_selected(node, x, y):
                return node

        subsection_index = inner_section_check.index(1)
        subsection_list = self.section_list[section_index][subsection_index]
        for node in subsection_list:
            if self.filter_node_selected(node, x, y):
                return node

        return node_selected


    # Handling mouse click events in rt canvas
    def mouse_clicked(self, x, y):
        self.clear_subtree_compare_canvas()
        if self.node_hover and self.filter_node_selected(self.node_hover,x,y):
            node_selected = self.node_hover
        else:
            node_selected = self.filter_node_from_section_list(x,y)
        # node_selected = self.rt.find_node(lambda node: self.filter_node_selected(node, x,y))

        if node_selected:
            self.adPy.select_subtree_from_tree(node_selected)

    def mouse_hover(self, x, y):
        current_time = time.time()
        if current_time - self.last_draw_time > 0.1:
            self.last_draw_time = current_time
            node_selected = self.filter_node_from_section_list(x, y)

            # if (x, y) != self.last_mouse_position:
            #     self.last_mouse_position = (x, y)

        else:
            return

        if not node_selected:
            self.node_hover = node_selected
            self[self.NODE_HOVER_LAYER].clear()
            self[self.NODE_HOVER_LAYER].flush()

        if node_selected != self.node_hover:
            self[self.NODE_HOVER_LAYER].clear()
            self[self.NODE_HOVER_LAYER].flush()

        if node_selected and node_selected != self.node_hover:
            self.node_hover = node_selected
            self.draw_hover_block()
            self.draw_node_details(node_selected)

    def write_block_label(self,subtree,layer_index):
        node_amt = len(subtree.root.leaf_nodes())
        label_str = f"{subtree.label} : {node_amt}"
        subtree.label_width = (len(label_str) * 6)
        label_x = subtree.block.topL.x + subtree.block.width - subtree.label_width
        label_y = subtree.block.topL.y + 12

        self[layer_index].fill_style = BLACK
        self[layer_index].fill_text(label_str, label_x, label_y)

        # Testing
        # self[layer_index].fill_text("layer: " + str(layer_index), subtree.block.botR.x + 5,
        #                             subtree.block.topL.y + 5)

    def draw_node_details(self,node):
        # support , length , exact match
        if node.is_leaf():
            support = '-'
        else:
            support = node.label

        length = node.edge.length if node.edge.length != None else 0

        if len(self.adPy.tc) > 0:
            exact_match = format(node.exact_match / len(self.adPy.tc) * 100, ".2f")
        else:
            exact_match = '- '

        self[self.NODE_HOVER_LAYER].begin_path()
        self[self.NODE_HOVER_LAYER].fill_style = BLACK
        self[self.NODE_HOVER_LAYER].font = LEAF_FONT

        if self.view_support:
            label_pos = self.tree_width + 50
        else:
            label_pos = self.tree_width + 200

        label_y = node.pos.y
        if node.pos.y + 40 > self.height:
            label_y = self.height - 40

        self[self.NODE_HOVER_LAYER].fill_text("Support : " + support, label_pos , label_y)
        self[self.NODE_HOVER_LAYER].fill_text("Length :  " + str(length), label_pos, label_y + 15)
        self[self.NODE_HOVER_LAYER].fill_text("Exact Match :  " + str(exact_match) + "%", label_pos , label_y + 30)

    def draw_hover_block(self):
        if not self.node_hover.is_leaf() and self.view_support:
            self.draw_rec(self.node_hover.pos.x + RT_X_INTERVAL , self.node_hover.pos.y - RT_X_INTERVAL + 2 , RT_X_INTERVAL,
                          RT_X_INTERVAL - 2, BLACK, layer_index=-2)
        else:
            self.draw_rec(self.node_hover.pos.x, self.node_hover.pos.y - RT_X_INTERVAL, RT_X_INTERVAL - 2, RT_X_INTERVAL - 2, BLACK, layer_index = -2)

        self[self.NODE_HOVER_LAYER].flush()

    def generate_subtree_block(self,node_selected):
        if node_selected.is_leaf():
            rec_x = node_selected.block.topL.x + RT_X_INTERVAL - 5
        else:
            rec_x = node_selected.mouse_range.topL.x + 2

        rec_y = node_selected.block.topL.y
        rec_width = self.tree_width - rec_x
        rec_height = node_selected.block.botR.y - rec_y

        new_block = Block(Point(rec_x, rec_y), Point(rec_x + rec_width, rec_y + rec_height))
        return new_block

    def draw_subtree_block(self,node_selected,new_subtree=None):

        new_block = self.generate_subtree_block(node_selected)
        new_block = self.adjust_block_size(new_block)

        new_subtree.set_block(new_block)

        # Get the available layer_index and draw block on the canvas
        layer_index = self.get_layer_index(new_subtree)
        self[layer_index].clear()

        self.draw_rec(new_block.topL.x,new_block.topL.y,new_block.width,new_block.height,new_subtree.color,
                      layer_index =layer_index)

        self.write_block_label(new_subtree,layer_index)

        # Mark the layer as occupied and record corresponding subtree's label
        self.layer_occupied[layer_index] = 1
        self[layer_index].label = new_subtree.label

        self[layer_index].flush()

        # Testing
        # self[layer_index].fill_text("layer: " + str(layer_index), new_block.botR.x + 5, new_block.topL.y + 5)

    def remove_subtree_block(self,subtree):
        # Clear corresponding layer
        clear_layer = self.sorted_layer_list.index(subtree.label)
        self[clear_layer].clear()

        # Manage subtree list and multicanvas layer
        del self.layer_block_list[subtree.label]
        self.sorted_layer_list.remove(subtree.label)

        self.rearrange_canvas_layer(clear_layer_index=clear_layer)

        self.layer_occupied[len(self.sorted_layer_list)] = 0


    def adjust_block_size(self,new_block):
        nested_index = None  # smallest block that new_block nested
        nesting_block_index = None  # largest block that new block nesting
        size_lock = False

        # 1. To adjust new block's size regardless of whether the tag will be overlapped
        if len(self.adPy.subtree_list) <= 1:
            return new_block

        self.adPy.subtree_list = sorted(self.adPy.subtree_list, key=lambda x: x.block_size,reverse=False)
        for index,tree in enumerate(self.adPy.subtree_list):

            # Check if new block nesting other existing subtree
            if new_block.check_nested_block(tree.block):
                nesting_block_index = index

                if not size_lock:
                    new_block.botR.x = tree.block.botR.x
                    new_block.calculate_width_height()

                if tree.block.topL.y == new_block.topL.y:
                    size_lock = True
                    new_block.botR.x += LABEL_MAX_WIDTH
                    new_block.width += LABEL_MAX_WIDTH

                continue

            # Check if new block nested beneath this subtree's block
            if tree.block.check_nested_block(new_block):
                nested_index = index

                new_block.botR.x = tree.block.botR.x
                new_block.calculate_width_height()

                # Check if y_coor is same → block's label overlapping
                if tree.block.topL.y == new_block.topL.y:
                    new_block.botR.x -= LABEL_MAX_WIDTH
                    new_block.width -= LABEL_MAX_WIDTH

                break

        # If new block is completely independent
        if nested_index == None:
            return new_block

        # 2. Check whether new block's label will overlap with existing block's → adjust existing block size
        if nested_index != None:
            # Shrink blocks within the new block
            interval = 0
            for i in range(nested_index, 0 ,-1):
                subtree = self.adPy.subtree_list[i-1]
                if subtree.block.topL.y == new_block.topL.y:
                    interval += 1
                    subtree.block.botR.x = new_block.botR.x - interval * LABEL_MAX_WIDTH
                    subtree.block.calculate_width_height()

                    redraw_layer = self.sorted_layer_list.index(subtree.label)
                    self[redraw_layer].clear()
                    self.draw_rec(subtree.block.topL.x, subtree.block.topL.y, subtree.block.width, subtree.block.height,
                                  subtree.color,
                                  layer_index=redraw_layer)

                    self.write_block_label(subtree, redraw_layer)
                    # self[redraw_layer].flush()

            return new_block
    # Sort multicanvas layer and return available layer's index
    # Layer of larger subtree should below the smaller subtree
    def get_layer_index(self,subtree):
        self.layer_block_list[subtree.label] = subtree.block_size  # { 'label' : block-size }
        self.sorted_layer_list = sorted(self.layer_block_list, key=lambda x: self.layer_block_list[x], reverse=True)

        # If no subtree selected
        if self.layer_occupied.count(1) == 0:
            return 0

        # If new_subtree is smaller than all existing subtree
        if all(subtree.block_size <= value for value in self.layer_block_list.values()):
            return self.layer_occupied.index(0)

        # If need to sort multicanvas layer - shift layers in list to the right
        canvas_tmp = Canvas(width = self.width,height = self.height)
        next_index = self.layer_occupied.index(0)
        for i in range(next_index, self.sorted_layer_list.index(subtree.label) , -1):
            canvas_tmp.clear()
            canvas_tmp.draw_image(self[i-1],0,0)
            self[i].clear()
            self[i].draw_image(canvas_tmp,0,0)

        self.layer_occupied[next_index] = 1
        return self.sorted_layer_list.index(subtree.label)

    # One subtree's block was removed, shift layers in list to the left
    def rearrange_canvas_layer(self,clear_layer_index):
        last_layer = 5

        canvas_tmp = Canvas(width=self.width, height=self.height)
        for i in range(clear_layer_index,last_layer):
            canvas_tmp.clear()
            if i < 4:
                canvas_tmp.draw_image(self[i+1], 0, 0)
            self[i].clear()
            self[i].draw_image(canvas_tmp, 0, 0)

    def reset_subtree_canvas(self):
        for i in range(0,5):
            self[i].clear()


    def clear_subtree_compare_canvas(self):
        self[self.SUBTREE_COMPARE_LAYER].clear()
        self[self.SUBTREE_COMPARE_LAYER].flush()

    def draw_subtree_compare_nodes(self,nodes_list):
        self.clear_subtree_compare_canvas()
        for node in nodes_list:
            rt_node = node.corr
            self.draw_nodes_dots(rt_node)

    def draw_nodes_dots(self,node):
        self[self.SUBTREE_COMPARE_LAYER].fill_style = BLACK
        self[self.SUBTREE_COMPARE_LAYER].stroke_style = BLACK
        radius = 4
        self[self.SUBTREE_COMPARE_LAYER].fill_arc(node.pos.x + RT_X_INTERVAL + 5 , node.pos.y - 5,
                                             radius, 0, 360)


        self[self.SUBTREE_COMPARE_LAYER].flush()





