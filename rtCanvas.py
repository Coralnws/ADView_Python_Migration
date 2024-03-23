from ipycanvas import Canvas,MultiCanvas
from AD_Py import *
from myCanvas import *
from Utils import *
import math

# Canvas for Reference Tree
class rtCanvas(MyCanvas):
    # Initial x,y coordinates
    x = 30
    y = 30

    # Subtree Block color
    subtree_color_list = [BLUE,ORANGE,GREEN,RED,PURPLE]

    subtree_label = ['A','B','C','D','E']

    def __init__(self,ad_Py,width,height,view_support,tc):
        super().__init__( 7, width = width, height = height)
        # Layer 7 - reference tree
        # Layer 6 - hover block
        self.ad_Py = ad_Py
        self.rt = ad_Py.rt
        self.view_support = view_support
        self.tree_width = 1

        # Manage Multicanvas
        self.layer_occupied = [0, 0, 0, 0, 0]
        self.layer_block_list = {}

        # Manage Subtree
        self.subtree_list = []
        self.subtree_label_used = [1, 1, 1, 1, 1]
        self.subtree_color_used = [1, 1, 1, 1, 1]

        self.node_hover = None

        self.draw_rt(node=self.rt.seed_node)
        self.on_mouse_down(self.mouse_clicked)
        self.on_mouse_move(self.mouse_hover)


    # Draw reference tree on rtcanvas
    def draw_rt(self,node=None,level=0):
        # Default : Tree root as reference tree's root
        if node is None:
            node = self[-1].tree.seed_node

        for child in node.child_node_iter():
            self.draw_rt(child, level + 1)

        # Set canvas style
        self[-1].fill_style = BLACK
        self[-1].font = LEAF_FONT

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

            self.draw_leaf_node(node)

            # Testing
            # self.draw_dots(node.pos.x + X_INTERVAL - 5 , node.pos.y - NODE_BLOCK_HEIGH  T)
            # self.draw_dots(node.pos.x + X_INTERVAL + len(str(node.taxon.label)) * 9, node.pos.y + 2)
            # self.draw_rec(node.pos.x + X_INTERVAL - 2,node.pos.y - NODE_BLOCK_HEIGHT,X_INTERVAL + len(str(node.taxon.label)) * 9,NODE_BLOCK_HEIGHT+2)

        else:
            child_nodes = node.child_nodes()

            # To determine the midpoint of vertical line connecting multiple children
            first_child = child_nodes[0] if child_nodes else None
            last_child = child_nodes[-1] if child_nodes else None

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
            self[-1].begin_path()
            self[-1].move_to(first_child.pos.x, first_child.pos.y - 5)
            self[-1].line_to(last_child.pos.x, last_child.pos.y - 5)
            self[-1].stroke()

            # Drawing horizontal branch
            self[-1].begin_path()
            self[-1].move_to(node.pos.x , node.pos.y - 5)
            if self.view_support:
                self[-1].line_to(self.x * (level+1), node.pos.y - 5)
            else:
                self[-1].line_to(RT_X_INTERVAL * (level + 1) + RT_X_INTERVAL, node.pos.y - 5)
            self[-1].stroke()

            if self.view_support:
                self[-1].fill_style = 'blue'
                if node.label:
                    self[-1].fill_text(node.label, node.pos.x + 3, node.pos.y - 8)
                self[-1].fill_style = BLACK

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

        # Testing
        # self.draw_dots(node.pos.x,node.pos.y - X_INTERVAL)
        # self.draw_dots(node.pos.x + X_INTERVAL, node.pos.y)
        # self.draw_dots(node.pos.x + X_INTER0L.y,8,node.range_botR.y - node.range_topL.y,BEIGE)

    # Drawing leaf node on canvas
    def draw_leaf_node(self,node):
        self[-1].fill_style = BLACK
        self[-1].begin_path()
        self[-1].move_to(node.pos.x, node.pos.y - 5)
        self[-1].line_to(node.pos.x + RT_X_INTERVAL, node.pos.y - 5)
        self[-1].stroke()
        self[-1].fill_text(node.taxon.label, node.pos.x + RT_X_INTERVAL, node.pos.y)

    # Filter function
    def filter_node_selected(self,node, x , y):
        # if node.mouse_range.check_inRange(Point(x,y)):
        if node.mouse_range.topL.x <= x <= node.mouse_range.botR.x and node.mouse_range.topL.y <= y <= node.mouse_range.botR.y:
            return True
        return False

    # Handling mouse click events in rt canvas
    def mouse_clicked(self, x, y):
        node_selected = self.rt.find_node(lambda node: self.filter_node_selected(node, x,y))
        if not node_selected:
            return

        # Testing
        # self.draw_rec(node_selected.block.topL.x,node_selected.block.topL.y,10,10,BLACK)

        self.draw_subtree_block(node_selected)

    def mouse_hover(self, x, y):
        node_selected = self.rt.find_node(lambda node: self.filter_node_selected(node, x,y))

        if not node_selected:
            self.node_hover = node_selected
            self[-2].clear()
            self[-2].flush()

        if node_selected != self.node_hover:
            self[-2].clear()
            self[-2].flush()

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

    def draw_node_details(self,node):
        # support , length , exact match
        if node.is_leaf():
            support = '-'
        else:
            support = node.label

        length = node.edge.length if node.edge.length != None else 0

        if len(self.ad_Py.tc) > 0:
            exact_match = format(node.exact_match / len(self.ad_Py.tc) * 100, ".2f")
        else:
            exact_match = '- '

        self[-2].begin_path()
        self[-2].fill_style = BLACK
        self[-2].font = LEAF_FONT

        if self.view_support:
            label_pos = self.tree_width + 50
        else:
            label_pos = self.tree_width + 200

        self[-2].fill_text("Support : " + support, label_pos , node.pos.y)
        self[-2].fill_text("Length :  " + str(length), label_pos, node.pos.y + 15)
        self[-2].fill_text("Exact Match :  " + str(exact_match) + "%", label_pos , node.pos.y + 30)

    def draw_hover_block(self):
        if not self.node_hover.is_leaf() and self.view_support:
            self.draw_rec(self.node_hover.pos.x + RT_X_INTERVAL , self.node_hover.pos.y - RT_X_INTERVAL + 2 , RT_X_INTERVAL,
                          RT_X_INTERVAL - 2, BLACK, layer_index=-2)
        else:
            self.draw_rec(self.node_hover.pos.x, self.node_hover.pos.y - RT_X_INTERVAL, RT_X_INTERVAL - 2, RT_X_INTERVAL - 2, BLACK, layer_index = -2)

        self[-2].flush()

    def draw_subtree_block(self,node_selected):
        # If the node was not selected → draw block and set node as selected_subtree
        if not hasattr(node_selected, 'selected') or node_selected.selected == False:
            # Ignore if 5 subtree had selected
            if len(self.subtree_list) >= 5:
                return

            # Calculate block position and size
            if node_selected.is_leaf():
                rec_x = node_selected.block.topL.x + RT_X_INTERVAL - 5
            else:
                rec_x = node_selected.mouse_range.topL.x + 2

            rec_y = node_selected.block.topL.y
            rec_width = self.tree_width - rec_x
            rec_height = node_selected.block.botR.y - rec_y

            new_block = self.adjust_block_size(Block(Point(rec_x, rec_y), Point(rec_x + rec_width, rec_y + rec_height)))

            # Choose block color
            color_index = self.subtree_color_used.index(1)
            self.subtree_color_used[color_index] = 0
            color = self.subtree_color_list[color_index]

            # Record new subtree - choose subtree label
            subtree_label_index = self.subtree_label_used.index(1)
            self.subtree_label_used[subtree_label_index] = 0
            label = self.subtree_label[subtree_label_index]

            # Create new subtree (class from myTree.py)
            new_subtree = myTree.Subtree(label = label, rt = self, root = node_selected, color = color, block = new_block)
            self.subtree_list.append(new_subtree)

            # Get the available layer_index and draw block on the canvas
            layer_index = self.get_layer_index(new_subtree)
            self[layer_index].clear()

            self.draw_rec(new_block.topL.x,new_block.topL.y,new_block.width,new_block.height,color,layer_index = layer_index)
            # self.tree_width = ori_tree_width

            self.write_block_label(new_subtree,layer_index)

            # Mark the layer as occupied and record corresponding subtree's label
            self.layer_occupied[layer_index] = 1
            self[layer_index].label = new_subtree.label

            # Link node with subtree and mark node as selected
            node_selected.subtree = new_subtree
            node_selected.selected = True

            self[layer_index].flush()

        else:
            # If the node was selected → erase block and remove node from selected_subtree

            # Clear corresponding layer
            clear_layer = self.sorted_layer_list.index(node_selected.subtree.label)
            self[clear_layer].clear()

            # Manage subtree list and multicanvas layer
            del self.layer_block_list[node_selected.subtree.label]
            self.sorted_layer_list.remove(node_selected.subtree.label)

            self.rearrange_canvas_layer(clear_layer_index = clear_layer)

            self.layer_occupied[len(self.sorted_layer_list)] = 0
            self.subtree_list.remove(node_selected.subtree)

            # Release Color
            color_index = self.subtree_color_list.index(node_selected.subtree.color)
            self.subtree_color_used[color_index] = 1

            # Release Subtree Label
            label_index = self.subtree_label.index(node_selected.subtree.label)
            self.subtree_label_used[label_index] = 1

            # Remove subtree from node and mark node as not-selected
            node_selected.subtree = None
            node_selected.selected = False

    def adjust_block_size(self,new_block):
        nested_index = None  # smallest block that new_block nested
        nesting_block_index = None  # largest block that new block nesting
        size_lock = False

        # 1. To adjust new block's size regardless of whether the tag will be overlapped
        self.subtree_list = sorted(self.subtree_list, key=lambda x: x.block.height,reverse=False)


        for index,tree in enumerate(self.subtree_list):

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
                subtree = self.subtree_list[i-1]
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
