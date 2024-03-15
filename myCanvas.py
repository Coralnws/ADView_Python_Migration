import ipycanvas as ipc
from ipycanvas import Canvas,MultiCanvas
import myTree as myTree
from Utils import *

# Use Multicanvas which has multiple canvas layer
class MyCanvas(MultiCanvas):
    tree = None
    line_dashes = [[2, 2], [5, 5], [8, 8], [10, 10]]

    def __init__(self,layer,width, height):
        super().__init__(layer,width = width, height = height)

    def draw_dots(self,x,y,layer_index = -1):
        self[layer_index].fill_style = RED
        self[layer_index].fill_circle(x, y, 1)
        self[layer_index].fill_style = BLACK

    def draw_solid_line(self,head,tail,color,layer_index = -1):
        self[layer_index].stroke_style = color
        self[layer_index].begin_path()
        self[layer_index].move_to(head.x, head.y)
        self[layer_index].line_to(tail.x,tail.y)
        self[layer_index].stroke()

    def draw_dotted_line(self,head,tail,color,line_dash_index=0,layer_index = -1):

        self[layer_index].stroke_style = color
        self[layer_index].set_line_dash(self.line_dashes[line_dash_index])
        self[layer_index].begin_path()
        self[layer_index].move_to(head.x, head.y)
        self[layer_index].line_to(tail.x,tail.y)
        self[layer_index].stroke()

        self[layer_index].set_line_dash([])

    # Draw rectangle on specific layer
    def draw_rec(self,x,y,width,height,color,layer_index = -1):
        self[layer_index].fill_style = color
        self[layer_index].fill_rect(x, y, width,height)

    def draw_frame(self,x,y,width,height,color,layer_index = -1):
        self[layer_index].stroke_style = color
        self[layer_index].stroke_rect(x, y, width, height)

# Canvas for Reference Tree
class rtCanvas(MyCanvas):
    # Initial x,y coordinates
    x = 30
    y = 30

    # Subtree Block color
    subtree_color_list = [BLUE,ORANGE,GREEN,RED,PURPLE]

    subtree_label = ['A','B','C','D','E']

    def __init__(self,tree,width,height,view_support,tc):
        super().__init__( 7, width = width, height = height)
        # Layer 7 - reference tree
        # Layer 6 - hover block
        self.tree = tree
        self.rt = tree.rt
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

        if len(self.tree.tc) > 0:
            exact_match = format(node.exact_match / len(self.tree.tc) * 100, ".2f")
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


# Canvas for Tree Collection
class tcCanvas(MyCanvas):
    id = None
    ad_list = []
    AD_LAYER = -1

    def __init__(self,tree,width=1100, height=1100,scale=1.0,max_ad=None):
        super().__init__(5, width = width, height = height)
        self.ad_list = []
        self.scale = scale
        self.max_ad = max_ad
        self.default_value()

        # Initialization
        self.tree = tree
        self.rt = tree.rt
        self.tc = tree.tc
        self.subtree_list = tree.subtree_list
        self.block_minimum_height = BLOCK_MINIMUM_HEIGHT[len(self.subtree_list) % MAX_SUBTREE]

        for subtree in self.subtree_list:
            subtree.get_leaf_nodes()

        for index,tc_tree in enumerate(self.tc):
            # print("tc_tree " + str(index))
            self.construct_ad(tc_tree=tc_tree, level=1)
            # self.ad_list[index].plot_tree()

            # print("=============")

        self.preprocess_ad_tree()
        self.draw_ad_tree_rec()

        # print("after construct ad")

        for index,ad_tree in enumerate(self.ad_list):
            if self.max_ad:
                if index >= max_ad:
                    break

            # print("ad_tree " + str(index))
            # print("topL.x:" + str(ad_tree.topL.x))
            # print("topL.y:" + str(ad_tree.topL.y))
            # print("botR.x: " + str(ad_tree.botR.x))
            # print("botR.y: " + str(ad_tree.botR.y))

            self.draw_ad_tree(ad_tree)
            self[self.AD_LAYER].flush()
            # print("draw_ad")

    def default_value(self):
        self.padding = DEFAULT_AD_PADDING  # Padding between ADs
        self.ad_per_row = DEFAULT_AD_PER_ROW
        self.ad_width = DEFAULT_AD_WIDTH
        self.ad_height = DEFAULT_AD_HEIGHT

    # y-padding = 5
    # y-gap between block = 3
    # block space = ad_height - 10 - (block_amount-1) * 3
    # block_height = block_node_amount/node_amount * block space
    # independent leaf node(from subtree) has same size as blank block

    # Traverse from root, and see whether descendant is target subtree/leaf node(need to compare all target leaf node, in case there have
    # any independent leaf node), if not descendant, append a blank block
    # if has descendant, continue traverse along the specific child
    # if current node is root of target subtree → append a subtree block
    # if found a target leaf node (means independent leaf node), append a leaf node block which same size as blank block

    # Function to construct AD_Tree for trees from tree collection
    # Traverse from root, check all child nodes relations, if not BLANK block, traverse along the child
    def construct_ad(self,tc_tree,tc_node=None,current_ad_node=None,level=0):

        if tc_node is None and current_ad_node is None: # When start from root
            tc_node = tc_tree.seed_node

            ad = myTree.AD_Tree(id=tc_tree.id,tc_tree=tc_tree,tc_canvas=self)  # Create a new AD_Tree which belong to current tc_tree
            ad.root = myTree.AD_Node(node_or_block=tc_node,x_level=level,type=ROOT)  # Set ad_tree's root
            current_ad_node = ad.root
            self.ad_constructing = ad
            self.ad_list.append(ad)

        for child in tc_node.child_node_iter():
            # print("Child:")
            # print(child)

            result = NON_DESCENDANT
            for subtree in self.subtree_list:
                result = self.check_relation(subtree=subtree,node=child,tree_index = tc_tree.id)

                if result == INDEPENDENT_LEAF:
                    # print("INDEPENDENT_LEAF")
                    # width,height,subtree,ad_tree,type,topL=None, botR=None
                    independent_leaf_block = self.individual_block(subtree=subtree,type=INDIVIDUAL_LEAF_BLOCK)
                    independent_leaf_block.color = subtree.color
                    independent_leaf_block.root = child
                    independent_leaf_block.taxa_list.add(child.taxon.label)
                    independent_leaf_block.taxa_count = 1
                    new_ad_node = myTree.AD_Node(node_or_block=independent_leaf_block, x_level=level, type=LEAF)
                    self.ad_constructing.insert_node(current_ad_node,new_ad_node)
                    # current_ad_node.children.append(new_ad_node)

                    break

                if result == SUBTREE_ROOT:
                    # print("SUBTREE")

                    subtree_block = self.ad_subtree_block(subtree=subtree)
                    subtree_block.color = subtree.color
                    subtree_block.root = child

                    subtree_block.taxa_list = {leaf.taxon.label for leaf in child.leaf_nodes()}
                    subtree_block.taxa_count = len(subtree_block.taxa_list.intersection(subtree.leaf_set))

                    if subtree_block.taxa_count == len(subtree.leaf_set):
                        subtree_block.exact_match = True
                    else:
                        subtree_block.exact_match = False

                    new_ad_node = myTree.AD_Node(node_or_block=subtree_block, x_level=level, type=LEAF)

                    self.ad_constructing.insert_node(current_ad_node, new_ad_node)
                    # current_ad_node.children.append(new_ad_node)

                    break

                if result == IS_DESCENDANT:
                    # print("IS_DESCENDANT")
                    new_ad_node = myTree.AD_Node(node_or_block=child, x_level=level, type=INTERNAL)
                    self.ad_constructing.insert_node(current_ad_node, new_ad_node)
                    # current_ad_node.children.append(new_ad_node)
                    self.construct_ad(tc_tree=tc_tree,tc_node=child,current_ad_node=new_ad_node,level=level+1)

                    break

            if result == NON_DESCENDANT:
                # print("NON_DESCENDANT")
                non_descendant_block = self.individual_block(subtree=None,type=INDIVIDUAL_BLOCK)
                non_descendant_block.color = BLANK
                non_descendant_block.root = child
                new_ad_node = myTree.AD_Node(node_or_block=non_descendant_block, x_level=level, type=LEAF)
                self.ad_constructing.insert_node(current_ad_node, new_ad_node)
                # current_ad_node.children.append(new_ad_node)


        # Subtree : label,rt,root,color,block

    # Function to check whether current branch/node related with target subtree/taxa
    def check_relation(self,subtree,node,tree_index):
        # subtree = subtree from subtree_list
        # node = node from tc_tree
        # tree_index = tc_tree's id

        # If node from tc is target leaf node in this subtree
        if node.is_leaf():
            if node.taxon.label in subtree.leaf_set:
                return INDEPENDENT_LEAF

        # If cunode is subtree's root
        if node == subtree.root.corr[tree_index]:
            return SUBTREE_ROOT

        leaf_nodes = {leaf.taxon.label for leaf in node.leaf_nodes()}

        # If target taxa is descendant of node
        if leaf_nodes.intersection(subtree.leaf_set):
            return IS_DESCENDANT

        # If current branch has no target taxa
        return NON_DESCENDANT

    # Function for ad_drawing preprocess
    # Draw AD in tc_canvas
    def preprocess_ad_tree(self):
        # Comfirm canvas scale
        # Preprocess ad_tree
        # 1. set position and size for all ad_tree
        # 2. set width and height for individual_block in every ad_tree
        # 3.
        self.default_value()

        # Calculate AD size, AD per row and padding
        self.ad_width = DEFAULT_AD_WIDTH * self.scale
        self.ad_height = DEFAULT_AD_HEIGHT * self.scale
        self.ad_per_row = ((self.width - (3 * self.padding)) // self.ad_width)

        while (self.width - ( self.ad_width * self.ad_per_row + 2 * self.padding )) < (self.ad_per_row-2) * self.padding :
            self.ad_per_row -= 1

        # padding_space = self.width - 4 * self.padding - (self.ad_per_row * self.ad_width)
        # self.padding = padding_space // (self.ad_per_row-1)

        # print("ad_width:" + str(self.ad_width))
        # print("ad_height:" + str(self.ad_height))
        # print("ad_per_row:" + str(self.ad_per_row))
        # print("padding:" + str(self.padding))

        ad_y = 20

        for index,ad_tree in enumerate(self.ad_list):
            # Set ad_tree position and size
                # Check if line break
            if index % self.ad_per_row == 0 and index != 0:
                ad_y += DEFAULT_AD_HEIGHT * self.scale + self.padding

                # Calculate AD_Tree position in row
            ad_index = index % self.ad_per_row
            ad_x = self.padding + (ad_index * self.ad_width) + (ad_index * self.padding)

            ad_tree.set_position_size(ad_x,ad_y)

            # Calculate block size in each ad_tree
            indv_block_width = DEFAULT_INDV_BLOCK_WIDTH * self.scale
            indv_block_height = DEFAULT_INDV_BLOCK_HEIGHT * self.scale

            indv_block_list,subtree_block_list = ad_tree.individual_subtree_block_list()

            for indv_block in indv_block_list:
                indv_block.width = indv_block_width
                indv_block.height = indv_block_height

            allocatable_space = ad_tree.allocatable_ad_space()
            taxa_cnt = ad_tree.ad_taxa_total()


            for subtree_block in subtree_block_list:
                subtree_block.height = self.block_minimum_height + ( len(subtree_block.belong_subtree.leaf_set) / taxa_cnt * allocatable_space)

    # Before draw_ad, must confirm these variable has been altered according scale
    # canvas.padding(padding might change when ad_per_row change)
    # AD_HEIGHT / AD_WIDTH / ad_per_row /
    # Some variable remain the same
    # ad_tree.padding / ad_tree.initial_x,initial_y /
    # Function to draw ad_tree in each AD
    # Statement in this function will only involve one ad_tree
    def draw_ad_tree(self,ad_tree=None, ad_node=None):

        if ad_node == None:
            ad_node = ad_tree.root

        for node in ad_tree.traverse_postorder():
            if node.type == LEAF:
                block = node.node_or_block

                # Calculate block's position and
                x = ad_tree.topL.x + ad_tree.x  + (ad_tree.x * node.x_level)  # ad_pos + padding
                y = ad_tree.topL.y + ad_tree.y

                block.topL = Point(x, y)
                block.botR = Point(ad_tree.botR.x - ad_tree.x, y + block.height)


                # Calculate branch position
                node.construct_branch(x,y + (block.height / 2))

                # self.branch_head = Point(self.topL.x - AD_BRANCH_LENGTH, self.topL.y + (self.height / 2))
                # self.branch_tail = Point(self.topL.x, self.topL.y + (self.height / 2))

                # Draw block
                if block.type == INDIVIDUAL_BLOCK:
                    block.width = DEFAULT_INDV_BLOCK_WIDTH * self.scale
                    if block.topL.x + block.width >= ad_tree.botR.x:
                        block.botR.x = ad_tree.botR.x - 5
                        block.width = block.botR.x - block.topL.x

                    self.draw_frame(block.topL.x, block.topL.y, block.width, block.height,
                                  GREY,
                                  layer_index=self.AD_LAYER)

                elif block.type == INDIVIDUAL_LEAF_BLOCK:
                    block.width = DEFAULT_INDV_BLOCK_WIDTH * self.scale
                    if block.topL.x + block.width >= ad_tree.botR.x:
                        block.botR.x = ad_tree.botR.x - 5
                        block.width = block.botR.x - block.topL.x

                    self.draw_exact_match_block(block.topL.x, block.topL.y, block.width, block.height,
                                  block_color=block.color,frame_color=GREY,
                                  layer_index=self.AD_LAYER)

                elif block.type == SUBTREE_BLOCK:
                    block.width = block.botR.x - block.topL.x

                    # Draw subtree block
                    if block.exact_match:
                        self.draw_exact_match_block(block.topL.x, block.topL.y, block.width, block.height,
                                      block_color=block.color,frame_color=BLACK,
                                      layer_index=self.AD_LAYER)
                    else:
                        self.draw_inexact_match_block(block.topL.x, block.topL.y, block.width, block.height,
                                                    block_color=block.color,frame_color=BLACK,
                                                    layer_index=self.AD_LAYER)

                    self.write_ad_block_label(block)


                ad_tree.y += block.height + ad_tree.padding

            if node.type == INTERNAL or node.type == ROOT:
                first_child = node.children[0]
                last_child = node.children[-1]

                # Check whether branch length of two children is different
                x = first_child.branch_head.x if first_child.branch_head.x < last_child.branch_head.x else last_child.branch_head.x
                node.unify_children_branches(x)
                # first_child.branch_head.x = x
                # last_child.branch_head.x = x

                y = first_child.branch_head.y + ((last_child.branch_head.y - first_child.branch_head.y) / 2)

                node.construct_branch(x, y)

                # Draw vertical branch linking two children
                self.draw_solid_line(Point(x, first_child.branch_head.y), Point(x, last_child.branch_head.y), BLACK,
                                     self.AD_LAYER)

                if node.type == ROOT:
                    # Draw horizontal branch
                    self.draw_solid_line(node.branch_head, node.branch_tail, BLACK, self.AD_LAYER)

                # Draw children's branches
                for child in node.children:
                    self.draw_solid_line(child.branch_head, child.branch_tail, BLACK, self.AD_LAYER)
    def write_ad_block_label(self,block):
        subtree_label = block.belong_subtree.label
        taxa_cnt = str(block.taxa_count)

        label_width = (len(taxa_cnt) * 10) * self.scale

        self[self.AD_LAYER].font = f'{13 * self.scale}px Times New Romans'
        # Write subtree label on top left corner
        self[self.AD_LAYER].fill_style = BLACK
        self[self.AD_LAYER].fill_text(subtree_label, block.topL.x + 5, block.topL.y + (12 * self.scale))

        # Write taxa count on top right corner
        self[self.AD_LAYER].fill_text(taxa_cnt, block.botR.x - label_width, block.topL.y + (12 * self.scale))

    def individual_block(self,subtree,type=INDIVIDUAL_BLOCK):
        #  Default size : 30x7
        #  Default : topL=Point(x,y),botR=Point(x + 30,y + 7)

        block =  AD_Block(width=DEFAULT_INDV_BLOCK_WIDTH * self.scale, height=DEFAULT_INDV_BLOCK_HEIGHT * self.scale, subtree=subtree, ad_tree=self.ad_constructing, type=type)
        return block


    def ad_subtree_block(self,subtree):
        # AD_Block height depends on node amount, width depends on x_level and tree_width
        block =  AD_Block(width=None,height=None,subtree=subtree,ad_tree=self.ad_constructing,type=SUBTREE_BLOCK)
        return block

    def draw_exact_match_block(self,x,y,width,height,block_color,frame_color,layer_index = -1):
        self.draw_rec(x,y,width,height,block_color,layer_index=layer_index)
        self.draw_frame(x, y, width, height, frame_color, layer_index=layer_index)

    def draw_inexact_match_block(self,x,y,width,height,block_color,frame_color,layer_index = -1):
        self.draw_rec(x, y, width, height, block_color, layer_index=layer_index)
        self.draw_dotted_line(Point(x,y),Point(x+width,y),frame_color,line_dash_index=1, layer_index=layer_index)
        self.draw_dotted_line(Point(x+width,y), Point(x + width, y+height), frame_color,line_dash_index=1, layer_index=layer_index)
        self.draw_dotted_line(Point(x, y), Point(x, y+height), frame_color, line_dash_index=1,layer_index=layer_index)
        self.draw_dotted_line(Point(x, y+height), Point(x + width, y + height), frame_color,line_dash_index=1, layer_index=layer_index)


    # For testing
    def draw_ad_tree_rec(self):
        for index,ad_tree in enumerate(self.ad_list):
            if self.max_ad:
                if index >= self.max_ad:
                    break
            self.draw_frame(ad_tree.topL.x,ad_tree.topL.y,ad_tree.width,ad_tree.height,BLACK,self.AD_LAYER)




'''
Full canvas width = 1100
AD width = 135
AD height = 150
x-gap/y-gap between AD = 20
'''


class Point():
    def __init__(self,x,y):
        self.x = x
        self.y = y

class Block():
    def __init__(self,topL,botR):
        self.topL = topL
        self.botR = botR
        self.branch_x = None

        if topL and botR:
            self.calculate_width_height()

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

    def check_nested_block(self,block):
        return self.check_in_range(block.topL) and block.botR.y <= self.botR.y

# Block for ad_tree
class AD_Block(Block):
    # id = None  # Block id in ad_tree

    def __init__(self, width,height,subtree,ad_tree,type,topL=None, botR=None):
        super().__init__(topL, botR)
        self.root = None  # To record this block represent which branch. If represent a subtree,root is subtree root
        self.taxa_list = set()  # Optional, taxa_list may ignore if has root(root is node from tc_tree,thus can get leaf_node via root)
        self.taxa_count = 0
        self.color = None
        self.belong_subtree = subtree
        self.belong_ad_tree = ad_tree

        self.type = type
        self.width = width
        self.height = height

    def construct_branch(self):
        self.branch_head = Point(self.topL.x - (AD_BRANCH_LENGTH * self.scale),self.topL.y + (self.height/2))
        self.branch_tail = Point(self.topL.x ,self.topL.y + (self.height/2))










