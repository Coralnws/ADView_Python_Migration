import ipycanvas as ipc
from ipycanvas import Canvas
import myTree as Tree
from Utils import *

class Point():
    def __init__(self,x,y):
        self.x = x
        self.y = y

class MyCanvas(Canvas):
    '''
    Additional Info:
        1. Subtree chosen
        2.
    '''
    tree = None
    node_name_list = []
    def __init__(self,width, height):
        super().__init__(width = width, height = height)

# Canvas for Reference Tree
class rtCanvas(MyCanvas):
    x = 30
    y = 30

    # Show support value
    view_support = False

    # Block color
    subtree_color_list = [BLUE,RED,ORANGE,BEIGE,PURPLE]
    subtree_color_used = [1,1,1,1,1]
    tree_initialize = 1

    
    def __init__(self,tree,width, height,view_support):
        super().__init__(width = width, height = height)
        self.tree = tree
        self.view_support = view_support
        self.tree_width = 1
        self.draw_rt(node = self.tree.seed_node)
        self.on_mouse_down(self.mouse_clicked)
        self.tree_initialize = 0

    # Draw reference tree on rtcanvas
    def draw_rt(self,node=None,level=0):
        # Default : Tree root as reference tree's root
        if node is None:
            node = self.tree.seed_node

        for child in node.child_node_iter():
            self.draw_rt(child, level + 1)

        # Set canvas style
        self.fill_style = BLACK
        self.font = LEAF_FONT

        if self.tree_initialize:
            node.level = level

        if node.is_leaf():
            if self.tree_initialize:
                if self.view_support:
                    x = self.x * node.level
                else:
                    x = X_INTERVAL * node.level + X_INTERVAL

                node.pos = Point(x,self.y)  # Node position

                self.y += Y_INTERVAL

                # Mouse click range of current leaf node
                node.range_topL = Point(node.pos.x, node.pos.y - NODE_BLOCK_HEIGHT )
                node.range_botR = Point(node.range_topL.x + X_INTERVAL + len(str(node.taxon.label)) * 12,
                                        node.range_topL.y + NODE_BLOCK_HEIGHT + 3)

                node.block_top = node.range_topL.y
                node.block_bottom = node.range_botR.y

                if node.range_botR.x > self.tree_width:
                    self.tree_width = node.range_botR.x + X_INTERVAL


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

            if self.tree_initialize:
                if self.view_support:
                    x = self.x * node.level
                else:
                    x = X_INTERVAL * node.level + X_INTERVAL

                tmp = first_child.pos.y + last_child.pos.y
                node.pos = Point(x, tmp / 2)

                # Mouse click range of current internal node
                node.range_topL = Point(node.pos.x, first_child.pos.y)
                node.range_botR = Point(node.range_topL.x + X_INTERVAL, last_child.pos.y - X_INTERVAL)


            # Draw vertical branch
            self.begin_path()
            self.move_to(first_child.pos.x, first_child.pos.y - 5)
            self.line_to(last_child.pos.x, last_child.pos.y - 5)
            self.stroke()

            # Drawing horizontal branch
            self.begin_path()
            self.move_to(node.pos.x , node.pos.y - 5)
            if self.view_support:
                self.line_to(self.x * (node.level+1), node.pos.y - 5)
            else:
                self.line_to(X_INTERVAL * (node.level + 1) + X_INTERVAL, node.pos.y - 5)
            self.stroke()

            if self.view_support:
                self.fill_style = 'blue'
                if node.label:
                    self.fill_text(node.label, node.pos.x + 3, node.pos.y - 8)
                self.fill_style = BLACK

            # Testing
            # self.draw_dots(node.pos.x + X_INTERVAL - 4, node.pos.y - X_INTERVAL)
            # self.draw_dots(node.pos.x + X_INTERVAL + 4, node.pos.y )
            # self.draw_rec(node.range_topL.x,node.range_topL.y,8,node.range_botR.y - node.range_topL.y,BEIGE)

        # To calculate subtree's block size
        if self.tree_initialize and node.parent_node:
            if hasattr(node.parent_node, 'block_top'):
                if node.block_top < node.parent_node.block_top:
                    node.parent_node.block_top = node.block_top
            else:
                node.parent_node.block_top = node.block_top

            if hasattr(node.parent_node, 'block_bottom'):
                if node.block_bottom > node.parent_node.block_bottom:
                    node.parent_node.block_bottom = node.block_bottom
            else:
                node.parent_node.block_bottom = node.block_bottom

    # Drawing leaf node on canvas
    def draw_leaf_node(self,node):
        self.fill_style = BLACK
        self.begin_path()
        self.move_to(node.pos.x, node.pos.y - 5)
        self.line_to(node.pos.x + X_INTERVAL, node.pos.y - 5)
        self.stroke()
        self.fill_text(node.taxon.label, node.pos.x + X_INTERVAL, node.pos.y)

    # Filter function
    def filter_node_selected(self,node, x , y):
        if node.range_topL.x <= x <= node.range_botR.x and node.range_topL.y <= y <= node.range_botR.y:
            return True
        return False

    # Handling mouse click events in rt canvas
    def mouse_clicked(self, x, y):
        # self.draw_dots(x,y)
        node_selected = self.tree.find_node(lambda node: self.filter_node_selected(node, x,y))

        self.select_leaf_node(node_selected)
        if hasattr(node_selected, 'selected'):
            node_selected.selected = not node_selected.selected
        else:
            node_selected.selected = True


    def draw_dots(self,x,y):
        self.fill_style = RED
        self.fill_circle(x, y, 1)
        self.fill_style = BLACK

    def draw_line(self,x1,x2,y1,y2):
        self.fill_style = RED
        self.begin_path()
        self.move_to(x1, y1)
        self.line_to(x2,y2)
        self.stroke()

        self.fill_style = BLACK

    def draw_rec(self,x,y,width,height,color):
        self.fill_style = color
        # fill_rect(x, y, width, height=None):
        # stroke_rect(x, y, width, height=None):
        self.fill_rect(x, y, width,height)

    # If a leaf node being clicked
    def select_leaf_node(self,node_selected):

        if node_selected.is_leaf():
            rec_x = node_selected.range_topL.x + X_INTERVAL - 5
            rec_y = node_selected.range_topL.y
            rec_width = self.tree_width - rec_x
            rec_height = node_selected.range_botR.y - node_selected.range_topL.y

        else:
            rec_x = node_selected.range_topL.x + 5
            rec_y = node_selected.block_top
            rec_width = self.tree_width - rec_x
            rec_height = node_selected.block_bottom - rec_y
            # node_selected.range_botR.y - node_selected.range_topL.y + X_INTERVAL


        if not hasattr(node_selected, 'selected') or node_selected.selected == False:
            # If the node was not selected → draw block and set node as selected_subtree
            color_index = self.subtree_color_used.index(1)
            self.subtree_color_used[color_index] = 0
            color = self.subtree_color_list[color_index]
            node_selected.color = color
            # width  = node_selected.range_botR.x - node_selected.range_topL.x
            self.draw_rec(rec_x,rec_y,rec_width,rec_height,color)
            self.draw_rt(node = node_selected)
        else:
            self.draw_rec(rec_x , rec_y, rec_width + 10, rec_height, BLANK)
            self.draw_rt(node=node_selected)

            color_index = self.subtree_color_list.index(node_selected.color)
            self.subtree_color_used[color_index] = 1


        # # If the node was not selected → draw block and set node as selected_subtree
        # if not hasattr(node_selected, 'selected') or node_selected.selected == False:
        #     color_index = self.subtree_color_used.index(1)
        #     self.subtree_color_used[color_index] = 0
        #     color = self.subtree_color_list[color_index]
        #     node_selected.color = color
        #     # width  = node_selected.range_botR.x - node_selected.range_topL.x
        #     self.draw_rec(node_selected.range_topL.x + X_INTERVAL - 5 , node_selected.range_topL.y, self.tree_width - node_selected.range_topL.x,
        #                   node_selected.range_botR.y - node_selected.range_topL.y, color)
        #     self.draw_leaf_node(node_selected)
        #
        # else:
        #     # If the node was selected → erase block and remove node from selected_subtree
        #     # width = (node_selected.range_botR.x - node_selected.range_topL.x) + 2
        #     self.draw_rec(node_selected.range_topL.x + X_INTERVAL - 9, node_selected.range_topL.y, self.tree_width - node_selected.range_topL.x + 9,
        #                   node_selected.range_botR.y - node_selected.range_topL.y, BLANK)
        #     self.draw_leaf_node(node_selected)
        #     color_index = self.subtree_color_list.index(node_selected.color)
        #     self.subtree_color_used[color_index] = 1


# Canvas for Tree Collection
class tcCanvas(MyCanvas):
    id = None
    def __init__(self,tree,width=800, height=600):
        super().__init__(width = width, height = height)
        self.tree = tree
        self.id = tree.id
