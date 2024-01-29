'''
Font size - leaf interval:
1. 18px - 20y
2. 16px - 15y
'''
import ipycanvas as ipc
from ipycanvas import Canvas
import myTree as Tree
from Utils import *

class MyCanvas(Canvas):
    '''
    Record this canvas belongs to which tree
    Additional Info:
        1. Subtree chosen
        2.
    '''
    tree = None
    node_name_list = []
    def __init__(self,width, height):
        super().__init__(width = width, height = height)


class rtCanvas(MyCanvas):  # Canvas for Reference Tree
    x = 30
    y = 30
    view_support = False
    subtree_color_list = [BLUE,RED,ORANGE,BEIGE,PURPLE]
    subtree_color_used = [1,1,1,1,1]
    
    def __init__(self,tree,width, height,view_support):
        super().__init__(width = width, height = height)
        self.tree = tree
        self.view_support = view_support
        self.draw_rt(node = self.tree.seed_node)
        self.on_mouse_down(self.mouse_clicked)

    def draw_rt(self,node=None,level=0):
        if node is None:
            node = self.tree.seed_node

        for child in node.child_node_iter():
            self.draw_rt(child, level + 1)

        self.fill_style = BLACK
        self.font = LEAF_FONT
        if node.is_leaf():
            if self.view_support:
                node.x = self.x * level
            else:
               node.x = X_INTERVAL * level + X_INTERVAL

            node.y = self.y
            self.draw_leaf_node(node)
            self.y += Y_INTERVAL

            node.x1 = node.x + X_INTERVAL - 2
            node.x2 = node.x1 + X_INTERVAL + len(str(node.taxon.label)) * 7.5
            node.y1 = node.y - NODE_BLOCK_HEIGHT - 2
            node.y2 = node.y1 + NODE_BLOCK_HEIGHT + 6

            # testing
            # self.draw_dots(node.x + X_INTERVAL - 5 , node.y - NODE_BLOCK_HEIGH  T)
            # self.draw_dots(node.x + X_INTERVAL + len(str(node.taxon.label)) * 9, node.y + 2)
            # self.draw_rec(node.x + X_INTERVAL - 2,node.y - NODE_BLOCK_HEIGHT,X_INTERVAL + len(str(node.taxon.label)) * 9,NODE_BLOCK_HEIGHT+2)

        else:
            if self.view_support:
                node.x = self.x * level
            else:
                node.x = X_INTERVAL * level + X_INTERVAL
            child_nodes = node.child_nodes()
            first_child = child_nodes[0] if child_nodes else None
            last_child = child_nodes[-1] if child_nodes else None
            self.begin_path()
            self.move_to(first_child.x, first_child.y - 5)
            self.line_to(last_child.x, last_child.y - 5)
            self.stroke()

            tmp = first_child.y + last_child.y
            
            node.y = tmp/2

            self.begin_path()
            self.move_to(node.x , node.y - 5)
            if self.view_support:
                self.line_to(self.x * (level+1), node.y - 5)
            else:
                self.line_to(X_INTERVAL * (level + 1) + X_INTERVAL, node.y - 5)
            self.stroke()

            if self.view_support:
                self.fill_style = 'blue'
                if node.label:
                    self.fill_text(node.label, node.x + 3, node.y - 8)
                self.fill_style = BLACK

            node.x1 = node.x + X_INTERVAL - 4
            node.x2 = node.x1 + X_INTERVAL
            node.y1 = node.y - X_INTERVAL
            node.y2 = node.y1 + X_INTERVAL


            # # testing
            # self.draw_dots(node.x + X_INTERVAL - 4, node.y - X_INTERVAL)
            # self.draw_dots(node.x + X_INTERVAL + 4, node.y )
            # self.draw_rec(node.x + X_INTERVAL - 4,node.y - X_INTERVAL,X_INTERVAL,X_INTERVAL)

    def draw_leaf_node(self,node):
        self.fill_style = BLACK
        self.begin_path()
        self.move_to(node.x, node.y - 5)
        self.line_to(node.x + X_INTERVAL, node.y - 5)
        self.stroke()
        self.fill_text(node.taxon.label, node.x + X_INTERVAL, node.y)


    def filter_node_selected(self,node, x , y):
        if node.x1 <= x <= node.x2 and node.y1 <= y <= node.y2:
            return True
        return False

    def mouse_clicked(self, x, y):
        # self.draw_dots(x,y)
        node_selected = self.tree.find_node(lambda node: self.filter_node_selected(node, x,y))


        if node_selected.is_leaf():
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
        self.fill_rect(x, y, width, height)

    def select_leaf_node(self,node_selected):
        if not hasattr(node_selected, 'selected') or node_selected.selected == False:
            color_index = self.subtree_color_used.index(1)
            self.subtree_color_used[color_index] = 0
            color = self.subtree_color_list[color_index]
            node_selected.color = color
            self.draw_rec(node_selected.x1, node_selected.y1, node_selected.x2 - node_selected.x1,
                          node_selected.y2 - node_selected.y1, color)
            self.draw_leaf_node(node_selected)
        else:
            self.draw_rec(node_selected.x1, node_selected.y1, (node_selected.x2 - node_selected.x1) + 2,
                          node_selected.y2 - node_selected.y1, BLANK)
            self.draw_leaf_node(node_selected)
            color_index = self.subtree_color_list.index(node_selected.color)
            self.subtree_color_used[color_index] = 1



class tcCanvas(MyCanvas): # Canvas for Tree Collection
    id = None
    def __init__(self,tree,width=800, height=600):
        super().__init__(width = width, height = height)
        self.tree = tree
        self.id = tree.id
