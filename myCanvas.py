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
    
    def __init__(self,tree,width, height,view_support):
        super().__init__(width = width, height = height)
        self.tree = tree
        self.view_support = view_support
        self.draw_rt(node = self.tree.seed_node)

    def draw_rt(self,node=None,level=0):
        if node is None:
            node = self.tree.seed_node

        for child in node.child_node_iter():
            self.draw_rt(child, level + 1)

        self.fill_style = 'black'
        self.font = LEAF_FONT
        if node.is_leaf():
            if self.view_support:
                node.x = self.x * level
            else:
               node.x = X_INTERVAL * level + X_INTERVAL

            node.y = self.y
            self.begin_path()
            self.move_to(node.x, node.y - 5)
            self.line_to(node.x + X_INTERVAL, node.y - 5)
            self.stroke()
            self.fill_text(node.taxon.label, node.x + X_INTERVAL,node.y)
            self.y += Y_INTERVAL

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
                self.fill_style = 'black'

            # # testing
            # self.draw_dots(node.x,node.y)

    def draw_dots(self,x,y):
        self.fill_style = 'red'
        self.fill_circle(x, y, 1)
        self.fill_style = 'black'

    def draw_line(self,x1,x2,y1,y2):
        self.fill_style = 'red'
        self.begin_path()
        self.move_to(x1, y1)
        self.line_to(x2,y2)
        self.stroke()

        self.fill_style = 'black'

    def draw_rec(self,x,y,width,height):
        self.fill_style = BLUE
        self.stroke_style = BLACK

        # fill_rect(x, y, width, height=None):
        # stroke_rect(x, y, width, height=None):
        self.fill_rect(x, y, width, height)
        self.stroke_rect(x, y, width, height)


class tcCanvas(MyCanvas): # Canvas for Tree Collection
    id = None
    def __init__(self,tree,width=800, height=600):
        super().__init__(width = width, height = height)
        self.tree = tree
        self.id = tree.id
