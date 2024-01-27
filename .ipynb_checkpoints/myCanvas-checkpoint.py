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
    
    def __init__(self,tree,width=800, height=2080):
        super().__init__(width = width, height = height)
        self.tree = tree
        self.draw_rt(node = self.tree.seed_node)

    def draw_nodes(self,node):
        if not node.is_leaf():
            self.fill_style = RT_NODE_COLOR
            self.fill_rect(node.x - 5, node.y - 5, 5, 5)

        self.fill_style = 'black'
        self.font = '13px Times New Romans'
        self.fill_text(node.name, node.x + 5, node.y + 5)

    def draw_rt(self,node=None,level=0):
        if node is None:
            node = self.tree.seed_node

        for child in node.child_node_iter():
            self.draw_rt(child, level + 1)

        self.fill_style = 'black'
        self.font = '16px Times New Roman'
        if node.is_leaf():
            node.x = self.x * level + X_INTERVAL
            node.y = self.y
            self.begin_path()
            self.move_to(node.x - 10, node.y - 5)
            self.line_to(node.x - 5, node.y - 5)
            self.stroke()
            self.fill_text(node.name, node.x,node.y)
            self.y += Y_INTERVAL
        else:
            node.x = self.x * level + X_INTERVAL
            child_count = len(node.children)
            child_one = node.children[0]
            child_last = node.children[child_count-1]
            self.begin_path()
            self.move_to(child_one.x - 10, child_one.y - 5)
            self.line_to(child_one.x - 10, child_last.y - 5)
            self.stroke()

            tmp = child_one.y + child_last.y
            
            node.y = tmp/2

            self.begin_path()
            self.move_to(node.x -10 , node.y - 5)
            self.line_to(self.x * (level + 1), node.y - 5)
            self.stroke()
            
            
            self.fill_style = 'blue'
            self.fill_text(node.support, node.x - 5, node.y - 8)
            self.fill_style = 'black'

class tcCanvas(MyCanvas): # Canvas for Tree Collection
    id = None
    def __init__(self,tree,width=800, height=600):
        super().__init__(width = width, height = height)
        self.tree = tree
        self.id = tree.id
