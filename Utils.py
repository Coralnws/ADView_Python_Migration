# CONSTANT
LEAF = "Leaf"
INTERNAL = 'Internal Node'
RT_NODE_COLOR = '#F38D76'
X_INTERVAL = 10
Y_INTERVAL = 20
FONT_SIZE = 16
NODE_BLOCK_HEIGHT = 13
INTERNAL_FONT = '13px Times New Romans'
LEAF_FONT = '16px Times New Romans'
CANVAS_MAX_WIDTH = 1000

# COLOR
BLANK = '#FFFFFF'
# RED = '#F38D76'
RED = '#BF726D'
# BLUE = '#33D5FF'
BLUE = '#82A8D5'
PURPLE = '#9A86BC'
GREEN = '#6F9D55'
BEIGE = '#F8DAAB'
BLACK = '#000000'
ORANGE = '#E8A553'

# AD view
INDIVIDUAL = "Individual"
CLUSTER = "Cluster"


# Subtree label
A = 1
B = 2
C = 3
D = 4
E = 5

class Point():
    def __init__(self,x,y):
        self.x = x
        self.y = y

class Block():
    def __init__(self,topL,botR):
        self.topL = topL
        self.botR = botR

    def check_inRange(self,point):
        return point.x >= self.topL.x and point.x <= self.botR.x and point.y >= self.topL.y and point.y <= self.botR.y

    def print_block(self):
        print("topL = (",self.topL.x,",",self.topL.y,")")
        print("botR = (",self.botR.x,",",self.botR.y,")")

    def get_size(self):
        return self.botR.y - self.topL.y

