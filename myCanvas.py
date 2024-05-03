from ipycanvas import Canvas,MultiCanvas
import ADpy as myTree
from Utils import *
import math

# Use Multicanvas which has multiple canvas layer
class MyCanvas(MultiCanvas):
    tree = None
    line_dashes = [[2, 2], [5, 5], [8, 8], [10, 10]]

    def __init__(self,layer,width, height):
        super().__init__(layer,width = width, height = height, sync_image_data=True)
        self.layer = layer

    def draw_dots(self,x,y,layer_index = -1,canvas=None):
        if not canvas:
            canvas = self[layer_index]
        canvas.fill_style = RED
        canvas.fill_circle(x, y, 1)
        canvas.fill_style = BLACK

    def draw_solid_line(self,head,tail,color,layer_index = -1,canvas=None):
        if not canvas:
            canvas = self[layer_index]
        canvas.stroke_style = color
        canvas.begin_path()
        canvas.move_to(head.x, head.y)
        canvas.line_to(tail.x,tail.y)
        canvas.stroke()

    def draw_dotted_line(self,head,tail,color,line_dash_index=0,layer_index = -1,canvas=None):
        if not canvas:
            canvas = self[layer_index]
        canvas.stroke_style = color
        canvas.set_line_dash(self.line_dashes[line_dash_index])
        canvas.begin_path()
        canvas.move_to(head.x, head.y)
        canvas.line_to(tail.x,tail.y)
        canvas.stroke()

        canvas.set_line_dash([])

    # Draw rectangle on specific layer
    def draw_rec(self,x,y,width,height,color,layer_index = -1,canvas=None):
        if not canvas:
            canvas = self[layer_index]
        canvas.fill_style = color
        canvas.fill_rect(x, y, width,height)

    def draw_frame(self,x,y,width,height,color,layer_index = -1,canvas=None):
        if not canvas:
            canvas = self[layer_index]
        canvas.stroke_style = color
        canvas.stroke_rect(x, y, width, height)

    def draw_dotted_frame(self,x,y,width,height,color,line_dash_index=0,layer_index = -1,canvas=None):
        if not canvas:
            canvas = self[layer_index]
        canvas.stroke_style = color
        canvas.set_line_dash(self.line_dashes[line_dash_index])
        canvas.stroke_rect(x, y, width, height)

        canvas.set_line_dash([])

    def draw_elided_branch(self,head,tail,color,layer_index = -1,canvas=None):
        self.draw_dotted_line(head,tail,color,line_dash_index=0,layer_index=layer_index,canvas=canvas)

        x = head.x + (tail.x - head.x) / 2
        y = tail.y - 5
        start_point = Point(x + 1,y)
        end_point = Point(x - 2,y + 10)

        self.draw_solid_line(head=start_point,tail=end_point,color=color,layer_index=layer_index,canvas=canvas)

    def crop_and_paste_canvas(self, paste_canvas,layer_index,copy_list,paste_list):
        x = copy_list[0]
        y = copy_list[1]
        width = copy_list[2]
        height = copy_list[3]

        crop = self[layer_index].get_image_data(x,y,width,height)

        x = paste_list[0]
        y = paste_list[1]

        paste_canvas.put_image_data(crop, x,y)
        # copy_canvas.observe(crop_and_paste_canvas, "image_data")

    def get_subtree(self,subtree_list,label):
        for subtree in subtree_list:
            if subtree.label == label:
                return subtree

        return None


class Point():
    def __init__(self,x,y):
        self.x = x
        self.y = y

class Block():
    def __init__(self,topL,botR):
        self.topL = topL
        self.botR = botR
        self.branch_x = None
        self.segment_list = []

        if topL and botR:
            self.calculate_width_height()
            self.topR = Point(self.botR.x, self.topL.y)
            self.botL = Point(self.topL.x, self.botR.y)

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

    def check_nested_block(self,block,select_tree=None):
        if not select_tree or select_tree == RT:
            return self.check_in_range(block.topL) and block.botR.y <= self.botR.y
        elif select_tree == TC:
            return self.check_in_range(block.botR) and block.botR.y <= self.botR.y

# Block for ad_tree
class AD_Block(Block):
    # id = None  # Block id in ad_tree

    def __init__(self, width,height,subtree,ad_tree,type,topL=None, botR=None):
        super().__init__(topL, botR)
        self.root = None  # To record this block represent which branch. If represent a subtree,root is subtree root
        self.taxa_list = set()
        self.subtree_taxa_count = 0   # subtree_taxa in block
        self.duplicate_subtree_taxa_count = []  # taxa_count of each subtree in same block
        self.exact_match = None
        self.color = None
        self.belong_subtree = subtree
        self.belong_ad_tree = ad_tree
        self.nested_tree = None
        self.is_elided = False
        self.is_subtree_duplicate = False

        self.type = type
        self.width = width
        self.height = height


class Missing_Taxa_Segment(Block):
    def __init__(self, subtree, ad_tree, missing_taxa,topL=None, botR=None):
        super().__init__(topL, botR)
        self.ad_tree = ad_tree
        self.subtree = subtree # Missing taxa from which subtree
        self.missing_taxa_list = missing_taxa # Include which taxa

class Subtree_Block_Segment(Block):
    def __init__(self, subtree,belong_block,topL=None, botR=None):
        super().__init__(topL, botR)
        self.subtree = subtree
        self.belong_block = belong_block
        self.belong_block.segment_list.append(self)




def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))

def rgb_to_hex(rgb_color):
    return '#{:02x}{:02x}{:02x}'.format(*rgb_color)

def lighten_color(hex_color, amount=10):
    rgb_color = hex_to_rgb(hex_color)
    lightened_rgb = tuple(min(255, c + amount) for c in rgb_color)
    return rgb_to_hex(lightened_rgb)




