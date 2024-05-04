import math
from ipycanvas import Canvas, hold_canvas
import ipywidgets as widgets
from IPython.display import display
from .myCanvas import *
from .myUtils import *

class TreeDistributionView:
    # Size
    subtree_button_width = 120
    subtree_button_height = 40

    subtree_button_x = 50
    subtree_button_y = 50

    spacer_height = 10
    spacer_width = 10
    button_gap = 20

    segment_bar_width = 750

    DEFAULT_CANVAS_WIDTH = 700 # or 700

    '''
    Canvas's height depends on its content
       - node_cnt * label_height
       - tree_cnt * label_height
       - ad_tree_height
    '''

    FONT_FAMILY = 'Times New Roman'
    FONT_SIZE = '16px'
    CHOSEN_FONT_SIZE = '20px'

    def __init__(self,adPy,tc_canvas_tmp,test=False):
        self.adPy = adPy
        self.tc_canvas_tmp = tc_canvas_tmp

        # Subtree button
        self.subtree_button_list = []
        self.subtree_button_style = []

        # Segment button
        self.segment_button_list = []
        self.segment_button_style = []

        # Nodes Canvas
        self.nodes_list_canvas = None

        # Tree id / Name canvas
        self.related_tc_tree_canvas = None

        # Cluster canvas
        self.cluster_canvas = None
        self.differentiate_inexact_match = True

        # Subtree block canvas
        self.rt_subtree_block_canvas = None

        # Layout
        self.subtree_button_hbox = None
        self.segment_button_hbox = None
        self.canvas_vbox = None

        # Vbox and Hbox used
        self.base_vbox = None
        self.canvas_hbox = None


        # Checkbox
        self.generate_checkbox_vbox()

        # Other Info
        self.tc_tree_list = None
        self.nodes_list = None
        self.cluster_list = []
        self.show_tree_name = False
        self.test = test
        self.empty_canvas = Canvas(width = self.DEFAULT_CANVAS_WIDTH,height = 50)
        self.canvas_vbox_height = 0
        self.agree_rt = False
        self.export_ready = False

        # Setup spacer
        self.vspacer = widgets.Output(layout=widgets.Layout(height=f'{self.spacer_height}px'))
        self.hspacer = widgets.Output(layout=widgets.Layout(width=f'{self.spacer_width}px'))
        self.button_spacer = widgets.Output(layout=widgets.Layout(width=f'{self.button_gap}px'))


        # Manage button event
        self.deselect_button = None
        self.subtree_button_clicked = None
        self.segment_button_clicked = None
        self.subtree_chosen = None

        # Initialization
        if test == True:
            self.testing_generate_subtree()
        else:
            self.subtree_list = self.adPy.subtree_list

        self.setup_subtree_button()


    # TESTING
    def testing_generate_subtree(self):
        self.subtree_list = []
        for i in range(5):
            subtree = Subtree(SUBTREE_LABEL_LIST[i], None, None, SUBTREE_COLOR_LIST[i])
            subtree.topology_list = {'Chlorokybus atmophyticus,Mesostigma viride,Spirotaenia minuta': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 63, 64, 65, 66, 67, 68], 'Chlorokybus atmophyticus,Spirotaenia minuta': [61]}
            self.subtree_list.append(subtree)

        # print("end testing function")


    # --------------------- Generate element / Setup Layout ----------------------
    def generate_checkbox_vbox(self):
        self.tree_name_checkbox = widgets.Checkbox(description='Show Tree Name',layout=widgets.Layout(width = 'auto'))
        self.tree_name_checkbox.observe(self.on_treename_checkbox_change, names='value')
        self.cluster_checkbox = widgets.Checkbox(description='Differentiate Inexact Match in Cluster',value=True,
                                                 layout=widgets.Layout(width = 'auto'))
        self.cluster_checkbox.observe(self.on_cluster_checkbox_change, names='value')

        self.checkbox_vbox = widgets.VBox(children=[self.tree_name_checkbox, self.cluster_checkbox],
                                          layout=widgets.Layout(width='auto'))

    def generate_subtree_button(self):
        # self.subtree_list = sorted(self.adPy.subtree_list, key=lambda x: x.label,reverse=False)
        # self.testing_generate_subtree()
        # self.subtree_list = self.adPy.subtree_list

        self.subtree_list = sorted(self.subtree_list, key=lambda x:x.label)
        for subtree in self.subtree_list:
            new_button_style = widgets.ButtonStyle(button_color=subtree.color)

            self.subtree_button_style.append(new_button_style)

            new_button = widgets.Button(description=subtree.label, style=new_button_style,
                                  layout=widgets.Layout(width=f'{self.subtree_button_width}px', height=f'{self.subtree_button_height}px'))
            new_button.style.font_weight = 'normal'
            new_button.style.font_family = self.FONT_FAMILY
            new_button.style.font_size = self.FONT_SIZE


            new_button.on_click(self.on_subtree_button_click)
            new_button.metadata = {'subtree': subtree}

            self.subtree_button_list.append(new_button)

    def generate_subtree_button_hgridbox(self,button_list):
        children = []
        for button in button_list:
            children.append(button)
            children.append(self.button_spacer)

        children.pop()

        tmp_hbox = widgets.HBox(children=children)

        title_label = widgets.Label(value='Tree Distribution by Reference Partition')
        title_label.style.font_weight = 'bold'
        title_label.style.font_family = self.FONT_FAMILY
        title_label.style.font_size = '18px'

        subtree_button_label = widgets.Label(value='Partition')
        subtree_button_label.style.font_weight = 'bold'
        subtree_button_label.style.font_family = self.FONT_FAMILY
        subtree_button_label.style.font_size = '18px'

        tmp_vbox = widgets.VBox(children=[title_label,subtree_button_label,tmp_hbox])

        children = [tmp_vbox]
        children.append(self.checkbox_vbox)

        # layout = widgets.Layout(grid_template_columns=self.DEFAULT_GRID_TEMPLATE,grid_gap=f'{self.DEFAULT_GRIP_GAP}px')
        self.subtree_button_hbox = widgets.HBox(children=children)

    def generate_segment_button_hbox(self):
        # Testing
        if self.adPy.tc:
            total_grid = len(self.adPy.tc)
        else:
            total_grid = 68

        segment_count = len(self.subtree_chosen.topology_list.keys())
        segment_allocatable_width = self.segment_bar_width - segment_count * MINIMUN_BUTTON_WIDTH
        color = LIGHT_GREY

        sorted_topology_list = sorted(self.subtree_chosen.topology_list.items(), key=lambda x: len(x[1][1]),
                                      reverse=True)
        self.segment_button_list = []
        for item in sorted_topology_list:
            tree_list = item[1][1]

            new_button_style = widgets.ButtonStyle(button_color=color)
            color = lighten_color(color,amount=20)

            button_width = MINIMUN_BUTTON_WIDTH + (len(tree_list)/total_grid * segment_allocatable_width)
            new_button = widgets.Button(description=str(len(tree_list)), style=new_button_style,
                                  layout=widgets.Layout(width=f'{button_width}px', height=f'{self.subtree_button_height}px',justify_content='flex-start'))
            new_button.style.font_weight = 'normal'
            new_button.style.font_family = self.FONT_FAMILY
            new_button.style.font_size = self.FONT_SIZE
            new_button.metadata = {'nodes_str' : item[0],'nodes_list': item[1][0],'tree_list': tree_list}
            new_button.on_click(self.on_segment_button_click)

            self.segment_button_list.append(new_button)
            self.segment_button_hbox = widgets.HBox(children=self.segment_button_list,layout=widgets.Layout(justify_content = 'flex-start'))

        subtree_button_label = widgets.Label(value='Distribution')
        subtree_button_label.style.font_weight = 'bold'
        subtree_button_label.style.font_family = self.FONT_FAMILY
        subtree_button_label.style.font_size = '18px'

        tmp_vbox = widgets.VBox(children=[subtree_button_label, self.segment_button_hbox])

        self.base_vbox.children = list(self.base_vbox.children) + [tmp_vbox]

        # display(self.segment_button_hbox)


        # for subtree in self.subtree_list:
        #     print(f"Subtree {subtree.label}")
        #     index = 1
        #     print(subtree.topology_list)
        #     for topology,tree_list in subtree.topology_list.items():
        #         print(f"Topology {index} : {len(tree_list)}")
        #         print(topology)
        #         index += 1

        # button_a_2.layout.grid_area = '1 / 1 / 1 / 4'
        # segment_list {'button_index_in_list' : col_num }

    def generate_canvas_vgridbox(self,canvas_list):
        # canvas_list.insert(0,self.tree_name_checkbox)
        children = []
        for canvas in canvas_list:
            children.append(canvas)
            children.append(self.vspacer)

        # children.append(self.empty_canvas)

        self.canvas_vbox = widgets.VBox(children=children)

        # display(self.canvas_vgridbox)

    # ---------------------- Draw element -------------------------
    def draw_subtree_button(self):
        self.setup_subtree_button_layout(self.subtree_button_list)

    def draw_base_vbox(self,layout_list):
        self.base_vbox = widgets.VBox(children=layout_list)
        display(self.base_vbox)

    def draw_canvas_hbox(self,layout_list):
        children = []
        for layout in layout_list:
            children.append(layout)

        self.canvas_hbox = widgets.HBox(children=children)
        self.base_vbox.children = list(self.base_vbox.children) + [self.vspacer,self.canvas_hbox]

    def draw_canvas_vgridbox(self):
    #     grid = widgets.GridBox([canvas1, canvas2, canvas3], layout=widgets.Layout(grid_template_columns="1fr", grid_gap='0px'))
        pass

    def draw_nodes_canvas(self):
        nodes_list = self.segment_button_clicked.metadata.get('nodes_str').split(',')

        tmp_canvas = Canvas(width=self.DEFAULT_CANVAS_WIDTH, height=len(nodes_list) * RT_Y_INTERVAL + 30)
        tmp_canvas.fill_style = BLACK
        tmp_canvas.font = f'18px {self.FONT_FAMILY}'

        pointer_x = 10
        pointer_y = 20
        tmp_canvas.fill_text("Nodes in Subtree: ", pointer_x, pointer_y)
        pointer_y += 5

        tmp_canvas.font = f'16px {self.FONT_FAMILY}'
        pointer_y += RT_Y_INTERVAL

        for node in nodes_list:
            tmp_canvas.fill_text(node, pointer_x, pointer_y)
            tmp_canvas.fill_text(node, pointer_x, pointer_y)
            pointer_y += RT_Y_INTERVAL

        self.nodes_list_canvas.height = pointer_y
        self.nodes_list_canvas.draw_image(tmp_canvas, 0, 0)

        self.canvas_vbox_height += pointer_y

        self.draw_canvas_frame(self.nodes_list_canvas)

    def draw_tc_tree_canvas(self):
        height = RT_Y_INTERVAL * len(self.tc_tree_list) + 30
        tmp_canvas = Canvas(width=self.DEFAULT_CANVAS_WIDTH, height=height)
        tmp_canvas.fill_style = BLACK
        tmp_canvas.font = f'18px {self.FONT_FAMILY}'

        pointer_x = 10
        pointer_y = 20
        tmp_canvas.fill_text("Tree Colletion:", pointer_x, pointer_y)
        pointer_y += 5

        tmp_canvas.font = f'16px {self.FONT_FAMILY}'
        pointer_y += RT_Y_INTERVAL
        str = ""
        for index,tc_tree_id in enumerate(self.tc_tree_list):
            str_drew = False
            if self.show_tree_name:
                if self.test:
                    str = f"{tc_tree_id} : Chlorokybus atmophyticus"
                else:
                    str = f"{tc_tree_id} : {self.adPy.tc[tc_tree_id-1].name}"

                tmp_canvas.fill_text(str, pointer_x, pointer_y)
                tmp_canvas.fill_text(str, pointer_x, pointer_y)

                pointer_y += RT_Y_INTERVAL
            else:
                str += f"{tc_tree_id}, "
                if index > 0 and index % 27 == 0 and not str_drew:
                    tmp_canvas.fill_text(str[:-2], pointer_x, pointer_y)
                    tmp_canvas.fill_text(str[:-2], pointer_x, pointer_y)
                    pointer_y += RT_Y_INTERVAL
                    str_drew = True
                    str = ""

        if not str_drew:
            tmp_canvas.fill_text(str[:-2], pointer_x, pointer_y)
            tmp_canvas.fill_text(str[:-2], pointer_x, pointer_y)

        self.related_tc_tree_canvas.height = pointer_y + RT_Y_INTERVAL
        self.related_tc_tree_canvas.draw_image(tmp_canvas, 0, 0)

        self.canvas_vbox_height += pointer_y + RT_Y_INTERVAL

        self.draw_canvas_frame(self.related_tc_tree_canvas)

    def draw_cluster_canvas(self):
        # Generate cluster from tcCanvas
        self.tc_canvas_tmp.differentiate_inexact_match = self.differentiate_inexact_match
        self.tc_canvas_tmp.cluster_from_ad_list()

        topology_list = []
        topology_set = set()
        ad_list = []
        for index, tc_tree_id in enumerate(self.tc_tree_list):
            tc_tree = self.adPy.tc[tc_tree_id-1]
            ad_tree = tc_tree.ad_tree
            topology_list.append(ad_tree.topology)

            topology_set.add(ad_tree.topology)

        for topology in topology_set:
            ad_list.append(topology.sample_ad_tree)

        self.tc_canvas_tmp.preprocess_ad_tree(ad_list)
        ad_height = DEFAULT_PADDING_BETWEEN_AD + (DEFAULT_AD_HEIGHT + CLUSTER_NUMBER_BAR_HEIGHT +
                                                  DEFAULT_PADDING_BETWEEN_BLOCK)
        self.cluster_list = []

        for topology in topology_set:
            ad_canvas = Canvas(width = DEFAULT_AD_WIDTH + 10 ,height = ad_height + 10)
            ad_canvas.tree_count = topology_list.count(topology)

            self.tc_canvas_tmp.draw_cluster_number_bar(topology.sample_ad_tree, topology_list.count(topology),
                                                       canvas=ad_canvas,agree_rt=topology.agree_rt)
            self.agree_rt = topology.agree_rt
            self.tc_canvas_tmp.draw_ad_tree(topology.sample_ad_tree, ad_canvas,cluster=True)
            ad_canvas.stroke_style = BLACK
            ad_canvas.stroke_rect(topology.sample_ad_tree.topL.x, topology.sample_ad_tree.topL.y,
                                  topology.sample_ad_tree.width, topology.sample_ad_tree.height)
            self.cluster_list.append(ad_canvas)


        # Draw cluster canvas

        canvas_height = math.ceil(len(self.cluster_list) / 4) * ad_height + 30
        self.cluster_canvas.height = canvas_height
        self.cluster_canvas.fill_style = BLACK
        self.cluster_canvas.font = f'18px {self.FONT_FAMILY}'

        pointer_x = 10
        pointer_y = 20
        self.cluster_canvas.fill_text("Cluster: ", pointer_x, pointer_y)
        pointer_y += 5

        for index,canvas in enumerate(self.cluster_list):
            if index > 0 and index % 4 == 0:
                pointer_y += ad_height
                pointer_x = 10

            self.cluster_canvas.draw_image(canvas,pointer_x,pointer_y)
            pointer_x += DEFAULT_AD_WIDTH + DEFAULT_PADDING_BETWEEN_AD

        self.canvas_vbox_height += canvas_height

        self.draw_canvas_frame(self.cluster_canvas)

    def draw_canvas_frame(self,canvas):
        canvas.stroke_style = GREY
        canvas.stroke_rect(0, 0, canvas.width, canvas.height)
        canvas.stroke_rect(2, 2, canvas.width - 4, canvas.height - 4)
    def draw_rt_subtree_block_canvas(self):
        # 先在reference tree画黑点，然后再crop出来
        rt_canvas = self.adPy.rt_canvas
        rt_canvas.clear_subtree_compare_canvas()
        nodes = self.segment_button_clicked.metadata.get('nodes_list')
        rt_canvas.draw_subtree_compare_nodes(nodes)

        self.rt_subtree_canvas_tmp = Canvas(width = rt_canvas.width,height = rt_canvas.height,sync_image_data=True)
        subtree_layer = rt_canvas.sorted_layer_list.index(self.subtree_chosen.label)

        self.rt_subtree_canvas_tmp.draw_image(rt_canvas[subtree_layer], 0, 0)
        self.rt_subtree_canvas_tmp.draw_image(rt_canvas[-1],0,0)
        self.rt_subtree_canvas_tmp.draw_image(rt_canvas[-3],0,0)

        block = self.subtree_chosen.block
        self.rt_subtree_block_canvas = Canvas(width=block.width, height=block.height)
        self.rt_subtree_block_vbox = widgets.VBox(children=[self.rt_subtree_block_canvas])

        self.rt_subtree_canvas_tmp.observe(self.crop_canvas, "image_data")

    def crop_canvas(self,*args, **kwargs):
        # subtree layer + rt_layer + subtree_compare_layer
        block = self.subtree_chosen.block
        self.rt_subtree_canvas_tmp.get_image_data()
        image_data = self.rt_subtree_canvas_tmp.get_image_data(block.topL.x, block.topL.y,
                                                                               block.width,block.height)

        self.rt_subtree_block_canvas.put_image_data(image_data,0,0)



    # ---------------------- Manage event -----------------------------
    def on_subtree_button_click(self,button):
        self.export_ready = False
        if self.subtree_button_clicked is not None:
            self.subtree_button_clicked.style.font_weight = 'normal'
            self.subtree_button_clicked.style.font_size = self.FONT_SIZE

            current_subtree_button_hbox = list(self.base_vbox.children)[0]

            self.base_vbox.children = [current_subtree_button_hbox]

        if button == self.subtree_button_clicked:
            self.subtree_button_clicked = None
            self.adPy.rt_canvas.clear_subtree_compare_canvas()
            return

        button.style.font_weight = 'bold'
        button.style.font_size = self.CHOSEN_FONT_SIZE
        self.subtree_button_clicked = button
        self.subtree_chosen = button.metadata.get('subtree')

        self.setup_subtree_segment()

    def on_segment_button_click(self,button):
        self.export_ready = False
        if self.segment_button_clicked is not None:
            self.segment_button_clicked.style.font_weight = 'normal'
            self.segment_button_clicked.style.font_size = self.FONT_SIZE

            self.canvas_vbox.children = []

            remain_children = list(self.base_vbox.children)[:2]
            self.base_vbox.children = remain_children

        if button == self.segment_button_clicked:
            self.segment_button_clicked = None
            self.adPy.rt_canvas.clear_subtree_compare_canvas()
            return

        button.style.font_weight = 'bold'
        button.style.font_size = self.CHOSEN_FONT_SIZE
        self.segment_button_clicked = button
        self.tc_tree_list = button.metadata.get('tree_list')



        self.setup_canvas_hbox()

    def on_treename_checkbox_change(self,change):
        if self.tree_name_checkbox.value:
            if not self.show_tree_name:
                self.show_tree_name = True
                self.related_tc_tree_canvas.clear()
                with hold_canvas(self.related_tc_tree_canvas):
                    self.draw_tc_tree_canvas()
        else:
            if self.show_tree_name:
                self.show_tree_name = False
                self.related_tc_tree_canvas.clear()
                with hold_canvas(self.related_tc_tree_canvas):
                    self.draw_tc_tree_canvas()

    def on_cluster_checkbox_change(self,change):
        if self.cluster_checkbox.value:
            self.differentiate_inexact_match = True
            self.cluster_canvas.clear()
            with hold_canvas(self.cluster_canvas):
                self.draw_cluster_canvas()

        else:
            self.differentiate_inexact_match = False
            self.cluster_canvas.clear()
            with hold_canvas(self.cluster_canvas):
                self.draw_cluster_canvas()


    # ---------------------- Setup Tree Distribution View ---------------
    def setup_subtree_button(self):
        self.generate_subtree_button()
        self.generate_subtree_button_hgridbox(self.subtree_button_list)
        self.draw_base_vbox([self.subtree_button_hbox])

    def setup_subtree_segment(self):
        self.generate_segment_button_hbox()
        # self.draw_base_vbox([self.segment_button_hbox])

    def setup_canvas_hbox(self):
        self.nodes_list_canvas = Canvas(width=self.DEFAULT_CANVAS_WIDTH, height=self.DEFAULT_CANVAS_WIDTH)
        with hold_canvas(self.nodes_list_canvas):
            self.draw_nodes_canvas()


        height = RT_Y_INTERVAL * len(self.tc_tree_list)
        self.related_tc_tree_canvas = Canvas(width=self.DEFAULT_CANVAS_WIDTH, height=height)
        with hold_canvas(self.related_tc_tree_canvas):
            self.draw_tc_tree_canvas()

        self.cluster_canvas = Canvas(width=self.DEFAULT_CANVAS_WIDTH, height=self.DEFAULT_CANVAS_WIDTH)
        with hold_canvas(self.cluster_canvas):
            self.draw_cluster_canvas()

        self.generate_canvas_vgridbox([self.related_tc_tree_canvas,self.cluster_canvas,self.nodes_list_canvas])

        self.draw_rt_subtree_block_canvas()

        self.draw_canvas_hbox([self.canvas_vbox, self.hspacer,self.rt_subtree_block_vbox])
        self.export_ready = True



