from ipycanvas import Canvas, hold_canvas
from myCanvas import *
from myUtils import *
import math

import copy

# Canvas for Tree Collection
class tcCanvas(MyCanvas):
    # Bookmark - tc
    id = None
    ad_list = []
    rt_ad = None

    # layer -1 for subtree comparison view
    # layer -2 for tree name
    # layer -3 forward for ad
    def __init__(self,adPy,view,ad_per_row=DEFAULT_AD_PER_ROW,width=1100, height=1100,scale=1.0,max_ad=None,
                 context_level=2,filter=INCLUDE,first_ad=None,last_ad=None,tree_id=None,tree_name=None,sort_by=ID,
                 escape_taxa_as_context_block=True,show_block_proportional=False,parameter_from_individual_ad=True,
                 subtree_independent=True,differentiate_inexact_match=True,layer=5,show_tree_name=False):
        super().__init__(layer = layer + 4, width = width, height = height)

        # Common attribute
        self.AD_LAYER = -3
        self.TREE_NAME_LAYER = -2
        self.COMPARISON_VIEW_LAYER = -1
        self.TREE_SELECTED_BLOCK_LAYER = 0

        self.adPy = adPy
        self.ad_list = []
        self.scale = scale
        self.max_ad = max_ad
        self.set_tc_canvas_default()
        self.ad_per_row = ad_per_row
        self.context_level = context_level
        self.scale_alter = False
        self.check_result_list = []
        self.filter_type = filter

        self.cluster_agree_rt = None

        self.view = view
        self.adPy = adPy
        self.rt = adPy.rt
        self.tc = self.adPy.sort_tc(tc_list=self.adPy.tc, sort_by=sort_by)

        self.escape_taxa_as_context_block = escape_taxa_as_context_block  # True if show all independent leaf
        self.show_block_proportional = show_block_proportional  # True if colored segment width proportional with block width
        self.parameter_from_individual_ad = parameter_from_individual_ad  # Set parameter of cluster ad same as individual ad canvas
        self.subtree_independent = subtree_independent


        self.duplicate_subtree = None

        # Individual_ad
        self.first_ad = first_ad
        self.last_ad = last_ad
        self.target_tc_tree_index = tree_id
        self.target_tc_tree_name = tree_name
        self.ad_drawn = 0
        self.sort_by = sort_by
        self.show_tree_name = show_tree_name
        self.tree_hover = None
        self.tree_selected = None

        self.topology_list = []
        self.topology_canvas_list = []
        self.prepare_subtree_list()

        # Construct rt ad and get topology
        check_elided = {}
        check_elided['node'] = None
        check_elided['context_level'] = 0
        self.construct_ad(tc_tree=self.adPy.rt, level=1, check_elided=check_elided)
        self.adPy.rt.ad_tree = self.ad_list[0]
        self.adPy.rt.topology_string = self.adPy.rt.ad_tree.ad_to_string(self,differentiate_inexact_match=False)

        self.ad_list = []
        # Cluster_ad
        self.differentiate_inexact_match = differentiate_inexact_match

        if view == AD_INDIVIDUAL:
            self.individual_ad()
            self.on_mouse_down(self.mouse_clicked)
            self.on_mouse_move(self.mouse_hover)

        elif view == AD_CLUSTER:
            self.cluster_ad()

        elif view == TREE_DISTRIBUTION:
            # self.prepare_subtree_list()
            self.generate_ad_list()

    #------------------- Individual ad view ----------------------#
    def individual_ad(self):

        if not self.target_tc_tree_index and not self.target_tc_tree_name:
            self.condition_exist = False
        else:
            self.condition_exist = True

        self.generate_ad_list()

        self.construct_individual_ad_canvas()

    def prepare_subtree_list(self):
        for subtree in self.adPy.subtree_list:
            subtree.get_leaf_nodes()
            subtree.topology_list = {}

        self.adPy.subtree_list = sorted(self.adPy.subtree_list, key=lambda x: len(x.leaf_set), reverse=True)

    def generate_ad_list(self):
        for index, tc_tree in enumerate(self.tc):
            if self.view == AD_INDIVIDUAL:
                if self.condition_exist and self.check_index_name_fulfill():
                    break

                if self.filter_type == INCLUDE:
                    result = self.check_include_condition(index, tc_tree)
                else:
                    result = self.check_exclude_condition(tc_tree)

                if result == CONTINUE:
                    continue
                elif result == BREAK:
                    break

            check_elided = {}
            check_elided['node'] = None
            check_elided['context_level'] = 0
            self.construct_ad(tc_tree=tc_tree, level=1, check_elided=check_elided)

            self.ad_drawn += 1

    def check_index_name_fulfill(self):
        if self.target_tc_tree_index or self.target_tc_tree_name:
            return False

        # Condition fulfill, break process
        return True

    def check_include_condition(self, index, tc_tree):
        if self.first_ad and index < self.first_ad - 1:
            return CONTINUE
        if self.last_ad and index >= self.last_ad:
            return BREAK
        if self.max_ad and self.ad_drawn >= self.max_ad:
            return BREAK

        check_tree_index_name = TRUE
        if self.target_tc_tree_index:
            if tc_tree.id in self.target_tc_tree_index:
                self.target_tc_tree_index.remove(tc_tree.id)
                if len(self.target_tc_tree_index) == 0:
                    self.target_tc_tree_index = None
                return TRUE
            else:
                check_tree_index_name = CONTINUE

        if self.target_tc_tree_name:
            if tc_tree.name in self.target_tc_tree_name:
                self.target_tc_tree_name.remove(tc_tree.name)
                if len(self.target_tc_tree_name) == 0:
                    self.target_tc_tree_name = None
                return TRUE
            else:
                check_tree_index_name = CONTINUE

        return check_tree_index_name

    def check_exclude_condition(self, tc_tree):
        if self.first_ad and self.last_ad and tc_tree.id >= self.first_ad and tc_tree.id <= self.last_ad:
            return CONTINUE

        if self.max_ad and self.ad_drawn >= self.max_ad:
            return BREAK

        check_tree_index_name = TRUE
        if self.target_tc_tree_index:
            if tc_tree.id in self.target_tc_tree_index:
                check_tree_index_name = CONTINUE
            else:
                return TRUE

        if self.target_tc_tree_name:
            if tc_tree.name in self.target_tc_tree_name:
                check_tree_index_name = CONTINUE
            else:
                return TRUE

        return check_tree_index_name

    def reset_tc_canvas(self):
        # for layer in range(0,self.layer):
        self[self.AD_LAYER].clear()
        self.set_tc_canvas_default()

    def set_tc_canvas_default(self):
        self.padding = DEFAULT_PADDING_BETWEEN_AD  # Padding between ADs
        self.ad_width = DEFAULT_AD_WIDTH * self.scale
        self.ad_height = DEFAULT_AD_HEIGHT * self.scale


    # ------------- Drawing Related ---------------- #
    def check_canvas_height(self,ad_list):
        if len(ad_list) > 0:
            self.height = ad_list[-1].botR.y + 20
        else:
            self.height = 10

    def construct_individual_ad_canvas(self):

        result = self.preprocess_ad_tree(self.ad_list)

        if result == SPACE_INSUFFICIENT:
            self.adjust_display()
            return

        self.check_canvas_height(self.ad_list)

        if self.scale_alter:
            self.display_content_exceed_error()
            # self.height += (self.height * self.scale) * 2

        self.cluster_from_ad_list()
        self.generate_topology_canvas()

        # self.draw_ad_tree_rec(self.ad_list)

    def construct_ad(self,tc_tree,construct_rt=False,tc_node=None,current_ad_node=None,level=0,nested_subtree=None,nested_ad=None,check_elided={},nested_subtree_duplicate=False):

        if tc_node is None and current_ad_node is None: # When start from root
            tc_node = tc_tree.seed_node
            tc_node.index = 0

            ad = AD_Tree(id=tc_tree.id,tc_tree=tc_tree,tc_canvas=self)  # Create a new AD_Tree which belong to
            tc_tree.ad_tree = ad
            # current tc_tree
            ad.root = AD_Node(node_or_block=tc_node,x_level=level,type=ROOT,child_index=None,type_prior=None)

            if len(tc_tree.missing) > 0:
                ad.missing_taxa_list = tc_tree.missing

            # Set ad_tree's root
            current_ad_node = ad.root
            self.ad_constructing = ad
            if construct_rt:
                self.rt_ad = ad
            else:
                self.ad_list.append(ad)
            check_elided['node'] = current_ad_node
            check_elided['context_level'] = 0

            result = None
            for subtree in self.adPy.subtree_list:
                result = self.check_relation(subtree=subtree,node=tc_node,tree_index = tc_tree.id,nested_subtree_duplicate=nested_subtree_duplicate)
                if result == SUBTREE_ROOT:
                    self.generate_block_by_result(result=result, subtree=subtree, tc_tree=tc_tree, tc_node=tc_node,
                                                  current_ad_node=current_ad_node , level=level,nested_subtree=nested_subtree, nested_ad=nested_ad,check_elided=check_elided,nested_subtree_duplicate=nested_subtree_duplicate)

                    return

            # if result == NON_DESCENDANT:
            #     self.generate_block_by_result(result=result, subtree=None, tc_tree=tc_tree, tc_node=tc_node,
            #                                   current_ad_node=current_ad_node, level=level,
            #                                   nested_subtree=nested_subtree, nested_ad=nested_ad,
            #                                   check_elided=check_elided,
            #                                   nested_subtree_duplicate=nested_subtree_duplicate)
            #     return


        result_list = {}

        for index,child in enumerate(tc_node.child_node_iter()):
            # print("============")
            # print("Child:" ,end="")
            # print(child)
            child.index = index
            result = NON_DESCENDANT
            for subtree in self.adPy.subtree_list:
                if nested_subtree and subtree.label in nested_subtree:
                    continue

                result_tmp = self.check_relation(subtree=subtree,node=child,tree_index = tc_tree.id,nested_subtree_duplicate=nested_subtree_duplicate)

                if result_tmp < result:
                    result = result_tmp

                if result == SUBTREE_ROOT or result == INDEPENDENT_LEAF:
                    # subtree_prior = 5 - SUBTREE_LABEL_LIST.index(subtree.label)
                    # print(f"subtree_prior = {subtree_prior}")
                    # result -= subtree_prior
                    # print(f"result = {result}")
                    break

            result_list[child] = [result,subtree]
            if tc_tree.id > 0:
                self.check_result_list.append(result_list[child])

        # Sort by result
        sorted_list = dict(sorted(result_list.items(), key=lambda item: item[1][0]))
        check_result = { value[0] for value in sorted_list.values() }


        if len(check_result) == 1:
            # Elided stop here
            if check_elided['node'] and check_elided['context_level'] >= self.context_level:
                ori_child = check_elided['node'].children.pop()

                parent = check_elided['node']
                current_ad_node.child_index = ori_child.child_index
                current_ad_node.x_level = ori_child.x_level + 1
                current_ad_node.is_elided = True
                check_elided['node'] = None
                check_elided['context_level'] = 0

                if nested_ad:
                    nested_ad.insert_node(parent, current_ad_node)
                else:
                    self.ad_constructing.insert_node(parent, current_ad_node)

        for node,result in sorted_list.items():
            self.generate_block_by_result(result=result[0],subtree=result[1],tc_tree=tc_tree,tc_node=node,current_ad_node=current_ad_node,level=level,nested_subtree=nested_subtree,nested_ad=nested_ad,check_elided=check_elided,nested_subtree_duplicate=nested_subtree_duplicate)


    def generate_block_by_result(self,result,subtree,tc_tree,tc_node=None,current_ad_node=None,level=0,nested_subtree=None,nested_ad=None,check_elided={},nested_subtree_duplicate=False):
        if result == INDEPENDENT_LEAF:
            # print("INDEPENDENT_LEAF")
            # if self.escape_taxa_as_context_block and not nested_subtree_duplicate:
            #     return

            # print("INDEPENDENT_LEAF")
            # print("Subtree:" + subtree.label)
            # width,height,subtree,ad_tree,type,topL=None, botR=None
            independent_leaf_block = self.individual_block(subtree=subtree, type=INDIVIDUAL_LEAF_BLOCK)
            independent_leaf_block.color = subtree.color
            independent_leaf_block.root = tc_node
            independent_leaf_block.taxa_list.add(tc_node.taxon.label)
            independent_leaf_block.subtree_taxa_count = 1
            new_ad_node = AD_Node(node_or_block=independent_leaf_block, x_level=level, type=LEAF,
                                        child_index=tc_node.index,type_prior=result)


            if not self.escape_taxa_as_context_block and check_elided['node'] and check_elided['context_level'] >=self.context_level:
                ori_child = check_elided['node'].children.pop()
                current_ad_node.child_index = ori_child.child_index
                if nested_ad:
                    nested_ad.insert_node(check_elided['node'], current_ad_node)
                else:
                    self.ad_constructing.insert_node(check_elided['node'], current_ad_node)

                current_ad_node.x_level = ori_child.x_level + 1
                new_ad_node.x_level = current_ad_node.x_level + 1
                current_ad_node.is_elided = True
                check_elided['node'] = None
                check_elided['context_level'] = 0

            elif not self.escape_taxa_as_context_block:
                check_elided['node'] = None
                check_elided['context_level'] = 0

            if nested_ad:
                nested_ad.insert_node(current_ad_node, new_ad_node)
            else:
                self.ad_constructing.insert_node(current_ad_node, new_ad_node)


            # current_ad_node.children.append(new_ad_node)

        elif result == SUBTREE_ROOT or result == DUPLICATE_SUBTREE_ROOT:
            # print("SUBTREE_ROOT")
            # taxa_list = all taxa under corresponding subtree
            # subtree_taxa_count = taxa of reference subtree exist in this corresponding subtree
            # if corresponding subtree exactly same as reference subtree

            if result == DUPLICATE_SUBTREE_ROOT:
                result = self.check_subtree_combine(tc_node, tree_index=tc_tree.id)
                subtree_block = self.ad_subtree_block(subtree=self.duplicate_subtree,subtree_duplicate=True)
                subtree_block.color = []
                subtree_block.taxa_list = {leaf.taxon.label for leaf in tc_node.leaf_nodes()}
                # subtree_block.duplicate_subtree_taxa_count
                subtree_block.exact_match = []

                for index,duplicate_subtree in enumerate(self.duplicate_subtree):
                    # Calculate subtree topology
                    duplicate_subtree.check_and_set_topology(tc_node,tc_tree.id)

                    subtree_block.color.append(duplicate_subtree.color)
                    subtree_block.duplicate_subtree_taxa_count.append(len(subtree_block.taxa_list.intersection(duplicate_subtree.leaf_set)))

                    if subtree_block.taxa_list == duplicate_subtree.leaf_set:
                        subtree_block.exact_match.append(True)
                    else:
                        subtree_block.exact_match.append(False)

                for count in subtree_block.duplicate_subtree_taxa_count:
                    subtree_block.subtree_taxa_count += count



            else:
                subtree_block = self.ad_subtree_block(subtree=subtree,subtree_duplicate=False)
                subtree.check_and_set_topology(tc_node,tc_tree.id)

                subtree_block.color = subtree.color
                subtree_block.taxa_list = {leaf.taxon.label for leaf in tc_node.leaf_nodes()}
                subtree_block.subtree_taxa_count = len(subtree_block.taxa_list.intersection(subtree.leaf_set))

                if subtree_block.taxa_list == subtree.leaf_set:
                    subtree_block.exact_match = True
                else:
                    subtree_block.exact_match = False

            subtree_block.taxa_list = {leaf.taxon.label for leaf in tc_node.leaf_nodes()}
            subtree_block.root = tc_node

            new_ad_node = AD_Node(node_or_block=subtree_block, x_level=level, type=LEAF,
                                        child_index=tc_node.index,type_prior=result)

            if check_elided['node'] and check_elided['context_level'] >= self.context_level:
                ori_child = check_elided['node'].children.pop()

                parent = check_elided['node']
                if not subtree_block.exact_match and not self.escape_taxa_as_context_block:
                    current_ad_node.child_index = ori_child.child_index
                    if nested_ad:
                        nested_ad.insert_node(parent, current_ad_node)
                    else:
                        self.ad_constructing.insert_node(parent, current_ad_node)
                    current_ad_node.x_level = ori_child.x_level + 1
                    new_ad_node.x_level = current_ad_node.x_level + 1
                    current_ad_node.is_elided = True
                    check_elided['node'] = None
                    check_elided['context_level'] = 0

                    if nested_ad:
                        nested_ad.insert_node(current_ad_node, new_ad_node)
                    else:
                        self.ad_constructing.insert_node(current_ad_node, new_ad_node)
                else:
                    new_ad_node.child_index = ori_child.child_index
                    new_ad_node.x_level = ori_child.x_level + 1
                    new_ad_node.is_elided = True
                    subtree_block.is_elided = True
                    if nested_ad:
                        nested_ad.insert_node(parent, new_ad_node)
                    else:
                        self.ad_constructing.insert_node(parent, new_ad_node)

            else:
                check_elided['node'] = None
                check_elided['context_level'] = 0
                if nested_ad:
                    nested_ad.insert_node(current_ad_node, new_ad_node)
                else:
                    self.ad_constructing.insert_node(current_ad_node, new_ad_node)


            # current_ad_node.children.append(new_ad_node)

            if tc_node.is_leaf():
                return

            # self.ad_constructing.plot_tree()

            # Check nested subtree/block
            for check_subtree in self.adPy.subtree_list:
                # result = self.check_subtree_combine(tc_node, tree_index=tc_tree.id)
                if result == DUPLICATE_SUBTREE_ROOT:
                    if check_subtree in self.duplicate_subtree:
                        continue

                if check_subtree == subtree or (nested_subtree and check_subtree.label in nested_subtree):
                    continue

                if self.check_leaf_node_descendant(check_subtree, tc_node) == IS_DESCENDANT:
                    # print("has nested")
                    # print("outer subtree = " + subtree.label)
                    # print("inner subtree = " +  check_subtree.label)

                    if nested_subtree:
                        nested_subtree += subtree.label
                    else:
                        nested_subtree = subtree.label

                    self.ad_constructing.nested_cnt += 1
                    nested_ad = AD_Tree(id=self.ad_constructing.id, tc_tree=tc_tree, tc_canvas=self)
                    nested_ad.root = AD_Node(node_or_block=tc_node, x_level=0, type=ROOT,child_index=None,
                                                   type_prior=result)

                    new_check_elided = {}
                    new_check_elided['node'] = nested_ad.root
                    new_check_elided['context_level'] = 0

                    nested_ad.is_nested = True
                    new_ad_node.nested_tree = nested_ad
                    subtree_block.nested_tree = nested_ad

                    if self.check_corresponding_branch_duplicate(check_subtree,nested_subtree,tc_tree.id):
                        nested_subtree_duplicate=True

                    self.construct_ad(tc_tree=tc_tree, tc_node=tc_node, current_ad_node=nested_ad.root, level=1,
                                      nested_subtree=nested_subtree, nested_ad=nested_ad,check_elided=new_check_elided,
                                      nested_subtree_duplicate=nested_subtree_duplicate)

                    break

            # self.duplicate_subtree = None

        elif result == IS_DESCENDANT:
            # print("IS_DESCENDANT")
            new_ad_node = AD_Node(node_or_block=tc_node, x_level=level, type=INTERNAL,
                                        child_index=tc_node.index,type_prior=result)

            if check_elided['node'] != None:
                check_elided['context_level'] += 1
            else:
                check_elided['node'] = new_ad_node


            if nested_ad:
                nested_ad.insert_node(current_ad_node, new_ad_node)
            else:
                self.ad_constructing.insert_node(current_ad_node, new_ad_node)
            # current_ad_node.children.append(new_ad_node)
            self.construct_ad(tc_tree=tc_tree, tc_node=tc_node, current_ad_node=new_ad_node, level=current_ad_node.x_level + 1,
                              nested_subtree=nested_subtree, nested_ad=nested_ad,check_elided=check_elided,nested_subtree_duplicate=nested_subtree_duplicate)

        elif result == NON_DESCENDANT:
            # print("NON_DESCENDANT")
            non_descendant_block = self.individual_block(subtree=None, type=INDIVIDUAL_BLOCK)
            non_descendant_block.color = BLANK
            non_descendant_block.root = tc_node
            new_ad_node = AD_Node(node_or_block=non_descendant_block, x_level=level, type=LEAF,
                                        child_index=tc_node.index,type_prior=result)

            if nested_ad:
                nested_ad.insert_node(current_ad_node, new_ad_node)
            else:
                self.ad_constructing.insert_node(current_ad_node, new_ad_node,)
            # current_ad_node.children.append(new_ad_node)

    def check_corresponding_branch_duplicate(self,first_subtree,subtree_label,tree_index):
        for char in subtree_label:
            for subtree in self.adPy.subtree_list:
                if subtree.label == char and subtree != first_subtree:
                    if first_subtree.root.corr[tree_index] == subtree.root.corr[tree_index]:
                        return True

    # Function to check whether current branch/node related with target subtree/taxa
    def check_relation(self,subtree,node,tree_index,nested_subtree_duplicate):
        # If node is subtree's root
        if node == subtree.root.corr[tree_index]:
            result = self.check_subtree_combine(node,tree_index=tree_index)
            if not self.subtree_independent and result == DUPLICATE_SUBTREE_ROOT:
                return DUPLICATE_SUBTREE_ROOT
            else:
                return SUBTREE_ROOT

        # If node from tc is target leaf node in this subtree
        if node.is_leaf():
            if node.taxon.label in subtree.leaf_set:
                return INDEPENDENT_LEAF

        if self.check_leaf_node_descendant(subtree,node) == IS_DESCENDANT:
            if self.escape_taxa_as_context_block:
                if nested_subtree_duplicate or self.check_subtree_in_descendant(node,tree_index=tree_index):
                    return IS_DESCENDANT
            else:
                return IS_DESCENDANT

        # If current branch has no target taxa
        return NON_DESCENDANT

    def check_leaf_node_descendant(self,subtree,node):
        leaf_nodes = {leaf.taxon.label for leaf in node.leaf_nodes()}

        if leaf_nodes.intersection(subtree.leaf_set):
            return IS_DESCENDANT

        else:
            return False

    # Check if node descendant has subtree
    def check_subtree_in_descendant(self,node,tree_index):
        for subtree in self.adPy.subtree_list:
            if self.adPy.is_descendant(parent_node=node,child_node=subtree.root.corr[tree_index]):
                return True

        return False

    def check_subtree_combine(self,node,tree_index):
        check_list = []

        for subtree in self.adPy.subtree_list:
            if node == subtree.root.corr[tree_index]:
                check_list.append(subtree)

        if len(check_list) > 1:
            self.duplicate_subtree = check_list

            return DUPLICATE_SUBTREE_ROOT
        else:
            return False

    # Function for ad_drawing preprocess
    def preprocess_ad_tree(self,ad_list):
        self.set_tc_canvas_default()

        for layer in range(0,self.layer):
            self[layer].ad_list = []

        # Calculate AD size, AD per row and padding
        self.ad_width = DEFAULT_AD_WIDTH * self.scale
        self.ad_height = DEFAULT_AD_HEIGHT * self.scale

        while (self.width - ( self.ad_width * self.ad_per_row + 2 * self.padding )) < (self.ad_per_row-2) * self.padding :
            self.ad_per_row -= 1
        # Testing
        # print("ad_width:" + str(self.ad_width))
        # print("ad_height:" + str(self.ad_height))
        # print("ad_per_row:" + str(self.ad_per_row))
        # print("padding:" + str(self.padding))

        ad_y = DEFAULT_PADDING_BETWEEN_AD
        if self.view == AD_CLUSTER or TREE_DISTRIBUTION:
            ad_y += CLUSTER_NUMBER_BAR_HEIGHT * self.scale + DEFAULT_PADDING_BETWEEN_BLOCK
        elif self.view == AD_INDIVIDUAL and self.show_tree_name:
                ad_y += TREE_NAME_LABEL_HEIGHT + DEFAULT_PADDING_BETWEEN_BLOCK

        current_row = 0
        for index,ad_tree in enumerate(ad_list):
            #
            # if index in [2,4,5,6,7,8,9,10,11,12,16,17,18,19]:
            #     ad_tree.missing_taxa_list = {'Monomastix opisthostigma'}
            if ad_tree.missing_taxa_list:
                for subtree in self.adPy.subtree_list:
                    missing_taxa_from_subtree = subtree.leaf_set.intersection(ad_tree.missing_taxa_list)
                    # print(f"missing_taxa_from_subtree {subtree.label}:")
                    # print(missing_taxa_from_subtree)
                    if len(missing_taxa_from_subtree) > 0:
                        missing_segment = Missing_Taxa_Segment(subtree,ad_tree,missing_taxa_from_subtree)
                        ad_tree.missing_taxa_segment.append(missing_segment)

            ad_tree.generate_block_list()
            # Set ad_tree position and size
            # Check if line break
            if index % self.ad_per_row == 0 and index != 0:
                current_row += 1
                ad_y +=DEFAULT_AD_HEIGHT * self.scale + self.padding
                if self.view == AD_CLUSTER:
                    ad_y += CLUSTER_NUMBER_BAR_HEIGHT * self.scale + DEFAULT_PADDING_BETWEEN_BLOCK
                elif self.view == AD_INDIVIDUAL and self.show_tree_name:
                    ad_y += TREE_NAME_LABEL_HEIGHT * self.scale + DEFAULT_PADDING_BETWEEN_BLOCK

            # Calculate AD_Tree position in row
            ad_index = index % self.ad_per_row
            ad_tree.index_in_row = ad_index
            ad_x = self.padding + (ad_index * self.ad_width) + (ad_index * (self.padding-5))

            ad_tree.set_position_size(ad_x,ad_y)
            ad_tree.located_layer = self.AD_LAYER

            # if self.default_value:
            ad_tree.set_default()

            # if index == 0:
            self.indv_block_width = DEFAULT_INDV_BLOCK_WIDTH * self.scale
            self.indv_block_height = DEFAULT_INDV_BLOCK_HEIGHT * self.scale
            self.nested_block_height = NESTED_BLOCK_MINIMUM_HEIGHT * self.scale
            self.block_minimum_height = BLOCK_MINIMUM_HEIGHT[len(self.adPy.subtree_list) % MAX_SUBTREE] * self.scale
            self.missing_taxa_height = DEFAULT_MISSING_TAXA_BLOCK_HEIGHT * self.scale

            # Testing
            # print("indv_block_width = " + str(self.indv_block_width))
            # print("indv_block_height = " + str(self.indv_block_height))
            # print("nested_block_height = " + str(self.nested_block_height))
            # print("block_minimum_height = " + str(self.block_minimum_height))

            sufficient_space = self.adjust_ad_block(ad_tree)
            # print("Sufficient space = " + str(sufficient_space))
            while sufficient_space < 0:
                # print("Space not enoughh")
                # reduction = abs(sufficient_space)
                # self.default_value = False

                indv_block_height_limit = False if self.indv_block_height > 4 else True
                nested_block_height_limit = False if self.nested_block_height > 4 else True
                block_height_limit = False if self.block_minimum_height > 7 else True
                padding_limit = False if ad_tree.padding > 2 else True
                missing_taxa_block_limit = False if ad_tree.padding > 7 else True

                if indv_block_height_limit and nested_block_height_limit and block_height_limit and padding_limit and missing_taxa_block_limit:
                    return SPACE_INSUFFICIENT

                if not indv_block_height_limit:
                    self.indv_block_height -= 1
                if not nested_block_height_limit:
                    self.nested_block_height -= 2
                if not block_height_limit:
                    self.block_minimum_height -= 5
                if not padding_limit:
                    ad_tree.padding -= 1
                if not missing_taxa_block_limit:
                    self.missing_taxa_height -= 2

                # Testing
                # print("indv_block_width = " + str(self.indv_block_width))
                # print("indv_block_height = " + str(self.indv_block_height))
                # print("nested_block_height = " + str(self.nested_block_height))
                # print("block_minimum_height = " + str(self.block_minimum_height))

                sufficient_space = self.adjust_ad_block(ad_tree)
                # print("Sufficient space = " + str(sufficient_space))
        return SPACE_SUFFICIENT

    def adjust_display(self):
        self.scale_alter = True
        self.scale += 0.1
        self.reset_tc_canvas()
        if self.view == AD_INDIVIDUAL:
            self.construct_individual_ad_canvas()
        else:
            self.construct_cluster_ad_canvas()

    def adjust_ad_block(self,ad_tree):
        # Calculate block size in each ad_tree
        space_required = (15 + ad_tree.padding) if ad_tree.is_nested else (ad_tree.y + 2 * ad_tree.padding *
                                                                           self.scale)

        if ad_tree.missing_taxa_list:
            space_required += self.missing_taxa_height + ad_tree.padding * self.scale

        indv_block_list, subtree_block_list,unnested_block_taxa = ad_tree.individual_subtree_block_list()
        for indv_block in indv_block_list:
            indv_block.width = self.indv_block_width
            indv_block.height = self.indv_block_height
            space_required += indv_block.height + ad_tree.padding * self.scale

        free_space = 0
        for index,subtree_block in enumerate(subtree_block_list):
            if subtree_block.nested_tree:
                subtree_block.nested_tree.height = self.adjust_ad_block(subtree_block.nested_tree)
                subtree_block.height = subtree_block.nested_tree.height
                # print(f"{subtree_block.belong_subtree.label} height = "+str(subtree_block.height))
                space_required += subtree_block.height + ad_tree.padding * self.scale

                if index == len(subtree_block_list)-1:
                    if ad_tree.is_nested:
                        space_required += ad_tree.padding * self.scale
                    else:
                        free_space = ad_tree.height - space_required
                        subtree_block.height += free_space
                        subtree_block.nested_tree.height = subtree_block.height

                        self.allocate_remaining_space(subtree_block.nested_tree)

            elif ad_tree.is_nested:
                subtree_block.height = self.nested_block_height
                space_required += subtree_block.height + ad_tree.padding * self.scale


            else:
                # Outer subtree block
                block_remaining = len(subtree_block_list) - index
                # print("block_remaining = " + str(block_remaining))
                minimum_height = block_remaining * self.block_minimum_height

                # print("space_required = " + str(space_required))
                # print("minimum_height = " + str(minimum_height))
                # print("This block is outer block,height = ", end="")

                if (self.ad_height - space_required - block_remaining * ad_tree.padding * self.scale) <= minimum_height:
                    # print("Return -minimum_height")
                    return -minimum_height
                else:
                    if free_space == 0:
                        free_space = self.ad_height - space_required - minimum_height - block_remaining * ad_tree.padding * self.scale

                    subtree_block.height = self.block_minimum_height + (
                            (subtree_block.subtree_taxa_count / unnested_block_taxa) * free_space)

                    # Testing
                    # print(subtree_block.height)

                    space_required += subtree_block.height + ad_tree.padding * self.scale

        ad_tree.space_required = space_required

        if ad_tree.is_nested:
            return space_required

        return self.ad_height - space_required

    def allocate_remaining_space(self,ad_tree):

        indv_block_list, subtree_block_list,taxa_total = ad_tree.individual_subtree_block_list()
        block = subtree_block_list[-1]
        if ad_tree.is_nested and block.nested_tree:
            free_space = ad_tree.height - ad_tree.space_required
            block.height += free_space

            if block.nested_tree:
                self.allocate_remaining_space(block.nested_tree)

    # Before draw_ad, must confirm these variable has been altered according scale
    # canvas.padding(padding might change when ad_per_row change)
    # AD_HEIGHT / AD_WIDTH / ad_per_row /
    # Some variable remain the same
    # ad_tree.padding / ad_tree.initial_x,initial_y /
    # Function to draw ad_tree in each AD
    # Statement in this function will only involve one ad_tree
    def paste_ad_tree_canvas(self,ad_tree):
        self.draw_ad_tree_rec(self.ad_list)

        self[ad_tree.located_layer].draw_image(ad_tree.topology.ad_canvas,ad_tree.topL.x,ad_tree.topL.y)

        if not ad_tree.is_nested and self.show_tree_name:
            self.write_tree_name(ad_tree)
            self.generate_tree_name_block(ad_tree)

            # Testing
            # self.draw_frame(x,y,self.ad_width,TREE_NAME_LABEL_HEIGHT,BLACK,ad_tree.located_layer)

        # Generate ad_tree.missing_taxa_block and draw the block
        if ad_tree.missing_taxa_list:
            x = ad_tree.topL.x + ad_tree.x * self.scale
            y = ad_tree.botR.y - 2 * ad_tree.padding * self.scale - self.missing_taxa_height
            width = ad_tree.botR.x - ad_tree.x * self.scale - x
            ad_tree.missing_taxa_block = Block(Point(x,y),Point(x+width,y+self.missing_taxa_height))

            self.draw_missing_taxa_block(ad_tree)

    def draw_ad_tree(self,ad_tree=None,canvas=None,cluster=False):
        with hold_canvas(self):
            ad_tree.y *= self.scale
            for node in ad_tree.traverse_postorder():
                if canvas and not ad_tree.is_nested:
                    if cluster:
                        ad_topL_x = 5
                        ad_topL_y = CLUSTER_NUMBER_BAR_HEIGHT + DEFAULT_PADDING_BETWEEN_BLOCK * 2 + 5

                        ad_tree.topL.x = ad_topL_x
                        ad_tree.topL.y = ad_topL_y
                    else:
                        ad_topL_x = 0
                        ad_topL_y = 0

                    ad_botR_x = self.ad_width
                    ad_botR_y = self.ad_height

                else:
                    ad_topL_x = ad_tree.topL.x
                    ad_topL_y = ad_tree.topL.y
                    ad_botR_x = ad_tree.botR.x
                    ad_botR_y = ad_tree.botR.y

                if node.type == LEAF:
                    block = node.node_or_block

                    scale_tmp = self.scale

                    if self.scale > 1.0:
                        scale_tmp = 1.0
                    elif self.scale < 0.6:
                        scale_tmp = 0.6

                    # scale_tmp = 0.7 if self.scale < 0.7 else self.scale

                    x = ad_topL_x + (ad_tree.x * scale_tmp) + (
                            ad_tree.x * node.x_level * scale_tmp)  # ad_pos + padding

                    # Calculate block position
                    if ad_tree.is_nested:
                        y = ad_topL_y + ad_tree.y + (12 * self.scale)
                    else:
                        y = ad_topL_y + ad_tree.y

                    block.topL = Point(x, y)

                    block.botR = Point(ad_botR_x - ad_tree.x * self.scale, y + block.height)

                    # Calculate branch position and generate node's branch
                    if block.is_elided:
                        node.construct_elided_branch(x, y + (block.height / 2), scale_tmp)
                    else:
                        node.construct_branch(x,y + (block.height / 2),scale_tmp)

                        # Draw block
                    if block.type == INDIVIDUAL_BLOCK or block.type == INDIVIDUAL_LEAF_BLOCK:
                        self.draw_indv_block(block=block,ad_tree=ad_tree,canvas=canvas)

                    elif block.type == SUBTREE_BLOCK:
                        self.draw_subtree_block(block=block,ad_tree=ad_tree,scale_tmp=scale_tmp,canvas=canvas,cluster=False)

                    if ad_tree.is_nested:
                        ad_tree.y += block.height + ad_tree.padding
                    else:
                        ad_tree.y += block.height + ad_tree.padding * scale_tmp

                if node.type == INTERNAL or node.type == ROOT:
                    self.draw_root_internal(node,layer_index=ad_tree.located_layer,canvas=canvas)

            if canvas is not None:
                return

        # if not ad_tree.is_nested and self.show_tree_name:
        #     self.write_tree_name(ad_tree)
        #     self.generate_tree_name_block(ad_tree)
        #
        #     # Testing
        #     # self.draw_frame(x,y,self.ad_width,TREE_NAME_LABEL_HEIGHT,BLACK,ad_tree.located_layer)
        #

        # Generate ad_tree.missing_taxa_block and draw the block
        if (self.view == AD_CLUSTER or self.view == TREE_DISTRIBUTION) and ad_tree.missing_taxa_list:
            x = ad_tree.topL.x + ad_tree.x * self.scale
            y = ad_tree.botR.y - 2 * ad_tree.padding * self.scale - self.missing_taxa_height
            width = ad_tree.botR.x - ad_tree.x * self.scale - x
            ad_tree.missing_taxa_block = Block(Point(x,y),Point(x+width,y+self.missing_taxa_height))

            self.draw_missing_taxa_block(ad_tree)

    def generate_tree_name_block(self,ad_tree):
        x = ad_tree.topL.x
        y = ad_tree.topL.y - TREE_NAME_LABEL_HEIGHT - DEFAULT_PADDING_BETWEEN_BLOCK

        ad_tree.tree_name_block = Block(Point(x, y), Point(x + self.ad_width, y + TREE_NAME_LABEL_HEIGHT +
                                                           DEFAULT_PADDING_BETWEEN_BLOCK))

    def draw_indv_block(self,block,ad_tree,canvas=None):
        block.width = DEFAULT_INDV_BLOCK_WIDTH * self.scale
        if block.topL.x + block.width >= ad_tree.botR.x:
            block.botR.x = ad_tree.botR.x - 5
            block.width = block.botR.x - block.topL.x

        color = BLANK if block.type == INDIVIDUAL_BLOCK else block.color

        self.draw_exact_match_block(block.topL.x, block.topL.y, block.width, block.height,
                                    block_color=color, frame_color=GREY,
                                    layer_index=ad_tree.located_layer,canvas=canvas)

    def draw_subtree_block(self,block,ad_tree,scale_tmp,canvas=None,cluster=False):
        block.width = block.botR.x - block.topL.x

        if block.is_subtree_duplicate:
            self.draw_duplicate_subtree_block(block,canvas=canvas)
        else:
            # Draw subtree block
            if block.exact_match:
                self.draw_exact_match_block(block.topL.x, block.topL.y, block.width, block.height,
                                            block_color=block.color, frame_color=BLACK,
                                            layer_index=ad_tree.located_layer,canvas=canvas)
            else:
                self.draw_inexact_match_block(block,canvas=canvas)

        if block.height > 10 and block.width > 17:
            self.write_ad_block_label(block,canvas=canvas)

        if block.nested_tree:
            block.nested_tree.located_layer = ad_tree.located_layer
            block.nested_tree.set_nested_tree_size(block)
            self.draw_ad_tree(ad_tree=block.nested_tree,canvas=canvas,cluster=cluster)

    def draw_missing_taxa_block(self,ad_tree):
        block = ad_tree.missing_taxa_block
        self.draw_dotted_frame(block.topL.x, block.topL.y, block.width, block.height,
                               color=BLACK, line_dash_index=1, layer_index=ad_tree.located_layer)

        layer_index = ad_tree.located_layer
        scale_tmp = self.scale
        if scale_tmp < 0.7:
            scale_tmp = 0.7

        # Draw segment
        if len(ad_tree.missing_taxa_segment) != 0:
            segment_height = block.height / len(ad_tree.missing_taxa_segment)
        else:
            segment_height = block.height

        for index,segment in enumerate(ad_tree.missing_taxa_segment):
            segment_y = block.topL.y + segment_height * index
            segment_width = self.get_color_segment_width(block.width,len(segment.missing_taxa_list),len(ad_tree.missing_taxa_list))

            self.draw_rec(block.topL.x, segment_y,segment_width,segment_height,segment.subtree.color,layer_index)

        if self.view == AD_INDIVIDUAL:
        # Write missing taxa count
            self[layer_index].font = f'{12 * scale_tmp}px Times New Romans'
            # Write subtree label on top left corner
            self[layer_index].fill_style = BLACK
            self[layer_index].text_align = RIGHT
            self[layer_index].fill_text(str(len(ad_tree.missing_taxa_list)), block.botR.x - 8,
                                        block.topL.y + (12 * scale_tmp))
            self[layer_index].text_align = LEFT

    def get_color_segment_width(self,width,segment,block):
        # width : total width
        # segment : taxa count in segment
        # block : total taxa count
        return (segment/block) * width

    def get_color_segment_height(self,height,segment,block):
        # height : total height
        # segment : taxa count in segment
        # block : total taxa count
        return (segment/block) * height

    def draw_duplicate_subtree_block(self,block,canvas=None):
        height = 0
        width = 0
        x = block.topL.x
        y = block.topL.y
        for i in range(len(block.belong_subtree)):
            if i == 0:
                side = UP
            elif i == (len(block.belong_subtree) - 1):
                side = DOWN
            else:
                side = MIDDLE

            if height == 0:
                height = self.get_color_segment_height(block.height,block.duplicate_subtree_taxa_count[i],block.subtree_taxa_count)
            else:
                height = block.height - height


            if block.exact_match[i]:
                self.draw_rec(x, y, block.width, height,
                              color=block.color[i],
                              layer_index=block.belong_ad_tree.located_layer,canvas=canvas)
                self.draw_rec_edge_absent(side=side,exact_match=True,x=x,y=y,width=block.width,height=height,
                                          layer_index=block.belong_ad_tree.located_layer,canvas=canvas)
            else:
                subtree_leaf = len(block.belong_subtree[i].leaf_set)

                width = self.get_color_segment_width(block.width,block.duplicate_subtree_taxa_count[i],
                                                     len(block.taxa_list))
                self.draw_rec(x, y, width, height,
                              color=block.color[i],
                              layer_index=block.belong_ad_tree.located_layer,canvas=canvas)
                self.draw_rec_edge_absent(side=side,exact_match=False,x=x,y=y,width=block.width,height=height,layer_index=block.belong_ad_tree.located_layer,canvas=canvas)

            y += height

    def draw_rec_edge_absent(self,side,exact_match,x,y,width,height,layer_index,canvas=None):
        left_vertical = [Point(x,y),Point(x,y + height)]
        right_vertical = [Point(x + width,y),Point(x + width,y + height)]

        if side == UP:
            horizontal = [Point(x,y),Point(x + width,y)]
        elif side == DOWN:
            horizontal = [Point(x,y + height),Point(x + width,y + height)]
        else:
            horizontal = None

        if exact_match:
            self.draw_solid_line(left_vertical[0],left_vertical[1],BLACK,layer_index=layer_index,canvas=canvas)
            self.draw_solid_line(right_vertical[0],right_vertical[1],BLACK,layer_index=layer_index,canvas=canvas)
            if horizontal:
                self.draw_solid_line(horizontal[0],horizontal[1],BLACK,layer_index=layer_index,canvas=canvas)
        else:
            self.draw_dotted_line(left_vertical[0],left_vertical[1],BLACK, line_dash_index=1,layer_index=layer_index,canvas=canvas)
            self.draw_dotted_line(right_vertical[0],right_vertical[1],BLACK, line_dash_index=1,layer_index=layer_index,canvas=canvas)
            if horizontal:
                self.draw_dotted_line(horizontal[0],horizontal[1],BLACK, line_dash_index=1,layer_index=layer_index,canvas=canvas)

    def draw_root_internal(self,node,layer_index,canvas=None):
        first_child = node.children[0]
        last_child = node.children[-1]

        # Check whether branch length of two children is different
        x = first_child.branch_head.x if first_child.branch_head.x < last_child.branch_head.x else last_child.branch_head.x
        node.unify_children_branches(x)
        # first_child.branch_head.x = x
        # last_child.branch_head.x = x

        y = first_child.branch_head.y + ((last_child.branch_head.y - first_child.branch_head.y) / 2)

        scale_tmp = 0.7 if self.scale < 0.7 else self.scale
        node.construct_branch(x, y, scale_tmp)

        # Draw vertical branch linking two children
        self.draw_solid_line(Point(x, first_child.branch_head.y), Point(x, last_child.branch_head.y), BLACK,
                             layer_index,canvas=canvas)

        if node.type == ROOT:
            # Draw horizontal branch
            self.draw_solid_line(node.branch_head, node.branch_tail, BLACK, layer_index,canvas=canvas)

        # Draw children's branches
        for child in node.children:
            if child.is_elided:
                self.draw_elided_branch(child.branch_head, child.branch_tail, BLACK, layer_index,canvas=canvas)
            else:
                self.draw_solid_line(child.branch_head, child.branch_tail, BLACK, layer_index,canvas=canvas)

    def write_tree_name(self,ad_tree,extend = False,x=None):
        tree_name = str(ad_tree.tc_tree.id) + " : " + ad_tree.tc_tree.name

        if not extend and len(tree_name) > TREE_NAME_MAX_LENGTH:
            tree_name = tree_name[:TREE_NAME_MAX_LENGTH] + "..."

        if extend:
            layer = self.TREE_NAME_LAYER
        else:
            layer = ad_tree.located_layer

        self[layer].fill_style = BLACK
        self[layer].font = f'{16 * self.scale}px Times New Romans'
        if extend:
            self[layer].fill_text(tree_name, x, ad_tree.topL.y - DEFAULT_PADDING_BETWEEN_BLOCK)
        else:
            self[layer].fill_text(tree_name, ad_tree.topL.x, ad_tree.topL.y - DEFAULT_PADDING_BETWEEN_BLOCK)

    def ad_tree_filter(self,x,y):
        for ad_tree in self.ad_list:
            if ad_tree.check_in_range(x,y):
                return ad_tree

    def mouse_clicked(self, x, y):
        # Find the clicked ad_tree
        tree_selected = self.ad_tree_filter(x,y)

        # self.draw_rec(x,y,10,10,RED,-1)
        if not tree_selected:
            return
        elif tree_selected == self.tree_selected:
            self.tree_selected = None
            self[self.TREE_SELECTED_BLOCK_LAYER].clear()
        elif tree_selected != self.tree_selected:
            self.tree_selected = tree_selected
            self[self.TREE_SELECTED_BLOCK_LAYER].clear()
            self[self.TREE_SELECTED_BLOCK_LAYER].flush()

            x_padding = 10

            if self.show_tree_name:
                y_padding = 25
            else:
                y_padding = 10


            self.draw_rec(tree_selected.topL.x - x_padding, tree_selected.topL.y - y_padding , tree_selected.width +
                        2 * x_padding,
                          tree_selected.height + y_padding + x_padding,
                          LIGHT_GREY,
                          self.TREE_SELECTED_BLOCK_LAYER)

        # Testing
        # self.draw_rec(tree_selected.topL.x,tree_selected.topL.y,tree_selected.width,tree_selected.height,GREY,-1)

    def mouse_hover(self, x, y):
        ad_tree_hovered = self.find_tree_name_block(x,y)

        if not ad_tree_hovered:
            self.tree_hover = ad_tree_hovered
            self[self.TREE_NAME_LAYER].clear()
            self[self.TREE_NAME_LAYER].flush()

        if ad_tree_hovered != self.tree_hover:
            self[self.TREE_NAME_LAYER].clear()
            self[self.TREE_NAME_LAYER].flush()

        if ad_tree_hovered and ad_tree_hovered != self.tree_hover:
            self.tree_hover = ad_tree_hovered
            block = ad_tree_hovered.tree_name_block


            tree_name = str(ad_tree_hovered.tc_tree.id) + " : " + ad_tree_hovered.tc_tree.name
            name_length = len(tree_name) * 7.5 * self.scale
            if ad_tree_hovered.index_in_row == 0:
                x = ad_tree_hovered.topL.x - ad_tree_hovered.padding
            elif ad_tree_hovered.index_in_row == (self.ad_per_row - 1):
                x = self.width - name_length - DEFAULT_PADDING_BETWEEN_AD - ad_tree_hovered.padding * 2
            else:
                ad_middle_x = ad_tree_hovered.topL.x + ad_tree_hovered.width/2
                x = ad_middle_x - (name_length / 2)

            self.draw_rec(x, block.topL.y, name_length , block.height, TREE_NAME_BG, self.TREE_NAME_LAYER)
            self.draw_frame(x, block.topL.y, name_length , block.height, BLACK, self.TREE_NAME_LAYER)
            self.write_tree_name(ad_tree_hovered, extend=True,x= x + 5)

    def find_tree_name_block(self,x,y):
        for ad_tree in self.ad_list:
            if ad_tree.tree_name_block.check_in_range(Point(x,y)):
                return ad_tree

    def write_ad_block_label(self,block,canvas=None):
        if not canvas and self.view == AD_INDIVIDUAL:
            layer_index = block.belong_ad_tree.located_layer
            canvas = self[layer_index]
        elif not canvas:
            layer_index = self.AD_LAYER
            canvas = self[layer_index]

        if block.is_subtree_duplicate:
            subtree_label = ""
            taxa_cnt = ""
            for index,subtree in enumerate(block.belong_subtree):
                subtree_label += subtree.label + '&'
            subtree_label = subtree_label[:-1]
        else:
            subtree_label = block.belong_subtree.label

        taxa_cnt = str(len(block.taxa_list))
        label_width = (len(taxa_cnt) * 12)

        subtree_label_width = (len(subtree_label) * 9)

        if label_width + subtree_label_width > block.width:
            return

        scale_tmp = self.scale

        if self.scale <= 0.7:
            scale_tmp = 0.8

        canvas.font = f'{12 * scale_tmp}px Times New Romans'

        # Write subtree label on top left corner
        canvas.fill_style = BLACK

        canvas.fill_text(subtree_label, block.topL.x + 5,block.topL.y + (12 * scale_tmp))

        if block.width > 28 and not self.view == AD_CLUSTER :
            # Write taxa count on top right corner
            canvas.fill_text(taxa_cnt, block.botR.x - label_width, block.topL.y + (12 * scale_tmp))

    def individual_block(self,subtree,type=INDIVIDUAL_BLOCK):
        #  Default size : 30x7
        #  Default : topL=Point(x,y),botR=Point(x + 30,y + 7)

        block =  AD_Block(width=DEFAULT_INDV_BLOCK_WIDTH * self.scale, height=DEFAULT_INDV_BLOCK_HEIGHT * self.scale, subtree=subtree, ad_tree=self.ad_constructing, type=type)
        return block

    def ad_subtree_block(self,subtree,subtree_duplicate):
        # AD_Block height depends on node amount, width depends on x_level and tree_width
        block =  AD_Block(width=None,height=None,subtree=subtree,ad_tree=self.ad_constructing,type=SUBTREE_BLOCK)
        block.is_subtree_duplicate = subtree_duplicate
        return block


    #------------------- Cluster ad view ----------------------#

    def cluster_ad(self):
        # 1.rt ad
        # 2.all ad_to_newick
        #   - compare with existing newick and record has how many different newick, as well as record this cluster has how many tree and tree_list, and see whether this cluster is cluster of rt as well
        # 3. draw_cluster ( ad just save the first ad tree)

        # Attribute initialization
        self.rt_ad = None
        self.topology_list = []
        self.ad_from_individual = False
        self.cluster_ad_list = []

        # self.prepare_subtree_list()

        self.generate_ad_list()
        # Check if has individual_ad_canvas
        # self.ad_individual_canvas = self.adPy.ad_individual_canvas
        # if self.ad_individual_canvas:
        #     if self.parameter_from_individual_ad:
        #         self.get_parameter_from_individual_canvas()
        #
        #         if self.ad_individual_canvas.max_ad and self.ad_individual_canvas.max_ad < len(self.tc):
        #             self.generate_ad_list()
        #         else:
        #             self.ad_list = self.ad_individual_canvas.ad_list
        #             self.ad_from_individual = True
        #
        #     elif self.check_ad_parameter(self.ad_individual_canvas) and not (self.ad_individual_canvas.max_ad and self.ad_individual_canvas.max_ad < len(self.tc)):
        #             self.ad_list = self.ad_individual_canvas.ad_list
        #             self.ad_from_individual = True
        #     else:
        #         self.generate_ad_list()
        # else:
        #     self.generate_ad_list()

        # Convert rt_ad to string
        # Convert ad_tree to string and save as AD_Topology

        self.cluster_from_ad_list()
        # self.generate_topology_canvas()

        self.construct_cluster_ad_canvas()
        # Testing
        # self.get_ad_string_list()
        # self.get_topology_list()

    def topology_list_to_ad_list(self):
        self.topology_list = sorted(self.topology_list, key=lambda x: x.tree_count, reverse=True)
        self.cluster_ad_list = []
        for index,topology in enumerate(self.topology_list):
            if topology.agree_rt:
                self.cluster_agree_rt = index
            self.cluster_ad_list.append(topology.sample_ad_tree)

    def cluster_from_ad_list(self):
        rt_topology = self.adPy.rt.topology_string
        self.topology_list = []
        for ad_tree in self.ad_list:
            self.inexact_block = []
            ad_string = ad_tree.ad_to_string(self,differentiate_inexact_match=self.differentiate_inexact_match)
            if ad_tree.missing_taxa_list is not None:
                ad_string = "[MISSING]" + ad_string
            topology_exist = self.check_topology_exist(ad_string)
            if topology_exist:
                if self.inexact_block:
                    for subtree_label in self.inexact_block:
                        topology_exist.set_block_inexact_match(subtree_label)
                topology_exist.add_ad(ad_tree)
                ad_tree.topology = topology_exist
            else:
                new_topology = AD_Topology(ad_string=ad_string, sample_ad_tree=ad_tree)

                if self.differentiate_inexact_match and "[inexact]" not in ad_string:
                    ad_string_tmp = ad_string.replace("[exact]", "")

                    if ad_string_tmp == rt_topology:
                        new_topology.agree_rt = True
                elif not self.differentiate_inexact_match:
                    ad_string_tmp = ad_string.replace("[inexact]", "")
                    ad_string_tmp = ad_string_tmp.replace("[exact]", "")

                    if ad_string_tmp == rt_topology:
                        new_topology.agree_rt = True

                self.topology_list.append(new_topology)
                ad_tree.topology = new_topology

    def construct_cluster_ad_canvas(self):
        self.max_ad = len(self.topology_list)
        ad_row = math.ceil(self.max_ad / self.ad_per_row)

        canvas_height = (ad_row * DEFAULT_AD_HEIGHT * self.scale) + (2 * DEFAULT_PADDING_BETWEEN_AD) + (
                (ad_row - 1) * DEFAULT_PADDING_BETWEEN_AD)
        self.height = canvas_height

        # Convert topology_list to ad_list
        self.topology_list_to_ad_list()

        result = self.preprocess_ad_tree(self.cluster_ad_list)

        if result == SPACE_INSUFFICIENT:
            self.adjust_display()
            return

        self.check_canvas_height(self.cluster_ad_list)

        if self.scale_alter:
            self.display_content_exceed_error()

        self.draw_ad_tree_rec(self.cluster_ad_list)

        for index,ad_tree in enumerate(self.cluster_ad_list):
            self.draw_cluster_number_bar(ad_tree,self.topology_list[index].tree_count,agree_rt=self.topology_list[
                index].agree_rt)
            self.draw_ad_tree(ad_tree)
            self[self.AD_LAYER].flush()

    def generate_topology_canvas(self):
        for index,topology in enumerate(self.topology_list):
            tmp_canvas = Canvas(width=self.ad_width,height=self.ad_height)
            self.draw_ad_tree(topology.sample_ad_tree,canvas=tmp_canvas)
            topology.ad_canvas = tmp_canvas
            self.topology_canvas_list.append(tmp_canvas)

    def check_topology_exist(self,ad_string):
        for topology in self.topology_list:
            if ad_string == topology.ad_string:
                return topology
        return False

    def draw_cluster_number_bar(self,ad_tree,tree_count,canvas=None,agree_rt=False):
        # Calculate position
        frame_width = ad_tree.width
        height = CLUSTER_NUMBER_BAR_HEIGHT * self.scale
        if canvas:
            x = 5
            y = 10
        else:
            x = ad_tree.topL.x
            y = ad_tree.topL.y - DEFAULT_PADDING_BETWEEN_BLOCK - height

        color_segment_width = self.get_color_segment_width(width=frame_width,segment=tree_count,block=len(self.tc))

        # Draw bar
        self.draw_rec(x, y, color_segment_width, height, LIGHT_GREY, layer_index=ad_tree.located_layer,canvas=canvas)
        self.draw_frame(x, y, frame_width, height, BLACK, layer_index=ad_tree.located_layer, canvas=canvas)

        if not canvas:
            canvas = self[self.AD_LAYER]
        # Write tree count
        if height > 5 and frame_width > 15:
            canvas.font = f'{12 * self.scale}px Times New Romans'
            canvas.fill_style = BLACK
            canvas.fill_text(str(tree_count), x + 5, y + (12 * self.scale))

        if agree_rt:
            canvas.fill_style = RED
            canvas.fill_text("R", x + ad_tree.width - 10, y + CLUSTER_NUMBER_BAR_HEIGHT - 3)
            canvas.fill_style = BLACK

    def draw_exact_match_block(self,x,y,width,height,block_color,frame_color,layer_index = -1,canvas=None):
        self.draw_rec(x,y,width,height,block_color,layer_index=layer_index,canvas=canvas)
        self.draw_frame(x, y, width, height, frame_color, layer_index=layer_index,canvas=canvas)

    def draw_inexact_match_block(self,block,canvas=None):
        # block.topL.x, block.topL.y, block.width, block.height,block_color=block.color, frame_color=BLACK,layer_index=self.AD_LAYER
        x = block.topL.x
        y = block.topL.y
        height = block.height
        width = block.width

        if self.show_block_proportional:
            segment_width = self.get_color_segment_width(block.width,block.subtree_taxa_count,len(block.taxa_list))
        else:
            segment_width = width

        self.draw_rec(x, y, segment_width, height, block.color, layer_index=block.belong_ad_tree.located_layer,
                      canvas=canvas)
        self.draw_dotted_line(Point(x,y),Point(x+width,y),BLACK,line_dash_index=1,
                              layer_index=block.belong_ad_tree.located_layer,canvas=canvas)
        self.draw_dotted_line(Point(x+width,y), Point(x + width, y+height), BLACK,line_dash_index=1,
                              layer_index=block.belong_ad_tree.located_layer,canvas=canvas)
        self.draw_dotted_line(Point(x, y), Point(x, y+height), BLACK, line_dash_index=1,
                              layer_index=block.belong_ad_tree.located_layer,canvas=canvas)
        self.draw_dotted_line(Point(x, y+height), Point(x + width, y + height), BLACK,line_dash_index=1,
                              layer_index=block.belong_ad_tree.located_layer,canvas=canvas)

    # General function
    def draw_ad_tree_rec(self,ad_list):
        for index,ad_tree in enumerate(ad_list):
            self.draw_frame(ad_tree.topL.x,ad_tree.topL.y,ad_tree.width,ad_tree.height,GREY,ad_tree.located_layer)

    def display_content_exceed_error(self):
        print("<Warning> : Excessive display content.")
        print(f"Scale has been automatically adjusted to {self.scale}.")

    def check_ad_parameter(self,canvas):
        if (self.context_level == canvas.context_level and
                self.escape_taxa_as_context_block == canvas.escape_taxa_as_context_block and self.subtree_independent ==
                canvas.subtree_independent and self.show_block_proportional == canvas.show_block_proportional and
                self.scale == canvas.scale):
            return True
        return False

    def get_ad_string_list(self):
        for index,ad_string in enumerate(self.ad_string_list):
            print("AD " + str(index) ,end=" :")
            print(ad_string)

    def get_topology_list(self):
        for ad_topology in self.topology_list:
            print(ad_topology.ad_string,end=" :")
            print(ad_topology.tree_count)

    def get_parameter_from_individual_canvas(self):
        self.scale = self.ad_individual_canvas.scale
        self.context_level = self.ad_individual_canvas.context_level
        self.escape_taxa_as_context_block = self.ad_individual_canvas.escape_taxa_as_context_block
        self.show_block_proportional = self.ad_individual_canvas.show_block_proportional
        self.subtree_independent = self.subtree_independent

