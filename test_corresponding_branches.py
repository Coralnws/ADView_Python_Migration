import pytest
import dendropy
from ADViewpy import *
import string
from myUtils import *

class TestCorr:
    nwk_str1 = "((D:0.1,(H:0.5,G:0.3)E:0.2)A:0.6,(B:0.2,C:0.3,F:0.1)D:0.4)E;"
    nwk_str2 = "((D:0.1,(T:0.5,P:0.3)E:0.2)A:0.6,(B:0.2,F:0.3,L:0.1)D:0.4)E;"
    node_list = []

    def generate_node(self):
        letters = string.ascii_uppercase

        for letter in letters:
            if letter > 'G':
                break

            new_taxon = dendropy.Taxon(label=letter)
            new_node = dendropy.Node(taxon=new_taxon,label=letter)
            self.node_list.append(new_node)

            if letter == 'A' or letter == 'B' :
                new_taxon = dendropy.Taxon(label=letter)
                new_node = dendropy.Node(taxon=new_taxon, label=letter)
                self.node_list.append(new_node)

    def test_similarity(self):
        expected = 1
        actual = 1

        taxon1 = dendropy.Taxon(label="node1")
        node1 = dendropy.Node(taxon=taxon1,label="node1")
        taxon2 = dendropy.Taxon(label="node2")
        node2 = dendropy.Node(taxon=taxon2, label="node2")

        self.generate_node()

        # same node : node_list[0,1] node_list[2,3]
        node1.add_child(self.node_list[0])
        node1.add_child(self.node_list[2])
        node1.add_child(self.node_list[4])

        node2.add_child(self.node_list[1])
        node2.add_child(self.node_list[3])
        node2.add_child(self.node_list[5])

        node_taxa = [leaf.taxon.label for leaf in node1.child_nodes()]

        target_set = set(node_taxa)
        missing = set()


        node_similarity =  ADViewpy.get_similarity(target_set={'A'}, node=self.node_list[0],
                                                   tc_node=self.node_list[1],missing=missing)
        if node_similarity != 1:
            actual = 0

        node_similarity = ADViewpy.get_similarity(target_set={'B'}, node=self.node_list[2],
                                                  tc_node=self.node_list[3], missing=missing)
        if node_similarity != 1:
            actual = 0

        node_similarity = ADViewpy.get_similarity(target_set={'C'}, node=self.node_list[4],
                                                  tc_node=self.node_list[5], missing=missing)
        if node_similarity != 0:
            actual = 0

        node_similarity = ADViewpy.get_similarity(target_set=target_set, node=node1,
                                                  tc_node=node2, missing=missing)
        if node_similarity != 0.5:
            actual = 0

        assert expected == actual

    def test_corr_same_tree(self):
        expected = 1
        actual = 1

        tree1 = dendropy.Tree.get(data=self.nwk_str1,schema='newick')
        tree2 = dendropy.Tree.get(data=self.nwk_str1, schema='newick')

        tree1.taxa = {leaf.taxon.label for leaf in tree1.seed_node.leaf_nodes()}
        tree2.taxa = {leaf.taxon.label for leaf in tree2.seed_node.leaf_nodes()}

        tree1.missing = tree2.taxa - tree1.taxa
        tree2.missing = tree1.taxa - tree2.taxa

        for node in tree1.postorder_node_iter():
            node_taxa = [leaf.taxon.label for leaf in node.leaf_nodes()]
            node.corr = None
            node.corr_similarity = 0

            target_set = set(node_taxa) - tree2.missing
            missing = tree2.missing - set(node_taxa)

            for tc_node in tree2.postorder_node_iter():
                similarity = ADViewpy.get_similarity(target_set=target_set,node=node,tc_node=tc_node,missing=missing)

                if similarity > node.corr_similarity:
                    node.corr = tc_node
                    node.corr_similarity = similarity

        for node in tree1.postorder_node_iter():
            corr = node.corr
            if node.is_leaf():
                if corr.taxon.label != node.taxon.label:
                    actual = 0
            else:
                tree1_node_taxa = [leaf.taxon.label for leaf in node.leaf_nodes()]
                tree2_node_taxa = [leaf.taxon.label for leaf in corr.leaf_nodes()]

                if tree1_node_taxa != tree2_node_taxa:
                    actual = 0

        assert expected == actual

    def test_corr_different_tree(self):
        expected = 1
        actual = 1

        tree1 = dendropy.Tree.get(data=self.nwk_str1,schema='newick')
        tree2 = dendropy.Tree.get(data=self.nwk_str2, schema='newick')

        tree1.taxa = {leaf.taxon.label for leaf in tree1.seed_node.leaf_nodes()}
        tree2.taxa = {leaf.taxon.label for leaf in tree2.seed_node.leaf_nodes()}

        tree1.missing = tree2.taxa - tree1.taxa
        tree2.missing = tree1.taxa - tree2.taxa

        for node in tree1.postorder_node_iter():
            node_taxa = [leaf.taxon.label for leaf in node.leaf_nodes()]
            node.corr = None
            node.corr_similarity = 0

            target_set = set(node_taxa) - tree2.missing
            missing = tree2.missing - set(node_taxa)

            for tc_node in tree2.postorder_node_iter():
                similarity = ADViewpy.get_similarity(target_set=target_set,node=node,tc_node=tc_node,missing=missing)

                if similarity > node.corr_similarity:
                    node.corr = tc_node
                    node.corr_similarity = similarity

        for node in tree1.postorder_node_iter():
            corr = node.corr
            if node.is_leaf():
                if node.taxon.label == 'G':
                    if corr != None:
                        actual = 0

            else:
                tree1_node_taxa = [leaf.taxon.label for leaf in node.leaf_nodes()]
                if corr:
                    tree2_node_taxa = [leaf.taxon.label for leaf in corr.leaf_nodes()]

                if tree1_node_taxa == ['F', 'G']:
                    if corr != None:
                        actual = 0

                elif tree1_node_taxa == ['B','C','F']:
                    if tree2_node_taxa != ['B','F','L']:
                        actual = 0


        assert expected == actual

class TestAD:
    def init_adviewpy(self,rt_file=None,tc_file=None):
        if not rt_file:
            rt_file = "Data/69species/astral.FAA.trim50genes.final.tre"
        if not tc_file:
            tc_file = "Data/69species/astral.FAA.trim50genes.final.tre"

        adviewpy = init(treefile=rt_file)
        outgroup = ["Uronema sp", "Monomastix opisthostigma", "Pyramimonas parkeae", "Nephroselmis pyriformis"]
        adviewpy.set_outgroup(outgroup)
        adviewpy.add_tree_collection(treefile=tc_file)
        adviewpy.select_subtree(["spirogyra", "mougeotia"])
        adviewpy.AD(view='AD Individual', scale=1.0, sort='id', context_level=2, show_tree_name=True)

        return adviewpy

    def check_elided_branches(self):
        expected = 1
        actual = 1

        adviewpy = self.init_adviewpy()

        # Check ad tree topology
        ad_tree = adviewpy.ad_individual_canvas.ad_list[0]
        # Should only has one level
        if len(ad_tree.root.children) != 2:
            actual = 0

        first_children = ad_tree.root.children[0]
        second_children = ad_tree.root.children[1]

        if not ((first_children.node_or_block.type == INDIVIDUAL_BLOCK and second_children.node_or_block.type ==
                 SUBTREE_BLOCK) or (
                        first_children.node_or_block.type == SUBTREE_BLOCK and second_children.node_or_block.type ==
                        INDIVIDUAL_BLOCK)):
            actual = 0

        if first_children.node_or_block.type == SUBTREE_BLOCK and not first_children.node_or_block.is_elided:
            actual = 0

        if second_children.node_or_block.type == SUBTREE_BLOCK and not second_children.node_or_block.is_elided:
            actual = 0

        assert expected == actual

    def test_check_relation(self):
        expected = 1
        actual = 1

        adviewpy = self.init_adviewpy()
        # Test check_relation() and elided branch
        node_level = 6
        for i in range(node_level):
            branch_one = adviewpy.ad_individual_canvas.check_result_list[(2 * i)]
            branch_two = adviewpy.ad_individual_canvas.check_result_list[(2 * i) + 1]

            if not ((branch_one[0] == 3  and branch_two[0] == 4) or (branch_one[0] == 4  and branch_two[0] == 3)):
                actual = 0

        branch_one = adviewpy.ad_individual_canvas.check_result_list[(2 * node_level)]
        branch_two = adviewpy.ad_individual_canvas.check_result_list[(2 * node_level) + 1]

        if not ((branch_one[0] == 4 and branch_two[0] == 0 and branch_two[1].label == 'A') or (branch_one[0] == 0 and
                                                                                               branch_one[1].label ==
                                                                                               'A' and branch_two[0]
                                                                                               == 4)):
            actual = 0

        assert expected == actual

        # Check

    def test_nested_tree(self):
        expected = 1
        actual = 1

        adviewpy = self.init_adviewpy()
        adviewpy.select_subtree(["mougeotia", "cushleckae"])
        adviewpy.AD(view='AD Individual', scale=1.0, sort='id', context_level=2, show_tree_name=True, tree_id=1)

        ad_tree = adviewpy.ad_individual_canvas.ad_list[0]
        subtree_root_node = ad_tree.root.children[0]
        if subtree_root_node.node_or_block.type != SUBTREE_BLOCK:
            subtree_root_node = ad_tree.root.children[1]

        subtree_block = subtree_root_node.node_or_block
        if not subtree_block.nested_tree:
            actual = 0

        if subtree_block.nested_tree and not subtree_block.nested_tree.is_nested:
            actual = 0

        ad_tree = subtree_block.nested_tree

        if len(ad_tree.root.children) != 2:
            actual = 0

        first_children = ad_tree.root.children[0]
        second_children = ad_tree.root.children[1]

        if not ((first_children.node_or_block.type == INDIVIDUAL_BLOCK and second_children.node_or_block.type ==
                 SUBTREE_BLOCK) or (
                        first_children.node_or_block.type == SUBTREE_BLOCK and second_children.node_or_block.type ==
                        INDIVIDUAL_BLOCK)):
            actual = 0

        if first_children.node_or_block.type == SUBTREE_BLOCK and (not first_children.node_or_block.is_elided or
        first_children.node_or_block.belong_subtree.label != 'B'):
            actual = 0

        if second_children.node_or_block.type == SUBTREE_BLOCK and (not second_children.node_or_block.is_elided or
                                                                    second_children.node_or_block.belong_subtree.label != 'B'):
            actual = 0



        assert expected == actual

    def test_ad_string(self):
        expected = 1
        actual = 1

        adviewpy = self.init_adviewpy(tc_file="Data/69species/tc_test_cluster.tre")

        ad_same = adviewpy.ad_individual_canvas.ad_list[1]
        ad_different = adviewpy.ad_individual_canvas.ad_list[0]

        # Check differentiate inexact match = True
        ad_same = adviewpy.ad_individual_canvas.ad_list[1]
        ad_different = adviewpy.ad_individual_canvas.ad_list[0]

        adviewpy.rt.topology_string = adviewpy.rt.ad_tree.ad_to_string(adviewpy.ad_individual_canvas,
                                                                       differentiate_inexact_match=True)
        ad_same.ad_string = ad_same.ad_to_string(adviewpy.ad_individual_canvas,differentiate_inexact_match=True)
        ad_different.ad_string = ad_different.ad_to_string(adviewpy.ad_individual_canvas,
                                                           differentiate_inexact_match=True)


        if ad_same.ad_string != adviewpy.rt.topology_string:
            actual = 0

        if ad_same.ad_string == ad_different.ad_string:
            actual = 0

        if ad_different.ad_string == adviewpy.rt.topology_string:
            actual = 0

        # Check differentiate inexact match = False
        adviewpy.rt.topology_string = adviewpy.rt.ad_tree.ad_to_string(adviewpy.ad_individual_canvas,
                                                                       differentiate_inexact_match=False)
        ad_same.ad_string = ad_same.ad_to_string(adviewpy.ad_individual_canvas, differentiate_inexact_match=False)
        ad_different.ad_string = ad_different.ad_to_string(adviewpy.ad_individual_canvas,
                                                           differentiate_inexact_match=False)

        if ad_same.ad_string != ad_different.ad_string != adviewpy.rt.topology_string:
            actual = 0


        assert expected == actual


if __name__ == '__main__':
    pytest.main([])