from myTree import *

mytree = myTree(treefile = "Data/69species/astral.FAA.trim50genes.final.tre")
# mytree.add_tree_collection(treefile = "Data/69species/tc.tre")

# 假设你有两个节点
node1 = mytree.rt.find_node_with_taxon_label("Uronema sp")
node2 = mytree.rt.find_node_with_taxon_label("Zamia vazquezii")

# distance_to_root = node2.level()
# print(f"{node2.label} 距离根节点的代数距离为: {distance_to_root}")
# for node in mytree.rt.seed_node.adjacent_nodes():
#     print(node)

root = mytree.rt.seed_node

print()

# for node in mytree.rt.postorder_node_iter():
#     print(node)
# mytree.add_tree_collection(treefile = "Data/69species/tc.tre",namefile="Data/69species/names.txt")
#
#
# print(node.corr)
# print(mytree.rt.taxa_list)

# for tc_tree in mytree.tc:
#     print(tc_tree.index)
#     print(tc_tree.name)
#     print(tc_tree.missing)

# for rt_node in rt_tree:
# 	for index,tc_tree in enumerate(tc):
#         targetset = rt_branch.taxa(set) - tc_tree.missing(set)   #branch.head_node.leaf_nodes()
#         for tc_node in tc_tree:  # post order
#             similarity = getSimilarity()
#             if similarity > rt_node.similarity:
#                 rt_node.corr[index] = tc_node
#                 break



# node.child_nodes() : only next generation nodes
# node.leaf_nodes() : get all leaf nodes
# .adjacent_nodes():