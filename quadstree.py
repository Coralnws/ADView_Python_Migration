import quads

tree = quads.QuadTree(
(50, 50),  # The center point
 100,  # The width
100,  # The height
 capacity=4
)
tree.insert((1,2))
tree.insert((10,20))
tree.insert((50,50))
tree.insert((20,40))
tree.insert((70,10))
tree.insert((3,88))
print(tree.nearest_neighbors((100,100),1))
