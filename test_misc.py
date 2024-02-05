def generate_coordinates_linear(n, left_boundary=0, right_boundary=100):
    # List to store node coordinates
    coordinates = []

    # Calculate coordinates for the nodes along the x-axis
    for i in range(n):
        x = left_boundary + (i / (n-1)) * (right_boundary - left_boundary)
        y = 0
        coordinates.append((x, y))

    return coordinates

# Example with 5 nodes
node_count = 5
node_coordinates = generate_coordinates_linear(node_count)

print("Node Coordinates:")
for i, (x, y) in enumerate(node_coordinates):
    print(f"Node {i + 1}: ({x:.2f}, {y:.2f})")
