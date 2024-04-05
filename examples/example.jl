using CutCellMesh # Mesh generation
using CutCellGeometry # Geometry definition

# CartesianGrid 
grid = CartesianGrid((10, 10), (0.1, 0.1))
mesh = generate_mesh(grid)
domain = ((minimum(mesh[1]), maximum(mesh[1])), (minimum(mesh[2]), maximum(mesh[2])))

# SignedDistanceFunction
circle_function = (x, y) -> sqrt(x^2 + y^2) - 1.0
square_function = (x, y) -> max(abs(x), abs(y)) - 1.0
identity_transform = (x, y, t) -> (x, y)

# Create Geometry
circle_sdf = SignedDistanceFunction(circle_function, identity_transform, domain, false)
square_sdf = SignedDistanceFunction(square_function, identity_transform, domain, false)

# Union of two signed distance functions
union_sdf = circle_sdf ⊔ square_sdf

# Transformation function
move_transform = (x, y, t) -> (x + t, y)

# Move the union of two signed distance functions
moving_sdf = SignedDistanceFunction(union_sdf.sdf_function, move_transform, union_sdf.domain, true)

t=0.0
plot_sdf(moving_sdf, domain, t)
readline()