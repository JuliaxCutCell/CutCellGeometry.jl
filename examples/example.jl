using CutCellMesh # Mesh generation
using CutCellGeometry # Geometry definition

# CartesianGrid 
grid = CartesianGrid((10, 10), (0.1, 0.1))
mesh = generate_mesh(grid)
domain = ((minimum(mesh[1]), maximum(mesh[1])), (minimum(mesh[2]), maximum(mesh[2])))

# SignedDistanceFunction
circle_function = (x, y, _=0) -> sqrt(x^2 + y^2) - 0.5
square_function = (x, y, _=0) -> max(abs(x), abs(y)) - 0.5
identity_transform = (x, y, t) -> (x, y)

# Create Geometry
circle_sdf = SignedDistanceFunction(circle_function, identity_transform, domain, false)
square_sdf = SignedDistanceFunction(square_function, identity_transform, domain, false)

# Union of two signed distance functions
union_sdf = circle_sdf ⊔ square_sdf

# Transformation function
move_transform = (x, y, t) -> ((x-0.25) + 2*t, (y-0.25))

# Move the union of two signed distance functions
moving_sdf = SignedDistanceFunction(union_sdf.sdf_function, move_transform, union_sdf.domain, true)

# Print the normal vector at a given point
println("Normal au point (0.5, 0.5) pour le cercle", normal(circle_sdf.sdf_function, (0.5, 0.5)))

# Calculer le vecteur tangentiel
tangent = tangent_vector(normal(circle_sdf.sdf_function, (0.5, 0.5)))
println("Vecteur tangentiel au point (0.5, 0.5) pour le cercle", tangent)

# Calculer la courbure pour le cercle
println("Courbure pour le cercle au point (0.5, 0.5)", curvature(circle_function.sdf_function, (0.5, 0.5)))

# Plot
t=0.0
plot_sdf(moving_sdf, domain, t)
readline()

# Geometric Moments
V, v_diag, bary, ax_diag, ay_diag = calculate_first_order_moments(circle_sdf.sdf_function, mesh)
w_diag, bx_diag, by_diag, border_cells_wx, border_cells_wy = calculate_second_order_moments(circle_sdf.sdf_function, mesh, bary)