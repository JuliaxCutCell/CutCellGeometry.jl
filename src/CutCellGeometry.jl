module CutCellGeometry

using CutCellMesh
using ForwardDiff
using Plots
Plots.default(show = true)
using LinearAlgebra
using Roots
using CartesianGeometry
using SparseArrays

# Include files
include("body.jl")
include("plot.jl")
include("utils.jl")
include("interface.jl")
include("vofi.jl")
include("primitives.jl")

# Export functions
export SignedDistanceFunction, evaluate_sdf, ⊔, ⊓, ⊖, complement
export plot_sdf, plot_cut_cells_levelset_intersections_and_midpoints
export compute_velocity, gradient_Phi, hessian_Phi, normal, tangent_vector, curvature
export evaluate_levelset, get_cut_cells, get_intersection_points, get_segment_midpoints
export calculate_first_order_moments, calculate_second_order_moments, first_aperture_1D, first_aperture_2D, first_aperture_3D, second_aperture_1D, second_aperture_2D, second_aperture_3D

# Exemple d'utilisation
# Cas 2D
grid = CartesianGrid((10, 10), (0.1, 0.1))
mesh = generate_mesh(grid)
domain = ((minimum(mesh[1]), maximum(mesh[1])), (minimum(mesh[2]), maximum(mesh[2])))

# Fonction de distance signée
cercl = cercle((0.5,0.5), 0.25)
square = rectangle((0.25, 0.25), (0.75, 0.75))
identity_transform = (x, y, t) -> (x, y)
circle_sdf = SignedDistanceFunction(cercl, identity_transform, domain, false)
square_sdf = SignedDistanceFunction(square, identity_transform, domain, false)

plot_sdf(circle_sdf, domain, 0.0)
readline()

circle_sdf = translate(circle_sdf, 0.1, 0.2, 0.3)
plot_sdf(circle_sdf, domain, 0.0)
readline()

# Moments
V, bary, As = calculate_first_order_moments(circle_sdf.sdf_function, mesh)
Ws, Bs = calculate_second_order_moments(circle_sdf.sdf_function, mesh, bary)

V, bary, As = calculate_first_order_moments(square_sdf.sdf_function, mesh)
Ws, Bs = calculate_second_order_moments(square_sdf.sdf_function, mesh, bary)

@show V
@show bary
@show As

# Cas 3D
grid = CartesianGrid((10, 10, 10), (0.1, 0.1, 0.1))
mesh = generate_mesh(grid)
domain = ((minimum(mesh[1]), maximum(mesh[1])), (minimum(mesh[2]), maximum(mesh[2])), (minimum(mesh[3]), maximum(mesh[3])))

# Fonction de distance signée
spher = sphere((0.5, 0.5, 0.5), 0.25)
cuboi = cuboid((0.25, 0.25, 0.25), (0.75, 0.75, 0.75))
identity_transform = (x, y, z, t) -> (x, y, z)
sphere_sdf = SignedDistanceFunction(spher, identity_transform, domain, false)
cuboid_sdf = SignedDistanceFunction(cuboi, identity_transform, domain, false)


# Moments
V, bary, As = calculate_first_order_moments(sphere_sdf.sdf_function, mesh)
Ws, Bs = calculate_second_order_moments(sphere_sdf.sdf_function, mesh, bary)

@show Ws
V, bary, As = calculate_first_order_moments(cuboid_sdf.sdf_function, mesh)
Ws, Bs = calculate_second_order_moments(cuboid_sdf.sdf_function, mesh, bary)

@show V
@show bary
end # module CutCellGeometry


