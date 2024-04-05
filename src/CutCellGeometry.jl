module CutCellGeometry
using CutCellMesh
using ForwardDiff
using Plots
Plots.default(show = true)
using LinearAlgebra
using Roots

# Include files
include("body.jl")
include("plot.jl")
include("utils.jl")
include("interface.jl")

# Export functions
export SignedDistanceFunction, evaluate_sdf, ⊔, ⊓, ⊖, complement
export plot_sdf, plot_cut_cells_levelset_intersections_and_midpoints
export compute_velocity, gradient_Phi, hessian_Phi, normal, tangent_vector, curvature
export evaluate_levelset, get_cut_cells, get_intersection_points, get_segment_midpoints

end # module CutCellGeometry


