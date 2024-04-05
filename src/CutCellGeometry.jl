module CutCellGeometry
using CutCellMesh
using ForwardDiff
using Plots
Plots.default(show = true)
using LinearAlgebra

# Include files
include("body.jl")
include("plot.jl")
include("utils.jl")

# Export functions
export SignedDistanceFunction, evaluate_sdf, ⊔, ⊓, ⊖, complement
export plot_sdf
export compute_velocity, gradient_Phi, hessian_Phi, normal, tangent_vector, curvature

end # module CutCellGeometry


