module CutCellGeometry
using CutCellMesh
using ForwardDiff
using Plots
Plots.default(show = true)
# Include files
include("body.jl")
include("plot.jl")

# Export functions
export SignedDistanceFunction, evaluate_sdf, ⊔, ⊓, ⊖, complement
export plot_sdf

end # module CutCellGeometry


