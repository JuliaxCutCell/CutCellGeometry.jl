"""
    plot_sdf(sdf, domain, t)

Plot the signed distance function `sdf` over the given `domain` at time `t`.

# Arguments
- `sdf`: A signed distance function.
- `domain`: The domain over which to plot the signed distance function. It should be a tuple of intervals for each dimension.
- `t`: The time at which to evaluate the signed distance function.

"""
function plot_sdf(sdf, domain, t)
    dim = length(domain)
    if dim == 1
        x_range = LinRange(domain[1][1], domain[1][2], 100)
        plot(x_range, x -> evaluate_sdf(sdf, t, x), color=:red, linewidth=2)
    elseif dim == 2
        x_range = LinRange(domain[1][1], domain[1][2], 100)
        y_range = LinRange(domain[2][1], domain[2][2], 100)
        contour(x_range, y_range, (x, y) -> evaluate_sdf(sdf, t, x, y), levels=[0.0], color=:red, linewidth=2)
    elseif dim == 3
        x_range = LinRange(domain[1][1], domain[1][2], 100)
        y_range = LinRange(domain[2][1], domain[2][2], 100)
        z_range = LinRange(domain[3][1], domain[3][2], 100)
        contour(x_range, y_range, z_range, (x, y, z) -> evaluate_sdf(sdf, t, x, y, z), levels=[0.0], color=:red, linewidth=2)
    else
        println("Dimension non prise en charge")
    end
end