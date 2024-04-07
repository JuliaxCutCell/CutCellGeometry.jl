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
        Plots.plot(x_range, x -> evaluate_sdf(sdf, t, x), color=:red, linewidth=2)
    elseif dim == 2
        x_range = LinRange(domain[1][1], domain[1][2], 100)
        y_range = LinRange(domain[2][1], domain[2][2], 100)
        Plots.contour(x_range, y_range, (x, y) -> evaluate_sdf(sdf, t, x, y), levels=[0.0], color=:red, linewidth=2)
    elseif dim == 3
        x_range = collect(LinRange(domain[1][1], domain[1][2], 100))
        y_range = collect(LinRange(domain[2][1], domain[2][2], 100))
        z_range = collect(LinRange(domain[3][1], domain[3][2], 100))
        S = [evaluate_sdf(sdf, t, x, y, z) for x in x_range, y in y_range, z in z_range]
    else
        println("Only 1D, 2D, and 3D domains are supported.")
    end
end

"""
    plot_cut_cells_levelset_intersections_and_midpoints(cut_cells, values, intersection_points, midpoints)

Plot the cut cells, level set intersections, and midpoints.

# Arguments
- `cut_cells`: An array of cut cells.
- `values`: The level set values.
- `intersection_points`: The intersection points.
- `midpoints`: The midpoints.

# Details
- The function creates a new plot and adds the level set as a contour line.
- It determines the dimension based on the first element of `cut_cells`.
- For each cut cell, it adds the corresponding shape to the plot.
- It adds the intersection points and midpoints to the plot.
- Finally, it displays the plot.

# Example
"""
function plot_cut_cells_levelset_intersections_and_midpoints(cut_cells, values, intersection_points, midpoints)
    # Créer un nouveau tracé
    plt = plot()

    # Déterminer la dimension à partir du premier élément de cut_cells
    dimension = length(Tuple(cut_cells[1]))

    if dimension == 1
        # En 1D, tracer une ligne à la valeur de zéro
        plot!(values, line=(:red, 2), label="Level Set")
    elseif dimension == 2
        # Ajouter la Level Set comme une ligne de contour
        contour!(values, levels=[0], color=:red)
    else
        println("3D plotting not supported")
    end
    # Ajouter chaque cellule coupée au tracé
    for cell in cut_cells
        if dimension == 1
            # En 1D, une cellule est un segment de ligne
            plot!(plt, [cell[1], cell[1]+1], color=:blue)
        elseif dimension == 2
            # En 2D, une cellule est un carré
            x = cell[1]
            y = cell[2]
            plot!(plt, [x, x+1, x+1, x, x], [y, y, y+1, y+1, y], fill = true, color = :blue, alpha=0.5)
        else
            # En 3D, une cellule est un cube (non représentable avec Plots.jl)
            println("3D plotting not supported")
        end
    end

    # Ajouter les points d'intersection au tracé
    for point in intersection_points
        if dimension == 1
            scatter!(plt, [point[1]], color=:green, markersize=4)
        elseif dimension == 2
            scatter!(plt, [point[1]], [point[2]], color=:green, markersize=4)
        else
            # En 3D, un point est un point dans l'espace (non représentable avec Plots.jl)
            println("3D plotting not supported")
        end
    end

    # Ajouter les points médians au tracé
    for point in midpoints
        if dimension == 1
            scatter!(plt, [point[1]], color=:yellow, markersize=4)
        elseif dimension == 2
            scatter!(plt, [point[1]], [point[2]], color=:yellow, markersize=4)
        else
            # En 3D, un point est un point dans l'espace (non représentable avec Plots.jl)
            println("3D plotting not supported")
        end
    end

    # Afficher le tracé
    plt
    display(plt)
    readline()

end