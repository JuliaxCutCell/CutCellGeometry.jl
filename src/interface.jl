"""
    evaluate_levelset(f::Function, mesh::Tuple{AbstractVector})

Evaluate a level set function `f` on a 1D mesh.

# Arguments
- `f::Function`: The level set function to evaluate.
- `mesh::Tuple{AbstractVector}`: The 1D mesh on which to evaluate the level set function.

# Returns
- An array containing the evaluated values of the level set function on the mesh.
"""
function evaluate_levelset(f::Function, mesh::Tuple{AbstractVector})
    x = mesh[1]
    [f(i) for i in x]
end

"""
    evaluate_levelset(f::Function, mesh::NTuple{2,AbstractVector})

Evaluate a level set function `f` on a 2D mesh.

# Arguments
- `f::Function`: The level set function to evaluate.
- `mesh::NTuple{2,AbstractVector}`: The 2D mesh on which to evaluate the level set function.

# Returns
- An array containing the evaluated values of the level set function on the mesh.
"""
function evaluate_levelset(f::Function, mesh::NTuple{2,AbstractVector})
    x, y = mesh
    [f(i, j) for i in x, j in y]
end

"""
    evaluate_levelset(f::Function, mesh::NTuple{3,AbstractVector})

Evaluate a level set function `f` on a 3D mesh.

# Arguments
- `f::Function`: The level set function to evaluate.
- `mesh::NTuple{3,AbstractVector}`: The 3D mesh on which to evaluate the level set function.

# Returns
- An array containing the evaluated values of the level set function on the mesh.
"""
function evaluate_levelset(f::Function, mesh::NTuple{3,AbstractVector})
    x, y, z = mesh
    [f(i, j, k) for i in x, j in y, k in z]
end

"""
    get_cut_cells(values::AbstractVector)

Compute the indices of the cut cells in a 1D array based on the sign change of adjacent elements.

# Arguments
- `values::AbstractVector`: The input 1D array.

# Returns
- `cut_cells`: An array of `CartesianIndex{1}` representing the indices of the cut cells.
"""
function get_cut_cells(values::AbstractVector)
    cut_cells = similar(values, CartesianIndex{1}, 0)

    for i in axes(values, 1)[begin:end-1]
        values[i] * values[i+1] < 0 &&
            push!(cut_cells, CartesianIndex(i))
    end

    cut_cells
end

"""
    get_cut_cells(values::AbstractMatrix)

Compute the indices of the cut cells in a 2D array based on the sign change of adjacent elements.

# Arguments
- `values::AbstractMatrix`: The input 2D array.

# Returns
- `cut_cells`: An array of `CartesianIndex{2}` representing the indices of the cut cells.
"""
function get_cut_cells(values::AbstractMatrix)
    cut_cells = similar(values, CartesianIndex{2}, 0)

    for j in axes(values, 2)[begin:end-1]
        for i in axes(values, 1)[begin:end-1]
            (values[i, j] * values[i+1, j] < 0 ||
             values[i, j] * values[i, j+1] < 0 ||
             values[i+1, j] * values[i+1, j+1] < 0 ||
             values[i, j+1] * values[i+1, j+1] < 0) &&
                push!(cut_cells, CartesianIndex(i, j))
        end
    end

    cut_cells
end

"""
    get_cut_cells(values::AbstractArray{<:Any,3})

Compute the indices of the cut cells in a 3D array based on the sign change of adjacent elements.

# Arguments
- `values::AbstractArray{<:Any,3}`: The input 3D array.

# Returns
- `cut_cells`: An array of `CartesianIndex{3}` representing the indices of the cut cells.
"""
function get_cut_cells(values::AbstractArray{<:Any,3})
    cut_cells = similar(values, CartesianIndex{3}, 0)

    for k in axes(values, 3)[begin:end-1]
        for j in axes(values, 2)[begin:end-1]
            for i in axes(values, 1)[begin:end-1]
                (values[i, j, k] * values[i+1, j, k] < 0 ||
                 values[i, j, k] * values[i, j+1, k] < 0 ||
                 values[i, j, k] * values[i, j, k+1] < 0 ||
                 values[i+1, j, k] * values[i+1, j+1, k] < 0 ||
                 values[i+1, j, k] * values[i+1, j, k+1] < 0 ||
                 values[i, j+1, k] * values[i, j+1, k+1] < 0 ||
                 values[i, j, k+1] * values[i+1, j, k+1] < 0 ||
                 values[i, j, k+1] * values[i, j+1, k+1] < 0 ||
                 values[i+1, j+1, k] * values[i+1, j+1, k+1] < 0) &&
                    push!(cut_cells, CartesianIndex(i, j, k))
            end
        end
    end

    cut_cells
end

"""
    get_intersection_points(values::AbstractVector{T}, cut_cells)

Compute the intersection points between a level set function and the edges of the cut cells.

# Arguments
- `values::AbstractVector{T}`: A vector containing the values of the level set function at each grid point.
- `cut_cells`: An iterable containing the indices of the cut cells.

# Returns
- `intersection_points`: A tuple of intersection points, where each point is represented as a tuple.
"""
function get_intersection_points(values::AbstractVector{T}, cut_cells) where {T}
    # Initialiser un tableau vide pour stocker les points d'intersection
    intersection_points = similar(cut_cells, Tuple{T}, 0)

    # Parcourir toutes les cellules coupées
    for index in cut_cells
        i, = Tuple(index)

        # Vérifier si la Level Set change de signe le long de cette arête
        if values[i] * values[i+1] < 0
            # Si c'est le cas, calculer le point d'intersection
            t = values[i] / (values[i] - values[i+1])
            x_intersect = i + t

            # Ajouter le point d'intersection à la liste
            push!(intersection_points, (x_intersect,))
        end
    end

    intersection_points
end

"""
    get_intersection_points(values::AbstractMatrix{T}, cut_cells)

Compute the intersection points between the level set function and the cell edges for a given set of cut cells.

# Arguments
- `values::AbstractMatrix{T}`: A matrix representing the level set function values.
- `cut_cells`: A collection of indices representing the cut cells.

# Returns
An array of tuples representing the intersection points.
"""
function get_intersection_points(values::AbstractMatrix{T}, cut_cells) where {T}
    # Initialiser un tableau vide pour stocker les points d'intersection
    intersection_points = similar(cut_cells, NTuple{2,T}, 0)

    # Parcourir toutes les cellules coupées
    for index in cut_cells
        i, j = Tuple(index)

        # Parcourir toutes les arêtes de la cellule
        for (di, dj) in ((0, 1), (1, 0), (0, -1), (-1, 0))
            # Vérifier si la Level Set change de signe le long de cette arête
            if 1<=i+di<= size(values, 1) && 1 <= j+dj <= size(values, 2)
                if values[i, j] * values[i+di, j+dj] < 0
                    # Si c'est le cas, calculer le point d'intersection
                    t = values[i, j] / (values[i, j] - values[i+di, j+dj])
                    x_intersect = j + t * dj
                    y_intersect = i + t * di

                    # Ajouter le point d'intersection à la liste
                    push!(intersection_points, (x_intersect, y_intersect))
                end
            end
        end
    end

    intersection_points
end

"""
    get_intersection_points(values::AbstractArray{T,3}, cut_cells) where {T}

Compute the intersection points between the level set function and the cut cells.

# Arguments
- `values::AbstractArray{T,3}`: A 3D array representing the level set function values.
- `cut_cells`: An iterable containing the indices of the cut cells.

# Returns
An array of tuples representing the intersection points.
"""
function get_intersection_points(values::AbstractArray{T,3}, cut_cells) where {T}
    # Initialiser un tableau vide pour stocker les points d'intersection
    intersection_points = similar(cut_cells, NTuple{3,T}, 0)

    # Parcourir toutes les cellules coupées
    for index in cut_cells
        i, j, k = Tuple(index)

        # Parcourir toutes les arêtes de la cellule
        for (di, dj, dk) in ((0, 1, 0), (1, 0, 0), (0, -1, 0), (-1, 0, 0), (0, 0, 1), (0, 0, -1))
            # Vérifier si la Level Set change de signe le long de cette arête
            if 1<=i+di<= size(values, 1) && 1 <= j+dj <= size(values, 2) && 1 <= k+dk <= size(values, 3)
                if values[i, j, k] * values[i+di, j+dj, k+dk] < 0
                    # Si c'est le cas, calculer le point d'intersection
                    t = values[i, j, k] / (values[i, j, k] - values[i+di, j+dj, k+dk])
                    x_intersect = j + t * dj
                    y_intersect = i + t * di
                    z_intersect = k + t * dk

                    # Ajouter le point d'intersection à la liste
                    push!(intersection_points, (x_intersect, y_intersect, z_intersect))
                end
            end
        end
    end

    intersection_points
end

"""
    get_segment_midpoints(values::AbstractVector, cut_cells, intersection_points)

Compute the midpoints of the intersection points on each cut cell.

# Arguments
- `values::AbstractVector`: The values vector.
- `cut_cells`: The indices of the cut cells.
- `intersection_points`: The intersection points.

# Returns
An array of midpoints.

# Examples
"""
function get_segment_midpoints(values::AbstractVector, cut_cells, intersection_points)
    # Initialiser un tableau vide pour stocker les points médians
    midpoints = similar(intersection_points, 0)

    # Parcourir toutes les cellules coupées
    for index in cut_cells
        i, = Tuple(index)

        # Récupérer les points d'intersection sur cette cellule
        cell_points = [point for point in intersection_points if point[1] >= i && point[1] <= i+1]

        # Calculer le point médian de tous les points d'intersection sur cette cellule
        x_mid = sum(point[1] for point in cell_points) / length(cell_points)
        midpoint = (x_mid,)

        # Ajouter le point médian à la liste
        push!(midpoints, midpoint)
    end

    midpoints
end

"""
    get_segment_midpoints(values::AbstractMatrix, cut_cells, intersection_points)

Compute the midpoints of line segments defined by the intersection points on cut cells.

# Arguments
- `values::AbstractMatrix`: A matrix of values.
- `cut_cells`: An array of indices representing the cut cells.
- `intersection_points`: An array of intersection points.

# Returns
An array of midpoints of line segments.

# Example
"""
function get_segment_midpoints(values::AbstractMatrix, cut_cells, intersection_points)
    # Initialiser un tableau vide pour stocker les points médians
    midpoints = similar(intersection_points, 0)

    # Parcourir toutes les cellules coupées
    for index in cut_cells
        i, j = Tuple(index)

        # Récupérer les points d'intersection sur cette cellule
        cell_points = [point for point in intersection_points if point[1] >= j && point[1] <= j+1 && point[2] >= i && point[2] <= i+1]

        # Calculer le point médian de tous les points d'intersection sur cette cellule
        x_mid = sum(point[1] for point in cell_points) / length(cell_points)
        y_mid = sum(point[2] for point in cell_points) / length(cell_points)
        midpoint = (x_mid, y_mid)

        # Ajouter le point médian à la liste
        push!(midpoints, midpoint)
    end

    midpoints
end

"""
    get_segment_midpoints(values::AbstractArray{<:Any,3}, cut_cells, intersection_points)

This function calculates the midpoints of segments in a three-dimensional space.

# Arguments
- `values::AbstractArray{<:Any,3}`: An array of values.
- `cut_cells`: A collection of indices representing the cut cells.
- `intersection_points`: A collection of points representing the intersections.

# Returns
- `midpoints`: An array of midpoints.

# Example
"""
function get_segment_midpoints(values::AbstractArray{<:Any,3}, cut_cells, intersection_points)
    # Initialiser un tableau vide pour stocker les points médians
    midpoints = similar(intersection_points, 0)

    # Parcourir toutes les cellules coupées
    for index in cut_cells
        i, j, k = Tuple(index)

        # Récupérer les points d'intersection sur cette cellule
        cell_points = [point for point in intersection_points if point[1] >= j && point[1] <= j+1 && point[2] >= i && point[2] <= i+1 && point[3] >= k && point[3] <= k+1]

        # Calculer le point médian de tous les points d'intersection sur cette cellule
        x_mid = sum(point[1] for point in cell_points) / length(cell_points)
        y_mid = sum(point[2] for point in cell_points) / length(cell_points)
        z_mid = sum(point[3] for point in cell_points) / length(cell_points)
        midpoint = (x_mid, y_mid, z_mid)

        # Ajouter le point médian à la liste
        push!(midpoints, midpoint)
    end

    return midpoints
end


"""
    get_edges(mesh)

Given a mesh `mesh`, this function returns the edges of the mesh.

# Arguments
- `mesh`: A 1D, 2D, or 3D mesh represented as an array of coordinates.

# Returns
- `edges`: An array of edges of the mesh.

# Examples
"""
function get_edges(mesh)
    if length(mesh) == 1  # 1D
        x = mesh[1]
        edges = [[x[i], x[i+1]] for i in 1:length(x)-1]
    elseif length(mesh) == 2  # 2D
        x, y = mesh
        edges = []
        for i in 1:length(x)-1
            for j in 1:length(y)-1
                push!(edges, [[x[i], y[j]], [x[i+1], y[j]]])
                push!(edges, [[x[i], y[j]], [x[i], y[j+1]]])
            end
        end
    else  # 3D
        x, y, z = mesh
        edges = []
        for i in 1:length(x)-1
            for j in 1:length(y)-1
                for k in 1:length(z)-1
                    push!(edges, [[x[i], y[j], z[k]], [x[i+1], y[j], z[k]]])
                    push!(edges, [[x[i], y[j], z[k]], [x[i], y[j+1], z[k]]])
                    push!(edges, [[x[i], y[j], z[k]], [x[i], y[j], z[k+1]]])
                end
            end
        end
    end
    return edges
end

"""
    get_front_positions(front_sdf::SignedDistanceFunction, mesh::Tuple{AbstractVector})

Compute the positions of the front points where the signed distance function (SDF) is zero for each edge in the given mesh.

# Arguments
- `front_sdf::SignedDistanceFunction`: The signed distance function representing the front.
- `mesh::Tuple{AbstractVector}`: The mesh represented as a tuple of abstract vectors.

# Returns
An array of front positions where the SDF is zero.

# Example
"""
function get_front_positions(front_sdf::SignedDistanceFunction, mesh)
    # Obtenir toutes les arêtes du maillage
    aretes = get_edges(mesh)

    # Créer une liste pour stocker les points P où Phi=0
    points_P = []

    # Calculer le point P où Phi=0 pour chaque arête
    for arete in aretes
        # Paramétriser l'arête
        f(t) = front_sdf.sdf_function((arete[1] .+ t .* (arete[2] .- arete[1]))...)
        try
            # Utiliser la méthode de la bissection pour trouver le point P
            t0 = find_zero(f, (0.0, 1.0), Bisection())
            x0 = arete[1] + t0 * (arete[2] - arete[1])
            
            # Ajouter le point P à la liste des points P
            push!(points_P, x0)
        catch e
            #println("Aucun point P trouvé sur l'arête entre les points $(arete[1]) et $(arete[2])")
        end
    end

    return points_P
end
