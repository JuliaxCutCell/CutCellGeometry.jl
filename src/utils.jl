"""
    compute_velocity(transform, x, t)

Compute the velocity at a given position and time using a transformation function.

# Arguments
- `transform`: A function that maps a position `(x, y)` and time `t` to a velocity vector.
- `x`: The position at which to compute the velocity.
- `t`: The time at which to compute the velocity.

# Returns
The velocity vector at the given position and time.
"""
function compute_velocity(transform, x, t)
    # Calculer le jacobien J
    J = ForwardDiff.jacobian(u -> [transform(u[1], u[2], t)...], [x[1], x[2]])
    @show J
    # Calculer la dérivée dot
    dot = ForwardDiff.derivative(u -> sum([transform(x[1], x[2], u)...]), t)
    @show dot
    # Calculer la vitesse
    velocity = -dot ./ J
    @show velocity
    return velocity
end

"""
    gradient_Phi(sdf, point)

Compute the gradient of a signed distance function (SDF) at a given point.

# Arguments
- `sdf`: A function representing the signed distance function.
- `point`: A tuple or array representing the coordinates of the point.

# Returns
The gradient of the SDF at the given point.
"""
function gradient_Phi(sdf, point)
    f = x -> sdf(x...)
    return ForwardDiff.gradient(f, [point...])
end

"""
    normal(sdf, point)

Compute the normalized normal vector at a given point on the surface defined by the signed distance function (sdf).

# Arguments
- `sdf`: A signed distance function representing the surface.
- `point`: The point at which to compute the normal vector.

# Returns
The normalized normal vector at the given point.
"""
function normal(sdf, point)
    return gradient_Phi(sdf, point) ./ norm(gradient_Phi(sdf, point))
end

"""
    tangent_vector(normal)

Generate a tangent vector orthogonal to the given normal vector.

# Arguments
- `normal`: The normal vector.

# Returns
The tangent vector orthogonal to the normal vector.
"""
function tangent_vector(normal)
    # Générer un vecteur aléatoire de la même longueur que le vecteur normal
    random_vector = rand(length(normal))
    
    # Calculer le vecteur tangentiel
    if length(normal) == 3
        tangent = cross(normal, random_vector)
    elseif length(normal) == 2
        tangent = [-normal[2], normal[1]]
    else
        tangent = normal
    end
    
    return tangent
end

"""
    hessian_Phi(sdf, point)

Compute the Hessian matrix of a signed distance function (sdf) at a given point.

# Arguments
- `sdf::Function`: The signed distance function.
- `point::Vector`: The point at which to compute the Hessian matrix.

# Returns
- `hessian::Matrix`: The Hessian matrix of the signed distance function at the given point.
"""
function hessian_Phi(sdf, point)
    f = x -> sdf(x...)
    return ForwardDiff.hessian(f, [point...])
end

"""
    curvature(sdf, point)

Compute the Gaussian curvature and mean curvature at a given point on the surface defined by the signed distance function (sdf).

# Arguments
- `sdf::Function`: The signed distance function that defines the surface.
- `point::Vector{Float64}`: The coordinates of the point on the surface.

# Returns
- `gaussian_curvature::Float64`: The Gaussian curvature at the given point.
- `mean_curvature::Float64`: The mean curvature at the given point.

# References
- [Goldman, R., & Kimmel, R. (2005). Curvature formulas for implicit curves and surfaces. Graphical Models, 67(4), 369-387.](https://u.cs.biu.ac.il/~katzmik/goldman05.pdf)
"""
function curvature(sdf, point)
    # Calculer le gradient de Phi au point
    gradient = gradient_Phi(sdf, point)
    
    # Calculer la hessienne de Phi au point
    hessian = hessian_Phi(sdf, point)
    adjoint_hessian = det(hessian) * inv(hessian)
    
    # Calculer la courbure gaussienne
    gaussian_curvature = gradient' * adjoint_hessian * gradient / norm(gradient)^4

    # Calculer la courbure moyenne
    mean_curvature = (gradient' * hessian * gradient - norm(gradient)^2 * tr(hessian))/(2*norm(gradient)^3)
    
    return gaussian_curvature, mean_curvature
end
