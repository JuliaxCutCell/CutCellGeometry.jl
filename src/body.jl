"""
    struct SignedDistanceFunction{T}

A struct representing a signed distance function.

# Fields
- `sdf_function::Function`: The function that computes the signed distance function.
- `transform::Function`: The coordinate transformation function.
- `domain::T`: The domain of definition for the function.
- `is_moving::Bool`: Indicates whether the geometry is in motion.

"""
struct SignedDistanceFunction{T}
    sdf_function::Function
    transform::Function
    domain::T
    is_moving::Bool
end


"""
    evaluate_sdf(sdf::SignedDistanceFunction, t, x...)

Evaluate the signed distance function (SDF) at the given coordinates.

Parameters
----------
sdf : SignedDistanceFunction
    The signed distance function to evaluate.
t : Float64
    The time parameter.
x : Float64...
    The coordinates at which to evaluate the SDF.

Returns
-------
Float64
    The value of the SDF at the given coordinates.

"""
function evaluate_sdf(sdf::SignedDistanceFunction, t, x...)
    if sdf.is_moving
        transformed_coordinates = sdf.transform(x..., t)
    else
        transformed_coordinates = x
    end
    return sdf.sdf_function(transformed_coordinates...)
end


"""
    ⊔(a::SignedDistanceFunction{T}, b::SignedDistanceFunction{T}) where T

Compute the minimum signed distance function between two signed distance functions `a` and `b`.

# Arguments
- `a::SignedDistanceFunction{T}`: The first signed distance function.
- `b::SignedDistanceFunction{T}`: The second signed distance function.

# Returns
A new `SignedDistanceFunction` object representing the minimum signed distance function between `a` and `b`.

"""
function ⊔(a::SignedDistanceFunction{T}, b::SignedDistanceFunction{T}) where T
    sdf(x::Vararg{Float64}) = min(a.sdf_function(x...), b.sdf_function(x...))
    SignedDistanceFunction(sdf, a.transform, a.domain, a.is_moving || b.is_moving)
end

"""
    ⊓(a::SignedDistanceFunction{T}, b::SignedDistanceFunction{T}) where T

Compute the union of two signed distance functions.

This function takes two `SignedDistanceFunction` objects `a` and `b` and returns a new `SignedDistanceFunction` object that represents the union of `a` and `b`. The union is computed by taking the maximum signed distance value at each point.

# Arguments
- `a::SignedDistanceFunction{T}`: The first signed distance function.
- `b::SignedDistanceFunction{T}`: The second signed distance function.

# Returns
- `SignedDistanceFunction{T}`: A new signed distance function representing the union of `a` and `b`.

"""
function ⊓(a::SignedDistanceFunction{T}, b::SignedDistanceFunction{T}) where T
    sdf(x::Vararg{Float64}) = max(a.sdf_function(x...), b.sdf_function(x...))
    SignedDistanceFunction(sdf, a.transform, a.domain, a.is_moving || b.is_moving)
end

"""
    ⊖(a::SignedDistanceFunction{T}, b::SignedDistanceFunction{T}) where T

Compute the signed distance function for the difference of two signed distance functions `a` and `b`.

# Arguments
- `a::SignedDistanceFunction{T}`: The first signed distance function.
- `b::SignedDistanceFunction{T}`: The second signed distance function.

# Returns
A new `SignedDistanceFunction` representing the signed distance function for the difference of `a` and `b`.

"""
function ⊖(a::SignedDistanceFunction{T}, b::SignedDistanceFunction{T}) where T
    sdf(x::Vararg{Float64}) = max(a.sdf_function(x...), -b.sdf_function(x...))
    SignedDistanceFunction(sdf, a.transform, a.domain, a.is_moving || b.is_moving)
end

"""
    complement(a::SignedDistanceFunction{T})

Create a new signed distance function that represents the complement of the input signed distance function `a`.

# Arguments
- `a::SignedDistanceFunction{T}`: The input signed distance function.

# Returns
A new `SignedDistanceFunction` object representing the complement of `a`.

"""
function complement(a::SignedDistanceFunction{T}) where T
    sdf(x::Vararg{Float64}) = -a.sdf_function(x...)
    SignedDistanceFunction(sdf, a.transform, a.domain, a.is_moving)
end