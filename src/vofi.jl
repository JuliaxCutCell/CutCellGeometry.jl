const T = Float64

"""
    calculate_first_order_moments(levelset, xyz)

Calculate the first-order moments of a levelset function.

# Arguments
- `levelset`: The levelset function.
- `xyz`: The coordinates of the points where the moments are calculated.

# Returns
- `V`: The zeroth-order moment.
- `bary`: The first-order moments.
- `As`: The second-order moments.

"""
function calculate_first_order_moments(levelset, xyz)
    # first-kind moments
    V, bary = integrate(Tuple{0}, levelset, xyz, T, zero)
    As = integrate(Tuple{1}, levelset, xyz, T, zero)
    
    return V, bary, As
end

"""
    first_aperture_1D(V, bary, As)

Compute the first aperture in 1D.

# Arguments
- `V`: Vector of volumes.
- `bary`: Barycentric coordinates.
- `As`: Vector of areas.

# Returns
- `v_diag`: Diagonal matrix of volumes.
- `ax_diag`: Diagonal matrix of areas.

# Examples
"""
function first_aperture_1D(V, bary, As)
    v_diag = spdiagm(0 => V) # Diagonal matrix of Volumes : V
    ax_diag = spdiagm(0 => As[1]) # Diagonal matrix of Ax

    return v_diag, ax_diag
end

"""
    first_aperture_2D(V, bary, As)

Compute the diagonal matrices `v_diag`, `ax_diag`, and `ay_diag` based on the input parameters.

# Arguments
- `V`: Vector of volumes.
- `bary`: Barycentric coordinates.
- `As`: Vector of areas.

# Returns
- `v_diag`: Diagonal matrix of volumes.
- `ax_diag`: Diagonal matrix of Ax.
- `ay_diag`: Diagonal matrix of Ay.
"""
function first_aperture_2D(V, bary, As)
    v_diag = spdiagm(0 => V) # Diagonal matrix of Volumes : V
    ax_diag = spdiagm(0 => As[1]) # Diagonal matrix of Ax
    ay_diag = spdiagm(0 => As[2]) # Diagonal matrix of Ay

    return v_diag, ax_diag, ay_diag
end

"""
    first_aperture_3D(V, bary, As)

Compute the first aperture in 3D.

# Arguments
- `V`: Vector of volumes.
- `bary`: Barycentric coordinates.
- `As`: Vector of areas.

# Returns
- `v_diag`: Diagonal matrix of volumes.
- `ax_diag`: Diagonal matrix of Ax.
- `ay_diag`: Diagonal matrix of Ay.
- `az_diag`: Diagonal matrix of Az.
"""
function first_aperture_3D(V, bary, As)
    v_diag = spdiagm(0 => V) # Diagonal matrix of Volumes : V
    ax_diag = spdiagm(0 => As[1]) # Diagonal matrix of Ax
    ay_diag = spdiagm(0 => As[2]) # Diagonal matrix of Ay
    az_diag = spdiagm(0 => As[3]) # Diagonal matrix of Az

    return v_diag, ax_diag, ay_diag, az_diag
end


"""
    calculate_second_order_moments(levelset, xyz, bary)

Calculate the second-order moments of a level set function.

# Arguments
- `levelset`: The level set function.
- `xyz`: The coordinates of the points where the moments are calculated.
- `bary`: The barycentric coordinates of the points where the moments are calculated.

# Returns
- `Ws`: The second-order moments of the level set function.
- `Bs`: The second-order moments of the gradient of the level set function.
"""
function calculate_second_order_moments(levelset, xyz, bary)
    # Moments (2nd order)
    Ws = integrate(Tuple{0}, levelset, xyz, T, zero, bary)
    Bs = integrate(Tuple{1}, levelset, xyz, T, zero, bary)

    return Ws, Bs
end

"""
    second_aperture_1D(Ws, Bs)

Compute the diagonal matrices `wx_diag` and `bx_diag` based on the input surfaces `Ws` and `Bs`.

# Arguments
- `Ws`: Array of surface values for Wx.
- `Bs`: Array of surface values for Bx.

# Returns
- `wx_diag`: Diagonal matrix of Wx.
- `bx_diag`: Diagonal matrix of Bx.
"""
function second_aperture_1D(Ws, Bs)
    Wx = Ws[1] # Surface in x : Wx
    Bx = Bs[1] # Surface in x : Bx

    wx_diag = spdiagm(0 => Wx) # Diagonal matrix of Wx
    bx_diag = spdiagm(0 => Bx) # Diagonal matrix of Bx

    return wx_diag, bx_diag
end

"""
    second_aperture_2D(Ws, Bs)

Constructs diagonal matrices from the given surface values.

# Arguments
- `Ws::Vector{Float64}`: Surface values in x and y directions for the W surface.
- `Bs::Vector{Float64}`: Surface values in x and y directions for the B surface.

# Returns
- `wx_diag::SparseMatrixCSC{Float64, Int64}`: Diagonal matrix of Wx.
- `wy_diag::SparseMatrixCSC{Float64, Int64}`: Diagonal matrix of Wy.
- `bx_diag::SparseMatrixCSC{Float64, Int64}`: Diagonal matrix of Bx.
- `by_diag::SparseMatrixCSC{Float64, Int64}`: Diagonal matrix of By.
"""
function second_aperture_2D(Ws, Bs)
    Wx = Ws[1] # Surface in x : Wx
    Wy = Ws[2] # Surface in y : Wy
    Bx = Bs[1] # Surface in x : Bx
    By = Bs[2] # Surface in y : By

    wx_diag = spdiagm(0 => Wx) # Diagonal matrix of Wx
    wy_diag = spdiagm(0 => Wy) # Diagonal matrix of Wy
    bx_diag = spdiagm(0 => Bx) # Diagonal matrix of Bx
    by_diag = spdiagm(0 => By) # Diagonal matrix of By

    return wx_diag, wy_diag, bx_diag, by_diag
end

"""
    second_aperture_3D(Ws, Bs)

Constructs diagonal matrices from the input arrays `Ws` and `Bs` and returns them.

# Arguments
- `Ws::Vector{Float64}`: Array of surface values in x, y, and z directions.
- `Bs::Vector{Float64}`: Array of surface values in x, y, and z directions.

# Returns
- `wx_diag::SparseMatrixCSC{Float64, Int64}`: Diagonal matrix constructed from `Wx`.
- `wy_diag::SparseMatrixCSC{Float64, Int64}`: Diagonal matrix constructed from `Wy`.
- `wz_diag::SparseMatrixCSC{Float64, Int64}`: Diagonal matrix constructed from `Wz`.
- `bx_diag::SparseMatrixCSC{Float64, Int64}`: Diagonal matrix constructed from `Bx`.
- `by_diag::SparseMatrixCSC{Float64, Int64}`: Diagonal matrix constructed from `By`.
- `bz_diag::SparseMatrixCSC{Float64, Int64}`: Diagonal matrix constructed from `Bz`.
"""
function second_aperture_3D(Ws, Bs)
    Wx = Ws[1] # Surface in x : Wx
    Wy = Ws[2] # Surface in y : Wy
    Wz = Ws[3] # Surface in z : Wz
    Bx = Bs[1] # Surface in x : Bx
    By = Bs[2] # Surface in y : By
    Bz = Bs[3] # Surface in z : Bz

    wx_diag = spdiagm(0 => Wx) # Diagonal matrix of Wx
    wy_diag = spdiagm(0 => Wy) # Diagonal matrix of Wy
    wz_diag = spdiagm(0 => Wz) # Diagonal matrix of Wz
    bx_diag = spdiagm(0 => Bx) # Diagonal matrix of Bx
    by_diag = spdiagm(0 => By) # Diagonal matrix of By
    bz_diag = spdiagm(0 => Bz) # Diagonal matrix of Bz

    return wx_diag, wy_diag, wz_diag, bx_diag, by_diag, bz_diag
end