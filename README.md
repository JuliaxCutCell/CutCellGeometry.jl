# CutCellGeometry

[![Build Status](https://github.com/fastaxx/CutCellGeometry.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fastaxx/CutCellGeometry.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Overview
CutCellGeometry.jl is a Julia package designed for geometry immersion into meshes. It offers functionalities for constructing geometries in 1D, 2D, and 3D spaces, as well as handling geometry motion.

## Usage
```
using CutCellMesh # Mesh generation
using CutCellGeometry # Geometry definition

# CartesianGrid 
grid = CartesianGrid((10, 10), (0.1, 0.1))
mesh = generate_mesh(grid)
domain = ((minimum(mesh[1]), maximum(mesh[1])), (minimum(mesh[2]), maximum(mesh[2])))

# SignedDistanceFunction
sphere_function = (x, y, z) -> sqrt(x^2 + y^2 + z^2) - 1.0
cube_function = (x, y, z) -> max(abs(x), max(abs(y), abs(z))) - 1.0
identity_transform = (x, y, z, t) -> (x, y, z)

# Create Geometry
sphere_sdf = SignedDistanceFunction(sphere_function, identity_transform, domain, false)
cube_sdf = SignedDistanceFunction(cube_function, identity_transform, domain, false)

# Union of two signed distance functions
union_sdf = sphere_sdf ⊔ cube_sdf

# Transformation function
move_transform = (x, y, z, t) -> (x + t, y, z)

# Move the union of two signed distance functions
moving_sdf = SignedDistanceFunction(union_sdf.sdf_function, move_transform, union_sdf.domain, true)

t=2.0
plot_sdf(moving_sdf, domain, t)
readline()
```

## ToDo
- Normals, Curvatures
- Parametric Curve
- Multiple Phase Initialization : N=3 =>Phi=0/0.5/1
- Driver for VOFI/CartesianGeometry(vlc) or reimplement (mouais non)
- Docs + Notebooks