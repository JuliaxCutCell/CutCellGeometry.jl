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
circle_function = (x, y) -> sqrt(x^2 + y^2) - 0.5
square_function = (x, y) -> max(abs(x), abs(y)) - 0.5
identity_transform = (x, y, t) -> (x, y)

# Create Geometry
circle_sdf = SignedDistanceFunction(circle_function, identity_transform, domain, false)
square_sdf = SignedDistanceFunction(square_function, identity_transform, domain, false)

# Union of two signed distance functions
union_sdf = circle_sdf ⊔ square_sdf

# Transformation function
move_transform = (x, y, t) -> ((x-0.25) + 2*t, (y-0.25))

# Move the union of two signed distance functions
moving_sdf = SignedDistanceFunction(union_sdf.sdf_function, move_transform, union_sdf.domain, true)

# Print the normal vector at a given point
println("Normal au point (0.5, 0.5) pour le cercle", normal(circle_sdf.sdf_function, (0.5, 0.5)))

# Calculer le vecteur tangentiel
tangent = tangent_vector(normal(circle_sdf.sdf_function, (0.5, 0.5)))
println("Vecteur tangentiel au point (0.5, 0.5) pour le cercle", tangent)

# Calculer la courbure pour le cercle
println("Courbure pour le cercle au point (0.5, 0.5)", curvature(circle_function.sdf_function, (0.5, 0.5)))

# Plot
t=0.0
plot_sdf(moving_sdf, domain, t)
readline()
```

## ToDo
- Marching Cubes algorithm to generate a mesh from the Signed Distance Function.
- Add Reference Function (Sphere, Torus, Cube ,...), Infinite 3D Primitives, Text, Images extruded
- Positioning Function (translate, scale, ...), Repetition
- Parametric Curve
- Multiple Phase Initialization : N=3 =>Phi=0/0.5/1
- Limits bounds domain when transform(t)
- Driver for VOFI/CartesianGeometry(vlc) or reimplement (mouais non)
- Docs + Notebooks