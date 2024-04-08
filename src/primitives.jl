# 1D 
front_1D(position::Float64) = (x, _=0) -> x - position

# Hypersphere
cercle(center::Tuple{Float64, Float64}, radius::Float64) = (x, y, _=0) -> sqrt((x - center[1])^2 + (y - center[2])^2) - radius
sphere(center::Tuple{Float64, Float64, Float64}, radius::Float64) = (x, y, z) -> sqrt((x - center[1])^2 + (y - center[2])^2 + (z - center[3])^2) - radius

# Hypercuboid
rectangle(lower_left::Tuple{Float64, Float64}, upper_right::Tuple{Float64, Float64}) = (x, y, _=0) -> max(max(lower_left[1] - x, x - upper_right[1]), max(lower_left[2] - y, y - upper_right[2]))
cuboid(lower_left::Tuple{Float64, Float64, Float64}, upper_right::Tuple{Float64, Float64, Float64}) = (x, y, z) -> max(max(lower_left[1] - x, x - upper_right[1]), max(max(lower_left[2] - y, y - upper_right[2]), max(lower_left[3] - z, z - upper_right[3])))

# HyperToroid
torus(center::Tuple{Float64, Float64, Float64}, major_radius::Float64, minor_radius::Float64) = (x, y, z) -> (sqrt((x - center[1])^2 + (y - center[2])^2) - major_radius)^2 + (z - center[3])^2 - minor_radius^2

# Hyperplan
plane(normal::Tuple{Float64, Float64}, d::Float64) = (x, y, _=0) -> dot(normal, [x, y]) - d
plane(normal::Tuple{Float64, Float64, Float64}, d::Float64) = (x, y, z) -> dot(normal, [x, y, z]) - d

# HyperCylinder
cylinder(center::Tuple{Float64, Float64}, radius::Float64) = (x, y, _=0) -> sqrt((x - center[1])^2 + (y - center[2])^2) - radius
cylinder(center::Tuple{Float64, Float64, Float64}, radius::Float64) = (x, y, z) -> sqrt((x - center[1])^2 + (y - center[2])^2) - radius