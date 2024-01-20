using LinearAlgebra, SparseArrays

println("This is preprocessing file for Mesoscale Concrete Diffusion Model")

include("preprocess.jl")

println("This is initilization of Mesoscale Concrete Diffusion Model")

include("initilization.jl")

function boundary_condition(concentration, affected_nodes)

    fixed_concentration = concentration .* ones(size(affected_nodes)[1], 1)

    affected_nodes = affected_nodes[: , 1]

    return [affected_nodes fixed_concentration]

end

#boundary_condition function output looks like the following:
#    node ids affected by fixed concentration   |   fixed concentration value

#inputs are the boundary condition and the unpartitioned stiffness_matrix

function partitioned_boundary_condition_function(stiffness_matrix, boundary_condition)

    removed_rows_and_columns = convert(Array{Int,1}, boundary_condition[:, 1])

    #steps to partition boundary condition -> remove rows affected by boundary condition in the stiffness matrix, 
    #multiply the columns affected by the boundary condition with the known value and subtract it from the rhs of the sys. of linear equation

    stiffness_matrix_removed_rows = stiffness_matrix[setdiff(1:end, removed_rows_and_columns), :]

    applied_boundary_condition = stiffness_matrix_removed_rows[: , removed_rows_and_columns] * boundary_condition[: , 2]

    return applied_boundary_condition

end

function partitioned_stiffness_matrix_function(stiffness_matrix, boundary_condition)

    #converting applied b.c. nodes into an array of integers to index / slice out...

    removed_rows_and_columns = convert(Array{Int,1}, boundary_condition[: , 1])

    partitioned_matrix = stiffness_matrix[setdiff(1:end, removed_rows_and_columns), setdiff(1:end, removed_rows_and_columns)]

    return partitioned_matrix

end

function partitioned_mass_matrix_function(mass_matrix, boundary_condition)

    #same as stiffness matrix partitioning function

    removed_rows_and_columns = convert(Array{Int,1}, boundary_condition[:, 1])

    partitioned_mass_matrix = mass_matrix[setdiff(1:end, removed_rows_and_columns), setdiff(1:end, removed_rows_and_columns)]

    return partitioned_mass_matrix

end


#lots of repetitive inputs into the functions below --> possibly write it all as one function later?

fixed_concentration_top = boundary_condition(1, top_face_nodes)

# partitioned_stiffness_matrix = partitioned_stiffness_matrix_function(stiffness_matrix, fixed_concentration_top)

# partitioned_mass_matrix = partitioned_mass_matrix_function(mass_matrix, fixed_concentration_top)

partitioned_boundary_condition = partitioned_boundary_condition_function(stiffness_matrix, fixed_concentration_top)

# println("the partitioned stiffness matrix size is: ", size(partitioned_stiffness_matrix))

# println("the partitioned mass matrix size is: ", size(partitioned_mass_matrix))

# println("the partitioned essential boundary condition size is: ", size(partitioned_boundary_condition))

#initialization of stuff

concentration_vector = spzeros(size(partitioned_boundary_condition)[1], 1)

#concentration_vector = zeros(total_degrees_of_freedom, 1)

flux_vector = spzeros(size(partitioned_boundary_condition)[1], 1)

#flux_vector =  zeros(total_degrees_of_freedom, 1)

#total time in seconds
time_end = 40

global time_step = 1

global time = 0

while time_end > time

    println("the time is: ", time)

    global old_concentration_vector = concentration_vector

    local stiffness_matrix, mass_matrix = stiffness_and_mass_matrix_2D("quad4", element_connectivity, nodal_coordinates, 0.001, 0.01, 0.99, time_step)

    #is creating an external flux vector necessary for this problem?

    local stiffness_matrix = stiffness_matrix + mass_matrix

    local partitioned_stiffness_matrix = partitioned_stiffness_matrix_function(stiffness_matrix, fixed_concentration_top)

    local partitioned_mass_matrix = partitioned_mass_matrix_function(mass_matrix, fixed_concentration_top)

    local partitioned_boundary_condition = partitioned_boundary_condition_function(stiffness_matrix, fixed_concentration_top)

    # assumption - there is no external flux vector --> set initial flux to zero?

    global flux_vector = spzeros(size(partitioned_boundary_condition)[1], 1)

    global flux_vector = flux_vector + partitioned_mass_matrix * old_concentration_vector - partitioned_boundary_condition

    global concentration_vector = partitioned_stiffness_matrix \ flux_vector

    local residual = flux_vector - partitioned_stiffness_matrix * concentration_vector

    global time = time + time_step

    println("concentration is:")

    println(concentration_vector)

    println("flux is:")

    println(flux_vector)

    println("the residual is: ", sum(residual))

end