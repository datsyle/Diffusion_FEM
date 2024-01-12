using Plots

function read_mesh(file_name::AbstractString, element_type)
    x_coords = Float64[]
    y_coords = Float64[]
    z_coords = Float64[]
    elems = Vector{Vector{Int}}()

    if element_type == 1
        open(file_name) do f
            for line in eachline(f)
                if startswith(line, '*')
                    continue
                end

                tokens = split(strip(line), ',')
                if length(tokens) == 5
                    elem_id, node_ids... = tokens
                    push!(elems, [parse(Int, n) for n in node_ids])
                elseif length(tokens) == 4
                    node_id, x, y, z = tokens
                    push!(x_coords, parse(Float64, x))
                    push!(y_coords, parse(Float64, y))
                    push!(z_coords, parse(Float64, z))
                end
            end
        end
    elseif element_type == 2
        open(file_name) do f
            is_elem_connectivity = 0
            for line in eachline(f)
                if startswith(line, '*')
                    continue
                end

                tokens = split(strip(line), ',')
                if is_elem_connectivity == 0
                    if length(tokens) == 4
                        node_id, x, y, z = tokens
                        push!(x_coords, parse(Float64, x))
                        push!(y_coords, parse(Float64, y))
                        push!(z_coords, parse(Float64, z))
                        if node_id == 1
                            is_elem_connectivity += 1
                        end
                    end
                end
            end
        end
    end

    return ( x_coords, y_coords, z_coords, elems)
end

#using read_mesh function to create 4 variables that contain all x-coordinates in mesh, all y-coordinates in mesh,
#all z-coordinates in mesh, and global nodal coordinates

#if you want to change the mesh file used -----> adjust mesh_information variable
#element type must be defined for mesh_information -----> 1 is for square elements

#jl_preprocessing_ex.inp

mesh_information = read_mesh("julia_sample_mesh.inp", 1)

node_ids = 1:length(mesh_information[1])

X_coordinates = mesh_information[1]

Y_coordinates = mesh_information[2]

Z_coordinates = mesh_information[3]

nodal_coordinates = [node_ids X_coordinates Y_coordinates Z_coordinates]

#print(reshape(nodal_coordinates[1])

#identifies nodes at the top boundary of the mesh <----- are loads applied at the nodes or at the faces of the elements in FEA?
#this variable only works on cubes, squares, and rectangles not a generalized form for all shapes

tolerance = 10^-7

top_face_nodes = nodal_coordinates[abs.(Y_coordinates .- maximum(Y_coordinates)) .< tolerance,:]

bottom_face_nodes = nodal_coordinates[abs.(Y_coordinates .- minimum(Y_coordinates)) .< tolerance,:]

right_face_nodes = nodal_coordinates[abs.(X_coordinates .- maximum(X_coordinates)) .< tolerance,:]

left_face_nodes = nodal_coordinates[abs.(X_coordinates .- minimum(X_coordinates)) .< tolerance,:]

#calling element connectivity array to a variable

element_connectivity = mesh_information[4]

#print("size of element connectivity is :", size(element_connectivity))

#total degrees of freedom of problem are 1 * number of nodes, with the degree being concentration <--- confirm later with professor berkin
#add in check that there are the same number of x_coords, y_coords, and z_coords

total_degrees_of_freedom = size(nodal_coordinates)[1]

total_number_of_elements = size(element_connectivity)[1]

#creating and saving the plot to check the geometry of the mesh

#mesh_plot = scatter!(X_coordinates, Y_coordinates)

#display(mesh_plot)

#savefig(mesh_plot, "mesh_plot_check.pdf")

println("The total degrees of freedom is: ",total_degrees_of_freedom)

println("The total number of elements is: ", total_number_of_elements)


