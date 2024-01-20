#using Ferrite

using LinearAlgebra, SparseArrays

include("preprocess.jl")

#creation of square stiffness matrix with size : (total dof x total dof)

stiffness_matrix = spzeros(total_degrees_of_freedom, total_degrees_of_freedom)

mass_matrix = spzeros(total_degrees_of_freedom, total_degrees_of_freedom)

#create function for stiffness & mass matrix:

function stiffness_and_mass_matrix_2D(element_type, element_connectivity, nodal_coordinates, diffusivity_constant, mass_transfer_rate, density, time_step)

    if element_type == "quad4"

        include("quad4_values.jl")
    
    end

    # stiffness_matrix = spzeros(total_degrees_of_freedom, total_degrees_of_freedom)

    # mass_matrix = spzeros(total_degrees_of_freedom, total_degrees_of_freedom)

    for element in element_connectivity

        #note to self: when going to 3-D, nodal_coordinates[element,1:2] ---> nodal_coordinates[element,:] to index all (x,y,z)
        #instead of just (x,y) for 2-D case
    
        element_node_coordinates = nodal_coordinates[element,2:3]
        
        for gauss_point in eachrow(gauss_points)

            #shape_function_N is supposed to be 1x4 matrix

            shape_function_N = Base.invokelatest(shape_function, gauss_point[1], gauss_point[2])

            #derivative_of_N is supposed to be 2x4 matrix
    
            derivative_of_N = Base.invokelatest(gradient_shape_function, gauss_point[1], gauss_point[2])

            #jacobian_of_N is supposed to be 2x2

            jacobian_of_N = transpose(derivative_of_N * element_node_coordinates)
    
            # B_matrix = ∂N/∂x∂y = [∂N1/∂x ∂N2/∂x ∂N3/∂x ∂N4/∂x]
            #                      [∂N1/∂y ∂N2/∂y ∂N3/∂y ∂N4/∂y]
        
            B_matrix = inv(jacobian_of_N) * derivative_of_N
        
            determinant_of_jacobian_N = det(jacobian_of_N)

            #dv is the volume element
        
            dv = determinant_of_jacobian_N * weight
            
            global element_stiffness = diffusivity_constant * (transpose(B_matrix) * B_matrix) * dv
            
            global mass_matrix_element = (1 / time_step) * density * mass_transfer_rate * (transpose(shape_function_N) * shape_function_N) * dv
    
            #println("sizes of the following matricies: ")
    
            #println("derivative of N:", size(derivative_of_N))
    
            #println("B_matrix:", size(B_matrix))
    
            #println("jaobian of N:", size(jacobian_of_N))
    
            #println("element stiffness:", size(element_stiffness))
    
            #println("mass matrix element:", size(mass_matrix_element))
    
            #println("dv:", size(dv))
    
        end
    
        global stiffness_matrix[element, element] = stiffness_matrix[element, element] + element_stiffness
    
        global mass_matrix[element, element] = mass_matrix[element, element] + mass_matrix_element
    
    end
    
    return [stiffness_matrix, mass_matrix]

end


#stiffness_matrix = stiffness_and_mass_matrix_2D("quad4", element_connectivity, nodal_coordinates, 1, 1, 0.5, 1)[1]

#mass_matrix = stiffness_and_mass_matrix_2D("quad4", element_connectivity, nodal_coordinates, 1, 1, 0.5, 1)[2]

#println(stiffness_matrix)

#println(mass_matrix)

#println("The size of the stiffness matrix is: ", size(stiffness_matrix))

#println("The size of the mass matrix is:      ", size(mass_matrix))





