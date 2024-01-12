#using Ferrite

using LinearAlgebra, SparseArrays

include("preprocess.jl")

diffusivity_constant = 1

density = 0.5

time_step = 1

mass_transfer_rate = 1

#application of boundary conditions, preset faces that are already calculated include top_face_nodes, bottom_face_nodes,
#right_face_nodes, left_face_nodes

function four_node_parametric_shape_function(xi, eta)

    N1 = 0.25*(1-xi)*(1-eta)

    N2 = 0.25*(1+xi)*(1-eta)

    N3 = 0.25*(1+xi)*(1+eta)

    N4 = 0.25*(1-xi)*(1+eta)

    return [N1 N2 N3 N4]

end

function gradient_four_node_parametric_shape_function(xi, eta)

    N1_partial_xi  = -0.25*(1-eta)

    N1_partial_eta = -0.25*(1-xi)

    N2_partial_xi  = 0.25*(1-eta)

    N2_partial_eta = -0.25*(1+xi)

    N3_partial_xi  = 0.25*(1+eta)

    N3_partial_eta = 0.25*(1+xi)

    N4_partial_xi  = -0.25*(1+eta)

    N4_partial_eta = 0.25*(1-xi)

    return [N1_partial_xi N2_partial_xi N3_partial_xi N4_partial_xi;
    N1_partial_eta N2_partial_eta N3_partial_eta N4_partial_eta]
    
     
end

#creation of square stiffness matrix with size : (total dof x total dof)

stiffness_matrix = spzeros(total_degrees_of_freedom, total_degrees_of_freedom)

mass_matrix = spzeros(total_degrees_of_freedom, total_degrees_of_freedom)

gauss_points_for_quadrilateral_element = [-1 -1; 1 -1; 1 1; -1 1] ./ sqrt(3)

weight = 1

#print(gauss_points_for_quadrilateral_element)

#print(size(gauss_points_for_quadrilateral_element))

for element in element_connectivity

    #print(size(element_stiffness))

    #note to self: when going to 3-D, nodal_coordinates[element,1:2] ---> nodal_coordinates[element,:] to index all (x,y,z)
    #instead of just (x,y) for 2-D case

    element_node_coordinates = nodal_coordinates[element,2:3]

    #print(element_node_coordinates)

    #print(size(element_node_coordinates))

    for gauss_point in eachrow(gauss_points_for_quadrilateral_element)

        #println(gauss_point)

        shape_function_N = four_node_parametric_shape_function(gauss_point[1], gauss_point[2])

        derivative_of_N = gradient_four_node_parametric_shape_function(gauss_point[1], gauss_point[2])

        #derivative_of_N is 2x4 matrix 

        jacobian_of_N = transpose(derivative_of_N * element_node_coordinates)

        # B_matrix = ∂N/∂x∂y = [∂N1/∂x ∂N2/∂x ∂N3/∂x ∂N4/∂x]
        #                      [∂N1/∂y ∂N2/∂y ∂N3/∂y ∂N4/∂y]

        #B_matrix = (1/2) .* (jacobian_of_N * derivative_of_N) <--- seems to be incorrect D:

        B_matrix = inv(jacobian_of_N) * derivative_of_N

        #check dimensions of jacobian_of_N (is supposed to be 2x2)

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

println("The size of the stiffness matrix is: ", size(stiffness_matrix))

println("The size of the mass matrix is:      ", size(mass_matrix))
    
