function shape_function(xi, eta)

    N1 = 0.25*(1-xi)*(1-eta)

    N2 = 0.25*(1+xi)*(1-eta)

    N3 = 0.25*(1+xi)*(1+eta)

    N4 = 0.25*(1-xi)*(1+eta)

    return [N1 N2 N3 N4]

end

function gradient_shape_function(xi, eta)

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

gauss_points = [-1 -1; 1 -1; 1 1; -1 1] ./ sqrt(3)

weight = 1
