include("solver.jl")

using WriteVTK

reshape(X_coordinates , 1 , :)

reshape(Y_coordinates , 1 , :)

reshape(Z_coordinates , 1 , :)


points = transpose(nodal_coordinates[: , 2 : 4])

println(size(points))

println(element_connectivity[1])

cells = [MeshCell(VTKCellTypes.VTK_QUAD , element_connectivity[1])]


vtk_grid("FEA_solutions" , points, cells) do vtk
    vtk["Concentration" , VTKPointData()]  = concentration_solution[: , 1]
    vtk["time", VTKFieldData()] = 1
end

