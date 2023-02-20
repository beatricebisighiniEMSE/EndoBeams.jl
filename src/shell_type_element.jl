#----------------------------------
# STRUCTURE
#----------------------------------

struct ShellMyElement{T}
    ielem::Int # element number
    node1::Int # index of node 1   
    node2::Int # index of node 2
    node3::Int # index of node 3
    sparsity_map:: Array{Int,1} # sparsity map from local indices (beam matrix) to global indices (global sparse matrix) -> computed in the constructor of the sparse matrices
    F2D::Mat33{T} # 3x3 deformation gradient that does not include the out of plane component
    F3D::Mat33{T} # 3x3 deformation gradient
    J::T #Surface Jacobian
    SPK::Mat33{T} # 3x3 Second Piola - Kirchhoff stress tensor
    Cauchy::Mat33{T} # 3x3 Cauchy Stress tensor
    CauchyDev::Mat33{T} # 3x3 Cauchy Stress tensor Deviatoric
    CauchyCyl::Mat33{T} # 3x3 Cauchy Stress tensor cylindrical coordinates
    CauchyCylDev::Mat33{T} # 3x3 Cauchy Stress tensor cylindrical coordinates Deviatoric
end

#----------------------------------
# CONSTRUCTORS
#----------------------------------

"Constructor of the beams (StructArray)"
function shell_constructor_elements(allnodes, connectivity)
        beams = StructArray(shell_constructor_single_element(
                i,
                allnodes[connectivity[i][1]], 
                allnodes[connectivity[i][2]], 
                allnodes[connectivity[i][3]]) 
            for i in 1:size(connectivity,1))
    return beams
end 

"Constructor of the one beam (MyBeam)"
function shell_constructor_single_element(i, node1, node2, node3,  T=Float64)
    
    ielem = i

    i1 = node1.i
    
    i2 = node2.i

    i3 = node3.i

    return ShellMyElement{T}(ielem, i1, i2, i3, zeros(Int, 324),zeros(T,3,3),zeros(T,3,3),0,zeros(T,3,3),zeros(T,3,3),zeros(T,3,3),zeros(T,3,3),zeros(T,3,3))
    
end 



