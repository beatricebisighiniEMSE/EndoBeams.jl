#----------------------------------
# STRUCTURE
#----------------------------------
struct Node{T}
   
    i::Int # index
    X₀::Vec3{T} # initial position

    # dof: total, displacement, angular
    idof_6::Vec6{Int}
    idof_disp::Vec3{Int}
    idof_ang::Vec3{Int}

    # current configuration of the node (@n+1)
    u::Vec3{T} # displacement
    u̇::Vec3{T} # velocity
    ü::Vec3{T} # acceleration
    w::Vec3{T} # angle (spin vecotr)
    ẇ::Vec3{T} # angular velocity
    ẅ::Vec3{T} # angular acceleration
    R::Mat33{T} # local rotation matrix
    ΔR::Mat33{T} # local rotation matrix variation

    # last configuration of the node (@n)
    uⁿ::Vec3{T}
    u̇ⁿ::Vec3{T}
    üⁿ::Vec3{T}
    wⁿ::Vec3{T}
    ẇⁿ::Vec3{T}
    ẅⁿ::Vec3{T}
    Rⁿ::Mat33{T}
    ΔRⁿ::Mat33{T}

    R_global_to_local::Mat33{T}  # rotation matrix from global  (carthesian) to local (cylindrical) coordinates

end

#----------------------------------
# CONSTRUCTOR
#----------------------------------

"""
nodes = constructor_nodes(X, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, plane, Rₑ⁰=nothing, T=Float64) 

Constructor of the nodes StructArray:
- `X`: nodes StructArray (created with constructor_nodes);
- `u⁰`: initial displacements;
- `u̇⁰`: initial velocities;
- `ü⁰`: initial accelerations;
- `w⁰`: initial rotations;
- `ẇ⁰`: initial rotation velocities;
- `ẅ⁰`: initial rotation acceleration;
- `plane`: plane used for the conversin in cylindrical coordinates in case of BCs expressed in cylindrical coordinates.
- `Rₑ⁰`: (not mandatory) initial rotation of the nodes.

Returns a StructArray{Node}, structure containing the information of the nodes. 
"""
function constructor_nodes(X, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, plane, Rₑ⁰=nothing, T=Float64) 

    nodes = StructArray(Node{T}(
            i, 
            X[i], 
            Vec6{Int}(6*(i-1).+(1,2,3,4,5,6)),
            Vec3(6*(i-1).+(1,2,3)), 
            Vec3(6*(i-1).+(4,5,6)), 
            Vec3(u⁰[3*(i-1)+1], u⁰[3*(i-1)+2], u⁰[3*(i-1)+3]), 
            Vec3(u̇⁰[3*(i-1)+1], u̇⁰[3*(i-1)+2], u̇⁰[3*(i-1)+3]), 
            Vec3(ü⁰[3*(i-1)+1], ü⁰[3*(i-1)+2], ü⁰[3*(i-1)+3]), 
            Vec3(w⁰[3*(i-1)+1], w⁰[3*(i-1)+2], w⁰[3*(i-1)+3]), 
            Vec3(ẇ⁰[3*(i-1)+1], ẇ⁰[3*(i-1)+2], ẇ⁰[3*(i-1)+3]), 
            Vec3(ẅ⁰[3*(i-1)+1], ẅ⁰[3*(i-1)+2], ẅ⁰[3*(i-1)+3]), 
            Rₑ⁰ isa AbstractVector ? Rₑ⁰[i] : ID3, 
            ID3,
            Vec3(u⁰[3*(i-1)+1], u⁰[3*(i-1)+2], u⁰[3*(i-1)+3]), 
            Vec3(u̇⁰[3*(i-1)+1], u̇⁰[3*(i-1)+2], u̇⁰[3*(i-1)+3]), 
            Vec3(ü⁰[3*(i-1)+1], ü⁰[3*(i-1)+2], ü⁰[3*(i-1)+3]), 
            Vec3(w⁰[3*(i-1)+1], w⁰[3*(i-1)+2], w⁰[3*(i-1)+3]), 
            Vec3(ẇ⁰[3*(i-1)+1], ẇ⁰[3*(i-1)+2], ẇ⁰[3*(i-1)+3]), 
            Vec3(ẅ⁰[3*(i-1)+1], ẅ⁰[3*(i-1)+2], ẅ⁰[3*(i-1)+3]), 
            Rₑ⁰ isa AbstractVector ? Rₑ⁰[i] : ID3,
            ID3,
            compute_local_to_global_matrix(i, X, plane[i])) for i in 1:size(X,1))

    return nodes

end



#----------------------------------
# UTILS 
#----------------------------------

#  Compute rotation matrix from cylindrical to carthesian coordinates
function compute_local_to_global_matrix(i, X, plane)

    if plane == "xy"

        Xi = X[i]
        Θi = atan(Xi[2], Xi[1])  
        return Mat33(cos(Θi), -sin(Θi), 0,  sin(Θi), cos(Θi), 0, 0, 0, 1)

    elseif plane == "yz"

        Xi = X[i]
        Θi = atan(Xi[3], Xi[2])
        return Mat33(1, 0, 0, 0, cos(Θi), -sin(Θi), 0, sin(Θi), cos(Θi))
    
    elseif plane == "xz"

        Xi = X[i]
        Θi = atan(Xi[3], Xi[1])
        return Mat33(cos(Θi), 0, -sin(Θi), 0, 1, 0, sin(Θi), 0, cos(Θi))

    end 
    
end
