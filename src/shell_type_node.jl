#----------------------------------
# STRUCTURE
#----------------------------------
struct ShellMyNode{T}
   
    i::Int # index
    pos::Vec3{T} # Initial nodal position
    u::Vec3{T} #nodal displacement solution
    v::Vec3{T} #nodal displacement solution
    a::Vec3{T} #nodal displacement solution
    pos_cur::Vec3{T} # Current nodal position
    idof::Vec3{Int}

    # # dof: total, displacement, angular
    # idof_6::Vec6{Int64}
    # idof_disp::Vec3{Int64}
    # idof_ang::Vec3{Int64}

    # # current configuration of the node (@n+1)
    # u::Vec3{T} # displacement
    # udt::Vec3{T} # velocity
    # udtdt::Vec3{T} # acceleration
    # w::Vec3{T} # angle (spin vecotr)
    # wdt::Vec3{T} # angular velocity
    # wdtdt::Vec3{T} # angular acceleration
    # R::Mat33{T} # local rotation matrix
    # Delt::Mat33{T} # local rotation matrix variation

    # # last configuration of the node (@n)
    # u_n::Vec3{T}
    # udt_n::Vec3{T}
    # udtdt_n::Vec3{T}
    # w_n::Vec3{T}
    # wdt_n::Vec3{T}
    # wdtdt_n::Vec3{T}
    # R_n::Mat33{T}
    # Delt_n::Mat33{T}

    # R_global_to_local::Mat33{T}  # rotation matrix from global  (carthesian) to local (cylindrical) coordinates

end

#----------------------------------
# CONSTRUCTOR
#----------------------------------
function shell_constructor_nodes(X, u, T=Float64) 
    
    nodes = StructArray(ShellMyNode{T}(
            i, 
            X[i], 
            u[i],
            u[i],
            u[i],
            X[i],
            [3*i-2, 3*i-1, 3*i])
        for i in 1:size(X,1))

    return nodes

end
