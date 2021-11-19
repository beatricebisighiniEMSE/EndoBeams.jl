#----------------------------------
# STRUCTURE
#----------------------------------
struct MyNode{T}
   
    i::Int # index
    pos::Vec3{T} # initial position

    # dof: total, displacement, angular
    idof_6::Vec6{Int64}
    idof_disp::Vec3{Int64}
    idof_ang::Vec3{Int64}

    # current configuration of the node (@n+1)
    u::Vec3{T} # displacement
    udt::Vec3{T} # velocity
    udtdt::Vec3{T} # acceleration
    w::Vec3{T} # angle (spin vecotr)
    wdt::Vec3{T} # angular velocity
    wdtdt::Vec3{T} # angular acceleration
    R::Mat33{T} # local rotation matrix
    Delt::Mat33{T} # local rotation matrix variation

    # last configuration of the node (@n)
    u_n::Vec3{T}
    udt_n::Vec3{T}
    udtdt_n::Vec3{T}
    w_n::Vec3{T}
    wdt_n::Vec3{T}
    wdtdt_n::Vec3{T}
    R_n::Mat33{T}
    Delt_n::Mat33{T}

    R_global_to_local::Mat33{T}  # rotation matrix from global  (carthesian) to local (cylindrical) coordinates

end

#----------------------------------
# CONSTRUCTOR
#----------------------------------

"Constructor of the nodes StructArray"
function constructor_nodes(X, u_0, udt_0, udtdt_0, w_0, wdt_0, wdtdt_0, plane, R0=Vector{Mat33{T}}(), Delt0=Vector{Mat33{T}}(), T=Float64) 
  
    if R0 == Vector{Mat33{T}}() && Delt0 == Vector{Mat33{T}}() 
            nodes = StructArray(MyNode{T}(
            i, 
            X[i], 
            Vec6{Int64}(6*(i-1).+(1,2,3,4,5,6)), # Note: use tuples instead of vectors when possible (doesn't allocate)
            Vec3{Int64}(6*(i-1).+(1,2,3)), 
            Vec3{Int64}(6*(i-1).+(4,5,6)), 
            Vec3{T}(u_0[3*(i-1)+1], u_0[3*(i-1)+2], u_0[3*(i-1)+3]), 
            Vec3{T}(udt_0[3*(i-1)+1], udt_0[3*(i-1)+2], udt_0[3*(i-1)+3]), 
            Vec3{T}(udtdt_0[3*(i-1)+1], udtdt_0[3*(i-1)+2], udtdt_0[3*(i-1)+3]), 
            Vec3{T}(w_0[3*(i-1)+1], w_0[3*(i-1)+2], w_0[3*(i-1)+3]), 
            Vec3{T}(wdt_0[3*(i-1)+1], wdt_0[3*(i-1)+2], wdt_0[3*(i-1)+3]), 
            Vec3{T}(wdtdt_0[3*(i-1)+1], wdtdt_0[3*(i-1)+2], wdtdt_0[3*(i-1)+3]), 
            Mat33{T}(1, 0, 0, 0, 1, 0, 0, 0, 1), 
            Mat33{T}(1, 0, 0, 0, 1, 0, 0, 0, 1),
            Vec3{T}(u_0[3*(i-1)+1], u_0[3*(i-1)+2], u_0[3*(i-1)+3]), 
            Vec3{T}(udt_0[3*(i-1)+1], udt_0[3*(i-1)+2], udt_0[3*(i-1)+3]), 
            Vec3{T}(udtdt_0[3*(i-1)+1], udtdt_0[3*(i-1)+2], udtdt_0[3*(i-1)+3]), 
            Vec3{T}(w_0[3*(i-1)+1], w_0[3*(i-1)+2], w_0[3*(i-1)+3]), 
            Vec3{T}(wdt_0[3*(i-1)+1], wdt_0[3*(i-1)+2], wdt_0[3*(i-1)+3]), 
            Vec3{T}(wdtdt_0[3*(i-1)+1], wdtdt_0[3*(i-1)+2], wdtdt_0[3*(i-1)+3]), 
            Mat33{T}(1, 0, 0, 0, 1, 0, 0, 0, 1), 
            Mat33{T}(1, 0, 0, 0, 1, 0, 0, 0, 1),
            compute_local_to_global_matrix(i, X, plane[i]))
        for i in 1:size(X,1))
    elseif  R0 != Vector{Mat33{T}}() && Delt0 == Vector{Mat33{T}}() 
            nodes = StructArray(MyNode{T}(
            i, 
            X[i], 
            Vec6{Int64}(6*(i-1).+(1,2,3,4,5,6)), # Note: use tuples instead of vectors when possible (doesn't allocate)
            Vec3{Int64}(6*(i-1).+(1,2,3)), 
            Vec3{Int64}(6*(i-1).+(4,5,6)), 
            Vec3{T}(u_0[3*(i-1)+1], u_0[3*(i-1)+2], u_0[3*(i-1)+3]), 
            Vec3{T}(udt_0[3*(i-1)+1], udt_0[3*(i-1)+2], udt_0[3*(i-1)+3]), 
            Vec3{T}(udtdt_0[3*(i-1)+1], udtdt_0[3*(i-1)+2], udtdt_0[3*(i-1)+3]), 
            Vec3{T}(w_0[3*(i-1)+1], w_0[3*(i-1)+2], w_0[3*(i-1)+3]), 
            Vec3{T}(wdt_0[3*(i-1)+1], wdt_0[3*(i-1)+2], wdt_0[3*(i-1)+3]), 
            Vec3{T}(wdtdt_0[3*(i-1)+1], wdtdt_0[3*(i-1)+2], wdtdt_0[3*(i-1)+3]), 
            R0[i], 
            Mat33{T}(1, 0, 0, 0, 1, 0, 0, 0, 1),
            Vec3{T}(u_0[3*(i-1)+1], u_0[3*(i-1)+2], u_0[3*(i-1)+3]), 
            Vec3{T}(udt_0[3*(i-1)+1], udt_0[3*(i-1)+2], udt_0[3*(i-1)+3]), 
            Vec3{T}(udtdt_0[3*(i-1)+1], udtdt_0[3*(i-1)+2], udtdt_0[3*(i-1)+3]), 
            Vec3{T}(w_0[3*(i-1)+1], w_0[3*(i-1)+2], w_0[3*(i-1)+3]), 
            Vec3{T}(wdt_0[3*(i-1)+1], wdt_0[3*(i-1)+2], wdt_0[3*(i-1)+3]), 
            Vec3{T}(wdtdt_0[3*(i-1)+1], wdtdt_0[3*(i-1)+2], wdtdt_0[3*(i-1)+3]), 
            R0[i], 
            Mat33{T}(1, 0, 0, 0, 1, 0, 0, 0, 1),
            compute_local_to_global_matrix(i, X, plane[i]))
        for i in 1:size(X,1))
    else
            nodes = StructArray(MyNode{T}(
            i, 
            X[i], 
            Vec6{Int64}(6*(i-1).+(1,2,3,4,5,6)), # Note: use tuples instead of vectors when possible (doesn't allocate)
            Vec3{Int64}(6*(i-1).+(1,2,3)), 
            Vec3{Int64}(6*(i-1).+(4,5,6)), 
            Vec3{T}(u_0[3*(i-1)+1], u_0[3*(i-1)+2], u_0[3*(i-1)+3]), 
            Vec3{T}(udt_0[3*(i-1)+1], udt_0[3*(i-1)+2], udt_0[3*(i-1)+3]), 
            Vec3{T}(udtdt_0[3*(i-1)+1], udtdt_0[3*(i-1)+2], udtdt_0[3*(i-1)+3]), 
            Vec3{T}(w_0[3*(i-1)+1], w_0[3*(i-1)+2], w_0[3*(i-1)+3]), 
            Vec3{T}(wdt_0[3*(i-1)+1], wdt_0[3*(i-1)+2], wdt_0[3*(i-1)+3]), 
            Vec3{T}(wdtdt_0[3*(i-1)+1], wdtdt_0[3*(i-1)+2], wdtdt_0[3*(i-1)+3]), 
            R0[i], 
            Delt0[i],
            Vec3{T}(u_0[3*(i-1)+1], u_0[3*(i-1)+2], u_0[3*(i-1)+3]), 
            Vec3{T}(udt_0[3*(i-1)+1], udt_0[3*(i-1)+2], udt_0[3*(i-1)+3]), 
            Vec3{T}(udtdt_0[3*(i-1)+1], udtdt_0[3*(i-1)+2], udtdt_0[3*(i-1)+3]), 
            Vec3{T}(w_0[3*(i-1)+1], w_0[3*(i-1)+2], w_0[3*(i-1)+3]), 
            Vec3{T}(wdt_0[3*(i-1)+1], wdt_0[3*(i-1)+2], wdt_0[3*(i-1)+3]), 
            Vec3{T}(wdtdt_0[3*(i-1)+1], wdtdt_0[3*(i-1)+2], wdtdt_0[3*(i-1)+3]), 
            R0[i], 
            Delt0[i],
            compute_local_to_global_matrix(i, X, plane[i]))
        for i in 1:size(X,1))
    end

    return nodes

end

#----------------------------------
# UTILS 
#----------------------------------

#  Compute rotation matrix from cylindrical to carthesian coordinates
function compute_local_to_global_matrix(i, X, plane)

    if plane == "xy"

        Xi = X[i]
        thetai = atan(Xi[2], Xi[1])  
        return Mat33{T}(cos(thetai), -sin(thetai), 0,  sin(thetai), cos(thetai), 0, 0, 0, 1)

    elseif plane == "yz"

        Xi = X[i]
        thetai = atan(Xi[3], Xi[2])
        return Mat33{T}(1, 0, 0, 0, cos(thetai), -sin(thetai), 0, sin(thetai), cos(thetai))
    
    elseif plane == "xz"

        Xi = X[i]
        thetai = atan(Xi[3], Xi[1])
        return Mat33{T}(cos(thetai), 0, -sin(thetai), 0, 1, 0, sin(thetai), 0, cos(thetai))

    end 
    
end