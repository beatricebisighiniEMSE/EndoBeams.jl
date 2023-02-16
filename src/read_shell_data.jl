# This function reads the mesh file and returns: 
# nodal coordinates, and element connectivity.
function read_mesh_file(input_folder)
    
    fid = open(input_folder*"/Preprocessing.dat", "r")

    curr_line = split(readline(fid))
    type_geometry = parse(Int, curr_line[1])
    if type_geometry == 1
        println("using CRIMSON geometry")
    elseif type_geometry == 2
        println("using gmsh geometry")
    else
        error("unknown type of geometry")
    end

    # Read no. of nodes
    curr_line = split(readline(fid))
    nnodes = parse(Int, curr_line[1])
    println("nnodes=",nnodes)

    # Read no. of elements
    curr_line = split(readline(fid))
    nelem = parse(Int, curr_line[1])
    println("nelem=",nelem)
    
    # Read nodal coordinates
    i = 0
    pts = zeros(nnodes*3) # initialize coordinates
    
    # Read coordinates of each node
    while i<nnodes*3
        
        curr_line = split(readline(fid))
        n = length(curr_line)
        
        ini  = i + 1
        iend = i + n
        index = 1
        
        for j in ini:iend
            pts[j] = parse(Float64, curr_line[index]) 
            index +=1
            
        end 
        
        i = iend
        
    end
    
    # Re-arrange coordinates data 
    pts_vec = Vector{Vec3{T}}()

    for i in 1:3:nnodes*3

        x = pts[i]
        y = pts[i+1]
        z = pts[i+2]
        push!(pts_vec, Vec3(x,y,z))
    end 

    # Now, read element connectivity
    i = 0
    elem = zeros(nelem*3)

    while i<nelem*3
        curr_line = split(readline(fid))
        n = length(curr_line)
        
        ini  = i + 1
        iend = i + n
        index = 1
        
        for j in ini:iend
            
            elem[j] = parse(Int, curr_line[index]) # Numbering starts from 1
            index +=1
            
        end 
        
        i = iend
    end


    # Rearrange element connectivity
    elem_vec = Vector{Vec3{Int}}()

    for i in 1:3:3*nelem

        i1 = elem[i]
        i2 = elem[i+1]
        i3 = elem[i+2]
        push!(elem_vec, Vec3(i1,i2,i3))
    end 

     # Read no. of fixed nodes
     curr_line = split(readline(fid))
     n_fixed_dofs = parse(Int, curr_line[1])

    fixed_dofs = zeros(n_fixed_dofs)
    # Now read fixed dofs
    for idof=1:1:n_fixed_dofs
        curr_line = split(readline(fid))
        fixed_dofs[idof] = parse(Int, curr_line[1]) 
    end

    
    close(fid)
    
    return  pts_vec, elem_vec, fixed_dofs, type_geometry
    #return pts_vec


end

# This function reads the fiber directions, if needed.
function read_fiber_directions!(material_parameters,input_folder)
    fid = open(input_folder*"/fibers.dat", "r")

    readline(fid) # Skip first line
    readline(fid) # Skip 2nd line
    readline(fid) # Skip 3rd line

    # Read longitudinal fiber directions
    for ielem=1:1:nelem
        curr_line = split(readline(fid))
        for j=1:1:3 # 3 components
            material_parameters.mat_v_fiber1[j,ielem] = parse(T, curr_line[j]) 
        end
    end
    
    readline(fid) # Skip one line

    # Read circumferential fiber directions
    for ielem=1:1:nelem
        curr_line = split(readline(fid))
        for j=1:1:3 # 3 components
            material_parameters.mat_v_fiber2[j,ielem] = parse(T, curr_line[j]) 
        end
    end

    close(fid)

    return material_parameters

end

# This function reads the input file that contains computational parameters
function read_comp_file(nelem,input_folder,sim_choices,type_geometry)

    fid = open(input_folder*"/comp_parameters.dat", "r")

    # Read material type
    # mat_type=1 (Mooney-Rivlin); mat_type=5 (Four Fibers)
    curr_line = split(readline(fid))
    mat_type = parse(Int, curr_line[1])

    # Read tolerance and max. no of iterations
    curr_line = split(readline(fid))
    sim_choices.comp_tolerance = parse(Float64, curr_line[1])
    sim_choices.comp_maxit = parse(Int, curr_line[2])

    
    
    # If time-dependent, read relevant parameters
    # 1: Steady, 0: time-dependent
    curr_line = split(readline(fid))
    #sim_choices = @set sim_choices.flag_steady_analysis=parse(Int, curr_line[1])
    sim_choices.flag_steady_analysis = parse(Int, curr_line[1])
    if sim_choices.flag_steady_analysis == 0
        tint_tini = parse(Float64, curr_line[2])
        tint_tend = parse(Float64, curr_line[3])
        sim_choices.tint_dt = parse(Float64, curr_line[4])
        tint_ro_inf = parse(Float64, curr_line[5])

        # Calculate coefficients for time-dependent solver
        sim_choices.tint_aM = 1.0-(2.0*tint_ro_inf-1.0)/(tint_ro_inf+1.0)

        sim_choices.tint_aF = 1.0-tint_ro_inf/(tint_ro_inf+1.0)

        sim_choices.tint_gamma = 1.0/2.0 + sim_choices.tint_aM-sim_choices.tint_aF

        sim_choices.tint_beta = (1.0/4.0)*(1+sim_choices.tint_aM - sim_choices.tint_aF)^2.0
        
    end
    
    # Initialize the constructor for material properties
    material_parameters = constructor_material_parameters(mat_type,nelem)

    # Read material properties
    # Read the flag to check whether uniform or regionally varying properties are used
    curr_line = split(readline(fid))
    println(curr_line)
    flag_regional_variation = parse(Int, curr_line[1])

    if flag_regional_variation==0          
        # Read uniform material parameters
        curr_line = split(readline(fid)) # I think that this was missing as material parameters are two lines below the flag , M.A. 02/03/2022
        curr_line = split(readline(fid))
        println(curr_line)
        for ielem in 1:1:nelem
            for i_coeff in 1:1:material_parameters.num_mat_properties
                material_parameters.coefficients[ielem,i_coeff]=parse(Float64, curr_line[i_coeff])
            end

        end

    else
        println("Reading regionally varying parameters from the file")

        # Read number of points on main vessel at which properties are prescribed and number of branches
        curr_line = split(readline(fid))
        num_points_along_main_vessel = parse(Int, curr_line[1])
        num_branches = parse(Int, curr_line[2])    
        
        # Array to store centerline distance of points on main vessel at which properties are prescribed
        path_length_main_vessel = Array{Float64}(undef,num_points_along_main_vessel+2)

        # Array to store regionally varying material values that are read in
        aux_mat_properties = Array{Float64}(undef,num_points_along_main_vessel+num_branches,
                                            material_parameters.num_mat_properties)

        # Call a function to read regionally varying material parameters
        read_regional_variation!(fid,mat_type,path_length_main_vessel,aux_mat_properties,
                                num_points_along_main_vessel,num_branches,material_parameters)
        
    end
    println("Finished reading computational parameters")

    close(fid)

    if type_geometry == 1
        # Array to store the values that are read in
        wall_id = Array{Float64}(undef,nelem)
        # Read wall id for each element
        read_wall_id!(wall_id,nelem,input_folder)

        # Interpolate regionally varying data, if needed
        if flag_regional_variation==1

            # Read path lengths for each element as per the centerline of main vessel
            file_centerline = open(input_folder*"/path_length_centerline_1.dat", "r")
            readline(file_centerline) # Skip first line

            # Array to store the values that are read in
            elemental_path_length = Array{Float64}(undef,nelem)

            for i=1:nelem
                curr_line = split(readline(file_centerline))
                elemental_path_length[i] = parse(Float64, curr_line[1])
            end

            close(file_centerline)

            # Now, interpolate the data and assign it to each element
            # Starting point of main centerline
            path_length_main_vessel[1]=0.0
            # End point of main centerline
            path_length_main_vessel[num_points_along_main_vessel+2]=maximum(elemental_path_length) 

            for ielem in 1:1:nelem # Loop over all elements
                wall_id_current_element = wall_id[ielem]

                if wall_id[ielem] == 1 # This element belongs to main vessel
                    # Centerline this od the current element
                    position_along_main_centerline = elemental_path_length[ielem]
                    # Interpolate material properties for the current element
                    material_parameters = get_properties_main_vessel!(ielem,position_along_main_centerline,
                                            path_length_main_vessel,num_points_along_main_vessel,
                                            aux_mat_properties,material_parameters)
                else # This element belongs to a branch
                    ibranch = wall_id_current_element - 1 # Branch number on which this element lies

                    # Assign this element the material properties of the branch it belongs to
                    for iprop in 1:1:material_parameters.num_mat_properties
                        material_parameters.coefficients[ielem,iprop] = aux_mat_properties[num_points_along_main_vessel+ibranch,iprop]
                    end
                end

            end
        end



        
    end

    
    println("Done assigning material properties to all elements")

    # Read fiber directions
    if material_parameters.mat_type == 5
        material_parameters = read_fiber_directions!(material_parameters,input_folder)
    end

    return material_parameters

end

# This subroutine reads wall ids for different walls in the geometry
function read_wall_id!(wall_id,nelem,input_folder)
    file_wall_id = open(input_folder*"/wall_id.dat", "r")

    readline(file_wall_id) # Skip first line
    for i=1:1:nelem
        curr_line = split(readline(file_wall_id))
        wall_id[i] = parse(Int, curr_line[1]) 

    end


    close(file_wall_id)
end

# This function reads regionally varying material parameters
function read_regional_variation!(fid,mat_type,path_length_main_vessel,aux_mat_properties,
                                  num_points_along_main_vessel,num_branches,material_parameters)


   
    # Read properties at all points along the main vessel
    for i=1:1:num_points_along_main_vessel+num_branches

        curr_line = split(readline(fid))
        
        # Read material properties at the current point along main vessel
        for j=1:1:material_parameters.num_mat_properties
            aux_mat_properties[i,j] = parse(Float64, curr_line[j])
        end
    end

    println("Finished reading material properties")




        
    # Read centerline distances of the points at which material properties are prescribed
    for i=1:1:num_points_along_main_vessel

        curr_line = split(readline(fid))
        
        path_length_main_vessel[i+1] = parse(Float64, curr_line[1])

    end


end

# This function performs linear interpolation of material properties along the main vessel
function get_properties_main_vessel!(ielem,position_along_main_centerline,
                                    path_length_main_vessel,num_points_along_main_vessel,
                                    aux_mat_properties,material_parameters)

    # If element lies in the first block of material properties
    # i.e. between centerline position=0 and the first point of material property assignment
    if position_along_main_centerline <= path_length_main_vessel[2]

        # Assign this element the same material properties as the first prescribed point
        for iprop in 1:1:material_parameters.num_mat_properties
            material_parameters.coefficients[ielem,iprop] = aux_mat_properties[1,iprop]
        end

    # If element lies in the last block of material properties
    # i.e. at centerline position greater or equal to the last point of material property assignment  
    elseif position_along_main_centerline >= path_length_main_vessel[num_points_along_main_vessel+1]

        # Assign this element the same material properties as the last prescribed point
        for iprop in 1:1:material_parameters.num_mat_properties
            material_parameters.coefficients[ielem,iprop] = aux_mat_properties[num_points_along_main_vessel,iprop]
        end

    # For all other elements, linearly interpolate the values between the closest prescribed points on either side of element
    else
        # Identify the closest prescribed points for the current element
        i_position = 0
        for i in 2:1:num_points_along_main_vessel
            if position_along_main_centerline >= path_length_main_vessel[i] && position_along_main_centerline <= path_length_main_vessel[i+1]
                i_position = i
            end  
        end

        # Assign this element the interpolated material properties 
        for iprop in 1:1:material_parameters.num_mat_properties
            slope = (aux_mat_properties[i_position,iprop]-aux_mat_properties[i_position-1,iprop])/
                    (path_length_main_vessel[i_position+1]-path_length_main_vessel[i_position])

            intercept = (aux_mat_properties[i_position-1,iprop]*path_length_main_vessel[i_position+1]-
                         aux_mat_properties[i_position,iprop]*path_length_main_vessel[i_position])/
                         (path_length_main_vessel[i_position+1]-path_length_main_vessel[i_position])

            material_parameters.coefficients[ielem,iprop] = slope*position_along_main_centerline+intercept
        end

    end

    return material_parameters

end

function read_pressure_loading(input_folder)
    
    fid = open(input_folder*"/pressure_loading.dat", "r")

    curr_line = split(readline(fid))
    nsteps = parse(Int, curr_line[1])
 
    pressure = zeros(nsteps) # initialize coordinates
    
    for i = 1:nsteps
        curr_line = split(readline(fid))
        pressure[i] = parse(Float64, curr_line[1]) 
    end

    close(fid)
    
    return  pressure
    #return pts_vec


end
