using EndoBeams
using DelimitedFiles
using Parameters
using LinearAlgebra
using Dierckx
using StaticArrays

include("utils_stent.jl")
include("crimping.jl")
include("positioning_cl.jl")
include("positioning.jl")
include("deployment.jl")

# -------------------------------------------------------------------------------------------
# Read training data sampling
# -------------------------------------------------------------------------------------------

data = readdlm("stent/input/models_input.txt")
data_sphere = readdlm("stent/input/model_sphere_output.txt")
nsimulations = size(data, 1)
deploy_pos_vec = data[:,end]

# -------------------------------------------------------------------------------------------
# Create stent
# -------------------------------------------------------------------------------------------

@unpack nbWires, rStent, rCrimpedStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern = BraidedStent()
initial_positions_stent, connectivity_stent = compute_bs_geom_given_nbTotalCells(nbWires, rStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern)
initial_positions_stent_mat = reshape(reinterpret(Float64, initial_positions_stent), (3, length(initial_positions_stent)))'

# -------------------------------------------------------------------------------------------
# Settings
# -------------------------------------------------------------------------------------------

nb_iterations = 20
output_dir_crimping = "stent/output3D/outputCrimping3D/"

function compute_case_n(n)

    # -------------------------------------------------------------------------------------------
    # Creation of training dataset
    # ----------------------------------------- --------------------------------------------------

    println("--------------------------------------")
    println("CASE N.$n:")

    output_dir_case = "stent/output3D/case$n/"
    if isdir(output_dir_case)
        rm(output_dir_case, recursive=true)
    end
    mkdir(output_dir_case)
    output_dir_positioning_cl = output_dir_case*"outputPositioningCl3D/"
    mkdir(output_dir_positioning_cl)
    output_dir_positioning = output_dir_case*"outputPositioning3D/"
    mkdir(output_dir_positioning)
    output_dir_deployment = output_dir_case*"outputDeployment3D/"
    mkdir(output_dir_deployment)

    filename_cl = "stent/input/centerline_$n.vtk"
    filename_sdf = "stent/input/sdf_$(n).vtk"

    deploy_pos_100 = deploy_pos_vec[n+1]
    deploy_pos = get_init_pos_deploy_middle(filename_cl, deploy_pos_100, initial_positions_stent, output_dir_crimping)

    # -------------------------------------------------------------------------------------------
    # Geometrical positioning centerline
    # -------------------------------------------------------------------------------------------

    initial_cl_stent = positioning_cl(initial_positions_stent, connectivity_stent, nb_iterations, deploy_pos, filename_cl, output_dir_crimping, output_dir_positioning_cl)
    println("           CENTERLINE POSITIONING COMPLETED")

    # -------------------------------------------------------------------------------------------
    # Physical positioning stent
    # -------------------------------------------------------------------------------------------

    positioning(initial_positions_stent, connectivity_stent, nb_iterations, initial_cl_stent[1], output_dir_crimping, output_dir_positioning_cl, output_dir_positioning)
    println("           POSITIONING COMPLETED")

    crimped_positions_stent = initial_positions_stent .+  read_ics_vec(readdlm(output_dir_crimping * "u.txt"))
    positions_cl =  get_centerline_stent(crimped_positions_stent)
    set_origin_stent!(crimped_positions_stent, positions_cl[1]-initial_cl_stent[1])
    open(output_dir_positioning * "/crimped.txt", "w") do io
       writedlm(io, crimped_positions_stent)
    end
    
    disp = read_ics_vec(readdlm(output_dir_positioning * "u.txt"))
    write_vtk_configuration("stent/output3D/positioned_$n.vtk", crimped_positions_stent + disp, connectivity_stent)

    # -------------------------------------------------------------------------------------------
    # Deployment
    # -------------------------------------------------------------------------------------------

    deployment(initial_positions_stent_mat, initial_positions_stent, connectivity_stent, filename_sdf, output_dir_crimping, output_dir_positioning, output_dir_deployment, output_dir_case)

    disp = read_ics_vec(readdlm(output_dir_case * "u.txt"))
    write_vtk_configuration("stent/output3D/deployed_$n.vtk", initial_positions_stent + disp, connectivity_stent)

end

compute_case_n(1)

