using DelimitedFiles
using Parameters
using LinearAlgebra
using Dierckx
using StaticArrays
using LinearAlgebra
using StructArrays
using EndoBeams
using SignedDistanceField

# cd("stent")
include("utils_stent.jl")
include("crimping.jl")
include("positioning_cl.jl")
include("positioning.jl")
include("deployment.jl")

# -------------------------------------------------------------------------------------------
# Create stent
# -------------------------------------------------------------------------------------------

@unpack nbWires, rStent, rCrimpedStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern = Wallstent()
initial_positions_stent, connectivity_stent = compute_bs_geom_given_nbTotalCells(nbWires, rStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern)
initial_positions_stent_mat = reshape(reinterpret(Float64, initial_positions_stent), (3, length(initial_positions_stent)))'

output_dir_crimping = "stent/output3D/outputCrimping3D/"
crimping(rStent, rCrimpedStent, initial_positions_stent, connectivity_stent, output_dir_crimping)

output_dir_positioning_cl = "stent/output3D/outputPositioningCl3D/"
output_dir_positioning = "stent/output3D/outputPositioning3D/"
output_dir_deployment = "stent/output3D/outputDeployment3D/"

nb_iterations = 20
deploy_pos =  get_init_pos_deploy_middle(filename_cl, initial_positions_stent, output_dir_crimping)

filename_surf = "stent/input/stenosis_model.stl"
filename_cl = "stent/input/stenosis_cl.vtk"
filename_sdf = "stent/input/stenosis_sdf.vtk"

# -------------------------------------------------------------------------------------------
# Geometrical positioning centerline
# -------------------------------------------------------------------------------------------

initial_cl_stent = positioning_cl(initial_positions_stent, connectivity_stent, nb_iterations, deploy_pos, filename_cl, output_dir_crimping, output_dir_positioning_cl)

# -------------------------------------------------------------------------------------------
# Physical positioning stent
# -------------------------------------------------------------------------------------------

positioning(initial_positions_stent, connectivity_stent, nb_iterations, initial_cl_stent[1], output_dir_crimping, output_dir_positioning_cl, output_dir_positioning)

crimped_positions_stent = initial_positions_stent .+  read_ics_vec(readdlm(output_dir_crimping * "u.txt"))
positions_cl =  get_centerline_stent(crimped_positions_stent)
set_origin!(crimped_positions_stent, positions_cl[1]-initial_cl_stent[1])
open(output_dir_positioning * "/crimped.txt", "w") do io
    writedlm(io, crimped_positions_stent)
end

# -------------------------------------------------------------------------------------------
# Deployment
# -------------------------------------------------------------------------------------------

if isfile(filename_sdf)
    rm(filename_sdf)
end
s, x, y, z = sdfgen(filename_surf, 0.1; padding = 1, acceleration=:KDTree)
write_vtk_general_structured_mesh(filename_sdf, -s, 0.1, x, y, z)

deployment(initial_positions_stent_mat, initial_positions_stent, connectivity_stent, filename_sdf, output_dir_crimping, output_dir_positioning, output_dir_deployment)

