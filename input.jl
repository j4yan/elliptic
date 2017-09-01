arg_dict = Dict{Any, Any}(
"operator_type" => "SBPOmega",
"use_DG" => true,
"delta_t" => 1.0e-5,
"t_max" => 50.0,
"res_tol" => 1e-13,
"step_tol" => 1e-13,
"itermax" => 1,
"order" => 4,
# "smb_name" => "/users/yanj4/Downloads/meshfiles/square4x4.smb",
# "smb_name" => "/users/yanj4/Downloads/meshfiles/square8x8.smb",
# "smb_name" => "/users/yanj4/Downloads/meshfiles/square12x12.smb",
"smb_name" => "/users/yanj4/Downloads/meshfiles/square16x16.smb",
# "smb_name" => "/users/yanj4/Downloads/meshfiles/square32x32.smb",
# "smb_name" => "/users/yanj4/Downloads/meshfiles/square64x64.smb",
"dmg_name" => ".null",
"SRC_name" => "SrcZeros",

# "numBC" => 1,
# "BC1" => [0,1,2,3],
# "BC1_name" => "DirichletHicken2011",
# "SRC_name" => "SrcHicken2011",
# "exactSolution" =>"ExactHicken2011",
# "Diffusion" => "DiffnHicken2011", 
# "Functional" => "FuncHicken2011",

"numBC" => 1,
"BC1" => [0,1,2,3],
"Diffusion" => "poly2nd", 
"BC1_name" => "DirichletTrig",
"SRC_name" => "SrcTrigPoly2ndDiffnSlightUnsteady2",
"exactSolution" =>"ExactTrigSlightUnsteady2",
"IC" =>"ExactTrigSlightUnsteady2",
"Functional" => "volumeEnergy",
# "write_energy" => true,

"Flux_name" => "SAT0",
"Cip" => -1.0,
"run_type" => 5,
"jac_type" => 2,
"jac_method" => 2,
"write_eigs" => false,
"write_jac" => false,
"solve" => true,
# "TimeAdvance" => "CN",
# "real_time" => true,
"real_time" => false,
)
