arg_dict = Dict{Any, Any}(
  # "operator_type" => "SBPGamma",
  "operator_type" => "SBPOmega",
  "use_DG" => true,
  "delta_t" => 1.0e-3,
  "t_max" => 5.0,
  "res_tol" => 1e-12,
  "step_tol" => 1e-12,
  "itermax" => 1,
  "order" => 3,
  # "smb_name" => "/users/yanj4/Downloads/meshfiles/square4x4.smb",
  # "smb_name" => "/users/yanj4/Downloads/meshfiles/square8x8.smb",
  "smb_name" => "/users/yanj4/Downloads/meshfiles/square32x32.smb",
  "dmg_name" => ".null",
  "numBC" => 2,
  "BC1" => [0, 1, 2],
  "BC2" => [3],
  "Functional" => "volumeAverage",

  "exactSolution" =>"ExactExpTrig",
  "BC1_name" => "DirichletExpTrig",
  "BC2_name" => "NeumannExpTrig",
  "Diffusion" => "poly2nd", # lambda = [10 & 0 \\ 0 & 10]
  "SRC_name" => "SrcExpTrigPoly2ndDiffn",
  "exactFunctional" => 1.846230857168755189823348538541e-02,

  # "exactSolution" =>"ExactTrig",
  # "BC1_name" => "DirichletTrig",
  # "Diffusion" => "poly2nd", # lambda = [x^2+1 & xy \\ xy & y^2+1]
  # "SRC_name" => "SrcTrigPoly2ndDiffn",

  # "SRC_name" => "SRC0",

  # "exactSolution" =>"ExactTrig",
  # "BC1_name" => "DirichletTrig",
  # "Diffusion" => "poly6th", # lambda = [x^6+1 & xy \\ xy & y^6+1]
  # "SRC_name" => "SrcTrigPoly6thDiffn",

  # "exactSolution" =>"ExactPoly2nd",
  # "Diffusion" => "DiffnPoly0th"
  # "BC1_name" => "DirichletPolynial2nd",
  # "SRC_name" => "SrcPoly2nd",

  # "SRC_name" => "SRC0",		# This is for energy stability test 

  # "Flux_name" => "SAT",			# SAT-BR2
  "Flux_name" => "SAT0",		# SAT-SIPG
  # "Flux_name" => "Shahbazi",	# Shahbazi
  "Cip" => -3.0,
  "run_type" => 5,
  "jac_type" => 2,
  "jac_method" => 2,
  "write_eigs" => false,
  "write_jac" => false,
  "write_conditioning" => false,
  "write_energy" => true,
  "solve" => true,
  # "use_jac_precond" => true,
  # "TimeAdvance" => "BDF2",
  # "TimeAdvance" => "CN",
  # "TimeAdvance" => "SDIRK4",
  "real_time" => false
)


