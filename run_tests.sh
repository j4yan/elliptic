#! /bin/bash

inputlines[1]='arg_dict = Dict{Any, Any}(
	"operator_type" => "SBPOmega",
    "use_DG" => true,
    "delta_t" => 5.0e-4,
	"t_max" => 50.0,
    "res_tol" => 1e-12,
    "step_tol" => 1e-12,
    "itermax" => 100,'
inputlines[2]='    "order" => '
inputlines[3]='    "smb_name" => "src/solver/poisson/meshfiles'
inputlines[4]='    "dmg_name" => ".null",
	"numBC" => 1,
    "BC1" => [0, 1, 2, 3],
	"Functional" => "volumeAverage",
	"Diffusion" => "identity",
	"BC1_name" => "DirichletTrigonometric",
	"SRC_name" => "SrcTrigonometricIdentityDiffusion",
	"exactSolution" =>"ExactTrigonometric",
    "Flux_name" => "Shahbazi",
	"Cip" => -2.0,
	"run_type" => 5,
    "jac_type" => 2,
	"write_eigs" => false,
    "write_jac" => false,
    "solve" => true,
    "use_jac_precond" => true,
    "real_time" => false
)'

meshfiles[1]="square4x4"
meshfiles[2]="square8x8"
meshfiles[3]="square16x16"
meshfiles[4]="square32x32"
meshfiles[5]="square64x64"
mfile_ext=".smb"

solver_dir=/users/yanj4/.julia/v0.4/PDESolver/src/solver/poisson_solver
meshfiles_dir=$solver_dir/meshfiles
inputfile="input.jl"

for p in `seq 1 4`
do 
	dir="p$p"
	mkdir -p $dir 
	cd $dir
	for file in "${meshfiles[@]}" 
	do
		mkdir -p $file
		cd $file
		mfile="$file"0"$mfile_ext"
		cp $meshfiles_dir/$mfile .
		if [ -e $inputfile ]
		then
			rm -f $inputfile
		fi
		mfile=$file$mfile_ext
		echo "${inputlines[1]}" >> $inputfile	
		echo "${inputlines[2]}$p," >> $inputfile	
		echo "${inputlines[3]}/$mfile\"," >> $inputfile	
		echo "${inputlines[4]}" >> $inputfile	

		startup=\"$solver_dir/startup_poisson.jl\"
		input=\"$PWD/input.jl\"
		echo $startup
		run="julia -e 'include($startup); run_poisson($input)' | tee fout"
		eval $run
		cd ..
	done
	cd ..
done
