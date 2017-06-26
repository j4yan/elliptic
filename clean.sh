#! /bin/bash

# rm -r newton_*
# rm -r mesh_complete
# rm -r solution_ic
# rm -r output_init
# rm -r residual
# rm arg_dict_output.jl
# rm *.dat
# # rm fout
# rm *.out

# rm -r solution_done
# rm -r solution_error
# if [ $1 == "all" ]; then
  # rm -r adjoint_Cd
  # rm -r adjoint_Cl
# fi

find . -iname "newton_*" -exec rm -r {} \;
find . -iname "mesh_complete" -exec rm -r {} \;
find . -iname "solution_ic" -exec rm -r {} \;
find . -iname "adjoint_Cd" -exec rm -r {} \;
find . -iname "adjoint_Cl" -exec rm -r {} \;
find . -iname "solution_done" -exec rm -r {} \;
find . -iname "solution_error" -exec rm -r {} \;
find . -iname "residual" -exec rm -r {} \;
find . -iname "output_init" -exec rm -r {} \;
find . -iname "arg_dict_output.jl" -exec rm -r {} \;
find . -iname "IC_0.dat" -exec rm -r {} \;
find . -iname "load_balance_0.dat" -exec rm -r {} \;
find . -iname "meshlog_0.dat" -exec rm -r {} \;
find . -iname "solution_ic.dat" -exec rm -r {} \;
find . -iname "nrm_face.dat" -exec rm -r {} \;
find . -iname "l2norm.dat" -exec rm -r {} \;
find . -iname "log_0.dat" -exec rm -r {} \;
find . -iname "convergence.dat" -exec rm -r {} \;
find . -iname "energy.dat" -exec rm -r {} \;
find . -iname "fout" -exec rm -r {} \;
find . -iname "functional.dat" -exec rm -r {} \;
find . -iname "nohup.out" -exec rm -r {} \;
