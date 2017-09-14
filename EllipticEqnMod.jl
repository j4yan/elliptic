module EllipticEqnMod

push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))

using PDESolver  # setup LOAD_PATH to find all PDESolver components
using ArrayViews
using ODLCommonTools
using SummationByParts
using PdePumiInterface
using NonlinearSolvers
using Utils
using MPI
using Input

import ODLCommonTools.sview
export AbstractEllipticData, EllipticData, EllipticData_, run_elliptic


abstract AbstractEllipticData{Tsol, Tres} <: AbstractSolutionData{Tsol, Tres}
abstract EllipticData{Tsol, Tres, Tdim} <: AbstractEllipticData{Tsol, Tres}

# now that EllipticData is declared, include other files that use it
include(joinpath(Pkg.dir("PDESolver"), "src/solver/debug.jl"))
include("elliptic.jl")
include("time_advance.jl")
include("bc.jl")
include("flux.jl")
include("source.jl")
include("output.jl")
include("exact_soln.jl")
include("diffusion.jl")
include("functional.jl")
include("types.jl")
include("startup_func.jl")  # function to invoke the solver


global const PhysicsName = "Elliptic"
register_physics(PhysicsName, EllipticEqnMod, run_elliptic)

end # end module
