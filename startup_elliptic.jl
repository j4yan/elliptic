push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"),     "src/solver/elliptic"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"),     "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"),     "src/Debugging"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"),     "src/Utils"))

using PDESolver
using ODLCommonTools
using PdePumiInterface
using SummationByParts
using EllipticEqnMod
using NonlinearSolvers
using ArrayViews
using Utils
using MPI
using Debug

function run(fin::AbstractString)
  include(joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl"))

  if !MPI.Initialized()
    MPI.Init()
  end

  println("ARGS = ", ARGS)
  println("size(ARGS) = ", size(ARGS))
  opts = read_input(fin)
  dofpernode = 1
  sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, dofpernode)

  for elem = 1:mesh.numEl
    for n = 1:mesh.numNodesPerElement
      if mesh.jac[n, elem] <= 0.0
        error("Error: negative jacobian")
      end
    end
  end

  # myrank = mesh.myrank
  # TODO: input argument for dofpernode

  # create Elliptic equation
  # var_type = opts["variable_type"]
  # @bp
  eqn = EllipticData_{Tsol, Tres, 2, Tmsh}(mesh, sbp, opts)

  init(mesh, sbp, eqn, opts, pmesh)

  res_vec = eqn.res_vec
  q_vec = eqn.q_vec

  iterate(mesh, pmesh, sbp, eqn, opts)	
  # solve_euler(mesh, sbp, eqn, opts, pmesh)
  # call_nlsolver(mesh, sbp, eqn, opts, pmesh)

  if haskey(opts, "exactSolution")
    l2norm, lInfnorm = calcErrorL2Norm(mesh, sbp, eqn, opts)
    println("L2Norm = ", l2norm)
    println("LinfNorm = ", lInfnorm)
    fname = "l2norm.dat"
    f = open(fname, "w")
    println(f, l2norm)
  end

  if haskey(opts, "Functional")
    if haskey(opts, "exactFunctional")
      exactFunctional = opts["exactFunctional"]
    end
    functional_value = Array(Tsol, mesh.numDofPerNode)
    eqn.functional(mesh, sbp, eqn, opts, functional_value)
    println("functional = ", abs(real(functional_value[1]) - exactFunctional))
    fname = "functional.dat"
    f = open(fname, "w")
    println(f, abs(real(functional_value[1]) - exactFunctional))
  end

  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  saveSolutionToMesh(mesh, real(eqn.q_vec))
  writeVisFiles(mesh, "solution_done")
end
