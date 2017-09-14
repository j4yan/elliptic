
# Startup function for 1 dof elliptic/parabolic equation
"""
  This function invokes the solver for the elliptic equation, using the
  specified input file

  Inputs:
    input_file: a string containing the path to an input file, or just a file
                name.  If it is just a file name, it is taken to be in the
                users pwd.

  Outputs:
    mesh: an AbstractMesh object
    sbp: the SBP operator used in solving the equation
    eqn: the AbstractSolutionData object used to solve the equation
    opts: the options dictonary

"""
function run_elliptic(input_file::AbstractString)

  mesh, sbp, eqn, opts, pmesh = createObjects(input_file)
  solve_elliptic(mesh, sbp, eqn, opts, pmesh)

  return mesh, sbp, eqn, opts
end


"""
  This function creates and initializes the mesh, sbp, eqn, and opts objects

  Inputs:
    file_name: input file name

  Outputs:
    mesh: an AbstractMesh.  The concrete type is determined by the options
          dictionary
    sbp: an AbstractSBP.  The concrete type is determined by the options
         dictionary
    eqn: an EulerData object
    opts: the options dictionary
    pmesh: mesh used for preconditioning, can be same object as mesh
"""
function createObjects(input_file::AbstractString)

  if !MPI.Initialized()
    MPI.Init()
  end

  opts = read_input(input_file)  # read input file and gets default values
  Tdim = opts["dimensions"]
  dofpernode = 1
  sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, dofpernode)

  # test if the mesh is tangled
  for elem = 1:mesh.numEl
    for n = 1:mesh.numNodesPerElement
      if mesh.jac[n, elem] <= 0.0
        error("Error: negative jacobian in elem-", elem)
      end
    end
  end

  if opts["write_timing"]
    MPI.Barrier(mesh.comm)
    if mesh.myrank == 0
      f = open("timing.dat", "a+")
      println(f, mesh_time)
      close(f)
    end
  end

  eqn = EllipticData_{Tsol, Tres, 2, Tmsh}(mesh, sbp, opts)

  # Initialize the elliptic equation
  init(mesh, sbp, eqn, opts)

  return mesh, sbp, eqn, opts, pmesh
end


"""
  Given fully initialized mesh, sbp, eqn, opts, this function solves
  the elliptic equations.  The 4 object should be obtained from 
  createObjects().
  

  Specifically, it applies an initial condition and invokes a nonlinear
  solver according to the options dictionary.

  Inputs:
    mesh: an AbstractMesh
    sbp: an AbstractSBP
    eqn: an EllipticData
    opts: the options dictionary.  This must be the options dictionary returned
          by createObjects().  Changing values in the options dictionary after
          calling createObjects() results in undefined behavior.
    pmesh: mesh used for preconditioning, can be same object as mesh.
           default value of mesh

"""
function solve_elliptic(mesh::AbstractMesh, sbp, eqn::EllipticData, opts, pmesh=mesh)

  myrank = mesh.myrank

  fill!(eqn.res, 0.0)
  fill!(eqn.res_vec, 0.0)

  q_vec = eqn.q_vec

  saveSolutionToMesh(mesh, q_vec)
  writeVisFiles(mesh, "solution_ic")

  MPI.Barrier( mesh.comm)

  # call_nlsolver(mesh, sbp, eqn, opts, pmesh)
  if opts["run_type"] == 1 || opts["run_type"] == 30 || opts["write_eigs"]
    call_nlsolver(mesh, sbp, eqn, opts, pmesh)
  else
    iterate(mesh, pmesh, sbp, eqn, opts)	
  end

  postproc(mesh, sbp, eqn, opts)

  MPI.Barrier(mesh.comm)
  if opts["finalize_mpi"]
    MPI.Finalize()
  end

  return mesh, sbp, eqn, opts
end  # end function

"""
  This function does post processing.
  Typical post processing includes calculation of errors, norms of important quantities, writing of files. etc.

  **Inputs**:
   * mesh
   * sbp
   * eqn
   * opts
"""
function postproc(mesh, sbp, eqn, opts)

  ##### Do postprocessing ######
  println("\nDoing postprocessing")

  Tsol = eltype(eqn.q)
  Tres = eltype(eqn.res)
  Tmsh = eltype(mesh.dxidx)
  myrank = mesh.myrank

  if haskey(opts, "exact_soln_func")
    l2norm, lInfnorm = calcErrorNorm(mesh, sbp, eqn, opts)
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
    functional_value = zeros(Tsol, 3)
    eqn.functional(mesh, sbp, eqn, opts, functional_value)
    println("functional_V = ", real(functional_value[1]))
    println("functional_D = ", real(functional_value[2]))
    println("functional_N = ", real(functional_value[3]))
    fname = "steady_functional.dat"
    f = open(fname, "w")
    println(f, real(functional_value))
  end

  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  saveSolutionToMesh(mesh, real(eqn.q_vec))
  writeVisFiles(mesh, "solution_done")

  return nothing
end


