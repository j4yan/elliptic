push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"),     "src/solver/elliptic"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"),     "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"),     "src/Debugging"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"),     "src/Utils"))

using ODLCommonTools
using PdePumiInterface
using SummationByParts
using EllipticEqnMod
using NonlinearSolvers
using ArrayViews
using Utils
using MPI
using Debug
@debug function run_elliptic(fin::ASCIIString)


    include(joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl"))

    if !MPI.Initialized()
         MPI.Init()
    end
    # @bp
    println("ARGS = ", ARGS)
    println("size(ARGS) = ", size(ARGS))
    opts = read_input(fin)

    # flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
    flag = opts["run_type"]

    # timestepping parameters
    delta_t = opts["delta_t"]

    if flag == 1
        Tmsh = Float64
        Tsbp = Float64
        Tsol = Float64
        Tres = Float64
    elseif flag == 4  # use Newton method using finite difference
        Tmsh = Float64
        Tsbp = Float64
        Tsol = Float64
        Tres = Float64
    elseif flag == 5  # use complex step dR/du
        Tmsh = Float64
        Tsbp = Float64
        Tsol = Complex128
        Tres = Complex128
    else
        error("We should never get here!\n")
    end

    order = opts["order"]
	if opts["operator_type"] == "SBPOmega"
		reorder = false
		internal = true
		shape_type = 2
	elseif opts["operator_type"] == "SBPGamma"
		reorder = false
		internal = false
		shape_type = 3
	else
		opType = opts["operator_type"]
		throw(ArgumentError("unrecognized operator type $opType"))
		error("unrecognized operator type $opType")
	end
    sbp = TriSBP{Tsbp}(degree=order, reorder=reorder, internal=internal)  # create linear sbp operator
    ref_verts = [-1.0 1.0 -1.0; -1.0 -1.0 1.0]
    interp_op = SummationByParts.buildinterpolation(sbp, ref_verts)
    sbpface = TriFace{Float64}(order, sbp.cub, ref_verts.')

    # create linear mesh with 1 dof per node
    dmg_name = opts["dmg_name"]
    smb_name = opts["smb_name"]
    println(dmg_name)
    println(smb_name)
    println("constructing DG mesh")
    # @bp
    mesh = PumiMeshDG2{Tmsh}(dmg_name, smb_name, order, sbp, opts, interp_op, sbpface;
                dofpernode=1, coloring_distance=opts["coloring_distance"], shape_type=shape_type)
    println("constructing DG mesh completed")
    if (opts["jac_type"] == 3 || opts["jac_type"] == 4) && opts["use_jac_precond"]
        pmesh = PumiMeshDG2Preconditioning(mesh, sbp, opts;
        coloring_distance=opts["coloring_distance_prec"])
    else
        pmesh = mesh
    end
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

    #
    # do some test
    #
	# interpolateFace(mesh, sbp, eqn, opts, mesh.coords, eqn.xy_face)
    # evalElliptic(mesh, sbp, eqn, opts)
    # testOperatorGradientInterpolation(mesh, sbp, eqn, opts)

    # if flag == 1
		# rk4(evalElliptic, opts["delta_t"], opts["t_max"], mesh, sbp, eqn, opts,
            # res_tol=opts["res_abstol"], real_time=opts["real_time"])
    # elseif flag == 4 || flag == 5
        # newton(evalElliptic, mesh, sbp, eqn, opts, pmesh, itermax=opts["itermax"],
                   # step_tol=opts["step_tol"], res_abstol=opts["res_abstol"],
                   # res_reltol=opts["res_reltol"], res_reltol0=opts["res_reltol0"])
    # else
        # error("We should never get here!\n")
    # end
	iterate(mesh, pmesh, sbp, eqn, opts)	

    if haskey(opts, "exactSolution")
        l2norm, lInfnorm = calcErrorL2Norm(mesh, sbp, eqn, opts)
        fname = "l2norm.dat"
        println("l2norm = ", l2norm)
        f = open(fname, "w")
        println(f, l2norm, " ", lInfnorm)
    end

	if haskey(opts, "Functional")
		functional_value = Array(Tsol, mesh.numDofPerNode)
		eqn.functional(mesh, sbp, eqn, opts, functional_value)
        fname = "functional.dat"
		println("functional = ", abs(real(functional_value)))
        f = open(fname, "w")
		println(f, abs(real(functional_value)))
	end
	
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    saveSolutionToMesh(mesh, real(eqn.q_vec))
    writeVisFiles(mesh, "solution_done")
end
