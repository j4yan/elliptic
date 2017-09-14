
function BDF2{Tmsh, Tsol, Tres, Tdim, Tsbp}(mesh::AbstractMesh{Tmsh}, 
                                            sbp::AbstractSBP{Tsbp}, 
                                            eqn::EllipticData{Tsol, Tres, Tdim}, 
                                            opts) 		
  dt = opts["delta_t"]
  for el= 1:mesh.numEl
    jac = sview(mesh.jac, :, el)
    q = sview(eqn.q, :, :, el )
    q1 = sview(eqn.q1, :, :, el )
    q2 = sview(eqn.q2, :, :, el )
    for n = 1:mesh.numNodesPerElement 
      for dof = 1:mesh.numDofPerNode
        dqdt = (1.5*q[dof, n] - 2.0*q1[dof, n] + 0.5*q2[dof, n])/dt
        eqn.res[dof, n, el] += sbp.w[n]/jac[n]*dqdt
      end
    end
  end
  return nothing
end

function CN{Tmsh, Tsol, Tres, Tdim, Tsbp}(mesh::AbstractMesh{Tmsh}, 
                                          sbp::AbstractSBP{Tsbp}, 
                                          eqn::EllipticData{Tsol, Tres, Tdim}, 
                                          opts) 		
  dt = opts["delta_t"]

  for el= 1:mesh.numEl
    for n = 1:mesh.numNodesPerElement 
      for dof = 1:mesh.numDofPerNode
        eqn.res[dof, n, el] *= 0.5 
        eqn.res[dof, n, el] += 0.5*eqn.res1[dof, n, el] 
      end
    end
  end

  for el= 1:mesh.numEl
    jac = sview(mesh.jac, :, el)
    q = sview(eqn.q, :, :, el )
    q1 = sview(eqn.q1, :, :, el )
    for n = 1:mesh.numNodesPerElement 
      for dof = 1:mesh.numDofPerNode
        dqdt = (q[dof, n] - q1[dof, n])/dt
        eqn.res[dof, n, el] += sbp.w[n]/jac[n]*dqdt
      end
    end
  end
  return nothing
end

function SDIRK4{Tmsh, Tsol, Tres, Tdim, Tsbp}(mesh::AbstractMesh{Tmsh}, 
                                              sbp::AbstractSBP{Tsbp}, 
                                              eqn::EllipticData{Tsol, Tres, Tdim}, 
                                              opts) 		
  dt = opts["delta_t"]
  c = zeros(Float64, 5, 6)
  c[1,1] = -4.0
  c[1,2] = 4.0
  c[2,1] = 4.0
  c[2,2] = -8.0
  c[2,3] = 4.0
  c[3,1] = 52.0/25.0
  c[3,2] = -168.0/25.0
  c[3,3]= 16.0/25.0
  c[3,4]= 4.0
  c[4,1] = 16.0/17.0
  c[4,2] = -89.0/17.0
  c[4,3] = 25.0/34.0
  c[4,4] = -15.0/34.0
  c[4,5] = 4.0
  c[5,1] = -28.0/3.0
  c[5,2] = 37.0/3.0
  c[5,3] = 103.0/6.0
  c[5,4] = -275.0/2.0
  c[5,5] = 340.0/3.0
  c[5,6] = 4.0
  println(eqn.istage)
  println()
  for stage = 1:eqn.istage+1
    for el= 1:mesh.numEl
      jac = sview(mesh.jac, :, el)
      qs = sview(eqn.q_irk, :, :, el, stage)
      for n = 1:mesh.numNodesPerElement 
        for dof = 1:mesh.numDofPerNode
          eqn.res[dof, n, el] += sbp.w[n]/jac[n]*qs[dof,n]*c[eqn.istage, stage]/dt
        end
      end
    end
  end
  return nothing
end

function evalRes{Tmsh, Tsol, Tres, Tdim, Tsbp}(mesh::AbstractMesh{Tmsh}, 
                                               sbp::AbstractSBP{Tsbp}, 
                                               eqn::EllipticData{Tsol, Tres, Tdim}, 
                                               opts,
                                               t=0.0) 		
  evalElliptic(mesh, sbp, eqn, opts)

  if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "BDF2"
    BDF2(mesh, sbp, eqn, opts)
  end

  if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "CN"
    CN(mesh, sbp, eqn, opts)
  end
  if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "SDIRK4"
    SDIRK4(mesh, sbp, eqn, opts)
  end
end

"""
  time accurate advancement. Now it's totally a mess. Better to 
  separate BDF2 and CN, or figure out how to factor this funciton.
"""
function iterate{Tmsh, Tsol, Tres, Tdim, Tsbp}(mesh::AbstractMesh{Tmsh}, 
                                               pmesh::AbstractMesh{Tmsh},
                                               sbp::AbstractSBP{Tsbp}, 
                                               eqn::EllipticData{Tsol, Tres, Tdim}, 
                                               opts) 		
  dt = opts["delta_t"]
  ntimestep = opts["itermax"] 

  if haskey(opts, "write_energy") && opts["write_energy"] == true
    f = open("energy.dat", "w")

    eqn.calc_energy(mesh, sbp, eqn, opts, eqn.energy)
    println(f, 0, ", ", real(eqn.energy[1]))
  end

  if haskey(opts, "exact_soln_func")
    ferr = open("unsteady_error.dat", "w")
  end

  if haskey(opts, "Functional")
    ffunc = open("unsteady_functional.dat", "w")
  end
  #
  # initialize q1 and q2
  #
  if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "BDF2"
    for el = 1:mesh.numEl
      for n = 1:mesh.numNodesPerElement
        for dof = 1:mesh.numDofPerNode
          eqn.q1[dof, n, el] = real(eqn.q[dof, n, el])
          eqn.q2[dof, n, el] = eqn.q1[dof, n, el]
        end
      end
    end
  end

  if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "CN"
    #
    # evaluate spatial residual so that we have res1
    #
    evalElliptic(mesh, sbp, eqn, opts)
    for el = 1:mesh.numEl
      for n = 1:mesh.numNodesPerElement
        for dof = 1:mesh.numDofPerNode
          eqn.q1[dof, n, el] = real(eqn.q[dof, n, el])
          eqn.res1[dof, n, el] = real(eqn.res[dof, n, el])
        end
      end
    end
  end

  newton_data, jac, rhs_vec = setupNewton(mesh, pmesh, sbp, eqn, opts, evalRes, alloc_rhs=false)
  ctx_residual = (evalRes, )
  jac_func = myphysicsJac
  rhs_func = physicsRhs

  step_tol = 1.0e-6
  res_abstol = opts["res_abstol"]
  res_reltol = opts["res_reltol"]
  res_reltol0=opts["res_reltol0"]
  for t = 1:ntimestep
    #
    # Initialize solution at 1st stage
    #
    if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "SDIRK4"
      for el = 1:mesh.numEl
        for n = 1:mesh.numNodesPerElement
          for dof = 1:mesh.numDofPerNode
            eqn.q_irk[dof, n, el, 1] = eqn.q[dof, n, el]
          end
        end
      end
    end

    eqn.params.t += dt 

    eqn.istage = 1
    for s = 1:eqn.nstages
      # for CN, the spatial residual should be stored
      if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "CN"
        eqn.params.t -= dt 
        evalElliptic(mesh, sbp, eqn, opts)
        for i = 1 : length(eqn.res)
          eqn.res1[i] = real(eqn.res[i])
        end
        eqn.params.t += dt 
      end 
      newtonInner(newton_data, mesh, sbp, eqn, opts, 
                  rhs_func, jac_func, jac, rhs_vec, 
                  ctx_residual, eqn.params.t,
                  itermax=1, step_tol=step_tol, 
                  res_abstol=res_abstol, 
                  res_reltol=res_reltol, 
                  res_reltol0=res_reltol0)

      #
      # update old variables, i.e., q1, q2, res1
      #
      if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "BDF2"
        for el = 1:mesh.numEl
          for n = 1:mesh.numNodesPerElement
            for dof = 1:mesh.numDofPerNode
              eqn.q2[dof, n, el] = eqn.q1[dof, n, el]
              eqn.q1[dof, n, el] = real(eqn.q[dof, n, el])
            end
          end
        end
      end
      if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "CN"
        for el = 1:mesh.numEl
          for n = 1:mesh.numNodesPerElement
            for dof = 1:mesh.numDofPerNode
              eqn.q1[dof, n, el] = real(eqn.q[dof, n, el])
              # eqn.res1[dof, n, el] = real(eqn.res[dof, n, el])
            end
          end
        end
      end
      if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "SDIRK4"
        for el = 1:mesh.numEl
          for n = 1:mesh.numNodesPerElement
            for dof = 1:mesh.numDofPerNode
              eqn.q_irk[dof, n, el, s+1] = eqn.q[dof, n, el]
            end
          end
        end
      end

      eqn.istage = eqn.istage + 1

    end


    # postprocessing
    if haskey(opts, "write_energy") && opts["write_energy"] == true
      eqn.calc_energy(mesh, sbp, eqn, opts, eqn.energy)
      println(f, t, " ", real(eqn.energy[1]))
      flush(f)
      if real(eqn.energy[1]) < 1.e-12
        break
      elseif real(eqn.energy[1]) > 1.e12
        break
      end
    end

    if haskey(opts, "exact_soln_func")
      l2norm, lInfNorm = calcErrorNorm(mesh, sbp, eqn, opts)
      println(ferr, t, "    ", l2norm)
      flush(ferr)
    end

    if haskey(opts, "Functional")
      functional_value = zeros(Tsol, 3)
      eqn.functional(mesh, sbp, eqn, opts, functional_value)
      println(ffunc, t, " ", real(functional_value[1]))
      flush(ffunc)
    end
  end
  if haskey(opts, "write_energy") && opts["write_energy"] == true
    close(f)
  end
  if haskey(opts, "Functional")
    close(ffunc)
  end
  if haskey(opts, "exact_soln_func")
    close(ferr)
  end
  return nothing
end

"""
  This funciton is a adaptor/wrapper of physicsJac. We need it to void recompute the jacobian matrix for this linear problem.
"""
function myphysicsJac(newton_data::NewtonData, 
                      mesh, 
                      sbp, 
                      eqn, 
                      opts, 
                      jac, 
                      ctx_residual, 
                      t=0.0; 
                      is_preconditioned::Bool=false)
  if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "SDIRK4"
    func = evalElliptic
  else
    func = ctx_residual[1]
  end

  epsilon = opts["epsilon"]::Float64
  pert = complex(0, epsilon)

  if eqn.params.recompute_jac
    res_copy = copy(eqn.res)  # copy unperturbed residual
    # calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, func, res_copy, pert, jac, t)
    physicsJac(newton_data, mesh, sbp, eqn, opts, jac, ctx_residual, t)
    eqn.params.recompute_jac = false
  end

  if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] != "SDIRK4"
    return nothing
  end

  if haskey(opts, "TimeAdvance") && opts["TimeAdvance"] == "SDIRK4"
    for i = 1 : length(jac.nzval)
      jac.nzval[i] *= eqn.params.irk_c[eqn.istage, eqn.istage]
    end

    dt = opts["delta_t"]

    for el= 1:mesh.numEl
      jac = sview(mesh.jac, :, el)
      for n = 1:mesh.numNodesPerElement 
        for dof = 1:mesh.numDofPerNode
          col = (el-1) * mesh.numNodesPerElement * mesh.numDofPerNode + (n-1)*mesh.numDofPerNode + dof
          for i = jac.colptr[col] : jac.colptr[col+1] - 1
            row = jac.rowval[i]
            if row == col
              jac.nzval[i] += sbp.w[n]/jac[n]/dt
            end
          end
        end
      end
    end
  end

  return nothing
end
