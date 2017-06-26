abstract AbstractFunctional

type volumeAverage <: AbstractFunctional
end
function call{Tmsh, Tsol, Tres, Tdim}(obj::volumeAverage,
                                      mesh::AbstractMesh{Tmsh},
                                      sbp::AbstractSBP,
                                      eqn::EllipticData{Tsol, Tres, Tdim},
                                      opts,
                                      val::Array{Tsol, 1})
  val[:] = 0.0
  nNodesPFace = mesh.numNodesPerFace
  nDofsPNode = mesh.numDofPerNode

  relax_coef = 1.0
  if haskey(opts, "unstable_coef")
    relax_coef = opts["unstable_coef"]
  end
  sbpface = mesh.sbpface
  wface = sview(sbpface.wface, :)
  w = sview(sbp.w, :)
  for elem = 1 : mesh.numEl
    jac = sview(mesh.jac, :, elem)
    q = sview(eqn.q, :, :, elem)
    for n = 1 : mesh.numNodesPerElement
      val[1] += w[n]/jac[n]*q[1, n]
    end
  end

  calcGradient(mesh, sbp, eqn, eqn.q, eqn.q_grad)
  interpolateBoundary(mesh, sbp, eqn, opts, eqn.q, eqn.q_bndry)
  # interpolateBoundary(mesh, sbp, eqn, opts, eqn.q_grad, eqn.q_grad_bndry)
  nrm = Array(Tmsh, mesh.numNodesPerFace, 2)
  area = Array(Tmsh, mesh.numNodesPerFace)
  Fv_face = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  for iBC = 1:mesh.numBC
    indx0 = mesh.bndry_offsets[iBC]
    indx1 = mesh.bndry_offsets[iBC+1] - 1

    bc_func = mesh.bndry_funcs[iBC]

    if isDirichlet(bc_func)
      for f = indx0:indx1
        bndry = mesh.bndryfaces[f]
        # geom info
        for n=1:mesh.numNodesPerFace
          dxidx = sview(mesh.dxidx_bndry, :, :, n, f)
          nrm_xi = sview(sbpface.normal, :, bndry.face)
          nrm[n,1] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
          nrm[n,2] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]
        end

        # viscous flux
        cmptFv_bndry(mesh, eqn, f, Fv_face)

        # contribution from Dirichlet term
        for n = 1:mesh.numNodesPerFace
          val[2] -= wface[n] * (nrm[n,1] * Fv_face[1, 1, n] + nrm[n,2] * Fv_face[2,1,n])
        end

        # we also have to consider the midification term
        q = sview(eqn.q_bndry, :, :, f)
        # penalty matrix
        pMat = zeros(Tsol, nDofsPNode, nNodesPFace, nNodesPFace)
        cmptBPMat(mesh, sbp, eqn, opts, f, pMat)
        gg = Array(Tsol, nDofsPNode)
        dq = Array(Tsol, nDofsPNode, nNodesPFace) 

        for n = 1:mesh.numNodesPerFace
          xy = sview(mesh.coords_bndry, :, n, f)
          # bc_func(xy, sview(gD, :, n))
          bc_func(xy, gg)
          # gD[:, n] = gg[:]
          for dof = 1:mesh.numDofPerNode
            dq[dof, n] = q[dof, n] - gg[dof]
          end
        end

        for n = 1 : mesh.numNodesPerFace
          val_j = 0.0
          for j = 1 : mesh.numNodesPerFace
            val_j += pMat[1, n, j] * dq[1, j] 
          end
          val[2] += val_j
        end

        val[2] *= relax_coef
      end

    elseif isNeumann(bc_func)
      for f = indx0:indx1
        bndry = mesh.bndryfaces[f]
        for n=1:mesh.numNodesPerFace
          dxidx = sview(mesh.dxidx_bndry, :, :, n, f)
          nrm_xi = sview(mesh.sbpface.normal, :, bndry.face)
          nrm[n,1] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
          nrm[n,2] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]
          area[n] = sqrt(nrm[n,1]*nrm[n,1] + nrm[n,2]*nrm[n,2])
        end
        q = sview(eqn.q_bndry, :, :, f)
        for n = 1:mesh.numNodesPerFace
          val[3] += wface[n]*area[n] * q[1, n] 
        end
      end
    else
      error("We should never get here!")
    end

  end

  return nothing
end



type volumeEnergy <: AbstractFunctional
end
function call{Tmsh, Tsol, Tres, Tdim}(obj::volumeEnergy,
                                      mesh::AbstractMesh{Tmsh},
                                      sbp::AbstractSBP,
                                      eqn::EllipticData{Tsol, Tres, Tdim},
                                      opts,
                                      val::Array{Tsol, 1})
  # @assert(length(val) == mesh.numDofPerNode)
  val[:] = 0.0

  w = sview(sbp.w, :)
  for elem = 1 : mesh.numEl
    jac = sview(mesh.jac, :, elem)
    q = sview(eqn.q, :, :, elem)
    for n = 1 : mesh.numNodesPerElement
      for dof = 1:mesh.numDofPerNode
        val[dof] += w[n]/jac[n]*q[dof, n]*q[dof,n]
      end
    end
  end

  return nothing
end

global const FunctionalDict = Dict{ASCIIString, AbstractFunctional}( 
  "volumeAverage" => volumeAverage(),
  "energy" => volumeEnergy()
)

function getFunctional(mesh::AbstractMesh,
                       sbp::AbstractSBP,
                       eqn::EllipticData,
                       opts)
  eqn.functional = FunctionalDict[opts["Functional"]]
end
