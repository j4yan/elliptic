export interpolateFace
"""
  Given variables q at elem_volume nodes, compute corresponding gradients

  **Input**
   * mesh
   * sbp
   * eqn
   * q: the variable (usually the solution) the gradient of which we are computing

   **Input/Output**
    * q_grad: the gradient of q

"""
function calcGradient{Tmsh, Tsol, Tres, Tdim, Tsbp}(mesh::AbstractDGMesh{Tmsh},
                                                    sbp::AbstractSBP{Tsbp},
                                                    eqn::EllipticData{Tsol, Tres, Tdim},
                                                    q::AbstractArray{Tsol, 3},
                                                    q_grad::AbstractArray{Tsol, 4})
  # @assert(size(q, 3) == mesh.numEl)
  # @assert(size(q, 2) == mesh.numNodesPerElement)
  # @assert(size(q, 1) == mesh.numDofPerNode)

  # @assert(size(q_grad, 4) == Tdim)
  # @assert(size(q_grad, 3) == mesh.numEl)
  # @assert(size(q_grad, 2) == mesh.numNodesPerElement)
  # @assert(size(q_grad, 1) == mesh.numDofPerNode)

  numElems = mesh.numEl
  numNodes = mesh.numNodesPerElement
  numVars = mesh.numDofPerNode
  numNodesPerElement = mesh.numNodesPerElement
  dim = Tdim
  # dim = Tdim

  q_grad[:,:,:,:] = 0.0

  Dx = Array(Tsbp, numNodesPerElement, numNodesPerElement, 2)
  for e=1:numElems
    # First compute Dx for this element
    calcDx(mesh, sbp, Int(e), Dx)
    # q_grad = (Dxq, Dyq), Dx = H^(-1)Qx
    for n=1:numNodesPerElement
      for ivar=1:numVars
        for d=1:dim
          for col=1:numNodesPerElement
            q_grad[ivar, n, e, d] += Dx[n,col,d]*q[ivar, col, e]
          end
        end
      end
    end
  end
end


"""
similar to [`interpolateBoundary`](@ref); interpolate some scalar variable (usually the solution) from elements to interfaces. We may not need it any more if the evaluation of residual is rewritten as volume/face-based.

  **Inputs**
   * mesh
   * sbp
   * eqn
   * opts
   * q: the variable to be interpolated
   * q_face: the variable on the boundary

"""
function interpolateFace{Tsol}(mesh::AbstractDGMesh,
                               sbp,
                               eqn,
                               opts,
                               q::Abstract3DArray,
                               q_face::AbstractArray{Tsol, 4})
  # interpolate solution
  interiorfaceinterpolate!(mesh.sbpface, mesh.interfaces, q, q_face)
end

function interpolateFace{Tsol, Tres, Tdim}(mesh::AbstractDGMesh,
                                           sbp,
                                           eqn::EllipticData{Tsol, Tres, Tdim},
                                           opts,
                                           grad::AbstractArray{Tsol, 4},
                                           grad_face::AbstractArray{Tsol, 5})
  # interpolate gradient of solution
  for d=1: Tdim
    dqdx = sview(grad, :, :, :, d)
    dqdx_face = sview(grad_face, :, :, :, :, d)
    interiorfaceinterpolate!(mesh.sbpface, mesh.interfaces, dqdx, dqdx_face)
  end
end

"""
  Calculate flux at edge cubature points using face-based form. Only works for BR2 and SIPG.
  all the flux can be categorized into two groups, ∫ v f dΓ, and ∫ ∇⋅F dΓ

  **Input**
   * mesh
   * sbp
   * eqn
   * opts
   * interfaces: interfaces on which we are computing fluxes

  **Input/Output**
   * xflux_face: x component of `F`
   * yflux_face: y component of `F`
   * flux_face: the scalar flux `f`

"""
function calcFaceFlux{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},
                                                     sbp::AbstractSBP,
                                                     eqn::EllipticData{Tsol, Tres, Tdim},
                                                     opts,
                                                     interfaces::AbstractArray{Interface,1},
                                                     xflux_face::AbstractArray{Tres, 3},
                                                     yflux_face::AbstractArray{Tres, 3},
                                                     flux_face::AbstractArray{Tres, 3})

  nfaces = length(interfaces)
  penalty_method = opts["Flux_name"]
  p = opts["order"]
  Cip = opts["Cip"]
  dq = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  sbpface = mesh.sbpface
  nrm = Array(Tmsh, mesh.numNodesPerFace, Tdim)
  nrm1 = Array(Tmsh, mesh.numNodesPerFace, Tdim)
  area = Array(Tmsh, mesh.numNodesPerFace)
  FvL = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  FvR = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  Fv_eL = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerElement)
  Fv_eR = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerElement)

  pMat = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerFace)
  #
  # |R 0| |Λxx Λxy| |R 0|^T
  # |0 R| |Λyx Λyy| |0 R|
  #
  sbpface = mesh.sbpface
  R = sview(sbpface.interp, :,:)
  RHR = Array(Tmsh, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerFace)
  BsqrtRHinvRtBsqrt = Array(Tmsh, mesh.numNodesPerFace, mesh.numNodesPerFace)
  HRBRH = Array(Tmsh, mesh.numNodesPerElement, mesh.numNodesPerElement)

  perm = zeros(Tmsh, sbp.numnodes, sbpface.stencilsize)
  Hinv = zeros(Tmsh, sbp.numnodes, sbp.numnodes)
  Bsqrt = zeros(Tmsh, sbpface.numnodes, sbpface.numnodes)
  for s = 1:sbpface.stencilsize
    perm[sbpface.perm[s, 1], s] = 1.0
  end
  for i = 1:sbp.numnodes
    Hinv[i,i] = 1.0/sbp.w[i]
  end
  for i = 1:sbpface.numnodes
    Bsqrt[i,i] = sqrt(sbpface.wface[i])
  end

  BsqrtRHinvRtBsqrt = Bsqrt*R.'*perm.'*Hinv*perm*R*Bsqrt 
  sigma = eigmax(BsqrtRHinvRtBsqrt)

  sigma_err = abs(real(sigma - eqn.params.const_delta))
  if sigma_err > 1.0e-13
    println("there is sth wrong with const_delta ", sigma, ", ", eqn.params.const_delta)
  end
  relax_coef = 1.0
  if haskey(opts, "unstable_coef")
    relax_coef = opts["unstable_coef"]
  end
  for f = 1:nfaces    # loop over faces
    face = interfaces[f]
    eL = face.elementL
    eR = face.elementR
    fL = face.faceL
    fR = face.faceR
    stencilSize = sbpface.stencilsize

    #
    # Compute geometric info on face
    #
    # dxidx = sview(mesh.dxidx_face, :, :, :, f)
    for n=1:mesh.numNodesPerFace
      dxidx = sview(mesh.dxidx_face, :, :, n, f)
      # norm vector in reference element
      nrm_xi = sview(mesh.sbpface.normal, :, fL)
      nrm[n,1] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
      nrm[n,2] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]

      area[n] = sqrt(nrm[n,1]*nrm[n,1] + nrm[n,2]*nrm[n,2])
      #
      # norm vector in physical domain without any scale, ie, |nrm|=1
      #
      nrm1[n,1] = nrm[n,1]/area[n]
      nrm1[n,2] = nrm[n,2]/area[n]
      # println(nrm[n,1],", ", nrm[n,2], ", ", area[n])
    end

    # compute Fv 
    cmptFv_interface(mesh, eqn, f, FvL, FvR)

    # First compute penalty
    cmptIPMat(mesh, sbp, eqn, opts, f, pMat)
    
    # solution jump
    qL = sview(eqn.q_face, :, 1, :, f)
    qR = sview(eqn.q_face, :, 2, :, f)
    for n = 1:mesh.numNodesPerFace
      for dof = 1:mesh.numDofPerNode
        dq[dof, n] = qL[dof, n] - qR[dof, n]
      end
    end

    # finally we start to compute fluxes
    for n=1:mesh.numNodesPerFace
      # term-2 vector [[q]]
      xflux_face[:, n, f] = -0.5*dq[:,n]*nrm[n,1]
      yflux_face[:, n, f] = -0.5*dq[:,n]*nrm[n,2]

      # term-3, {{∇q}}⋅norm
      flux_face[:, n, f] = -0.5*((FvL[1, :, n]+ FvR[1, :, n])*nrm[n,1] + (FvL[2, :, n]+ FvR[2, :, n])*nrm[n,2])

    end
    
    # term-4, penalty term 
    if penalty_method == "SAT0"  
      for n=1:mesh.numNodesPerFace
        for dof = 1:mesh.numDofPerNode
          flux_face[dof, n, f] += pMat[dof, n, n]*dq[dof, n] * relax_coef
        end
      end
    elseif penalty_method == "SAT" 
      for n=1:mesh.numNodesPerFace
        for j = 1:mesh.numNodesPerFace
          for dof = 1:mesh.numDofPerNode
            flux_face[dof, n, f] += pMat[dof, n, j]*dq[dof, j] * relax_coef
          end
        end
      end
    else
      error("We should never get here")	
    end
  end # end of loop over all interfaces

  return nothing
end # end of function calcFaceFlux


"""

  Similar to [`cmptFv_bndry`](@ref);
  Compute the diffusion flux on an interface. There are possibly two ways: 
    1). interpolate the solution from element to interface, and then compute
        the flux.
    2). compute the elememt flux and then interpolate the flux instead of solution to interface.
  The experiments show suble difference between two approaches. But the theoretical
  proof in the journal paper uses the second one.

  **Input**
   * mesh
   * eqn
   * iface: index of boundary face

  **Input/Output**
   * Fv_face: the diffusion flux
"""
function cmptFv_interface{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                                                  eqn::EllipticData{Tsol, Tres, Tdim},
                                                  iface::Int, 
                                                  FvL::AbstractArray{Tsol, 3},
                                                  FvR::AbstractArray{Tsol, 3})
  interface = mesh.interfaces[iface]
  sbpface = mesh.sbpface
  eL = interface.elementL
  eR = interface.elementR
  fL = interface.faceL
  fR = interface.faceR

  # prepareation for term-2, ie, Λ∇q. There are 2 ways to compute it:
  # 1) interpolate then compute
  # 2) compute then interpolate

  #
  # method-1
  #
  #
  # face variables 
  #
  # lambda_fL = sview(eqn.lambda_face, :, :, :, 1, :, iface) 
  # lambda_fR = sview(eqn.lambda_face, :, :, :, 2, :, iface) 
  # dqdx_fL = sview(eqn.q_grad_face, :, 1, :, iface, 1)
  # dqdy_fL = sview(eqn.q_grad_face, :, 1, :, iface, 2)
  # dqdx_fR = sview(eqn.q_grad_face, :, 2, :, iface, 1)
  # dqdy_fR = sview(eqn.q_grad_face, :, 2, :, iface, 2)

  # for n = 1:mesh.numNodesPerFace
    # for var = 1:mesh.numDofPerNode
      # FvL[1, var, n] = lambda_fL[1, 1, var, n]*dqdx_fL[var, n] + lambda_fL[1, 2, var, n]*dqdy_fL[var, n] 
      # FvL[2, var, n] = lambda_fL[2, 1, var, n]*dqdx_fL[var, n] + lambda_fL[2, 2, var, n]*dqdy_fL[var, n] 
      # FvR[1, var, n] = lambda_fR[1, 1, var, n]*dqdx_fR[var, n] + lambda_fR[1, 2, var, n]*dqdy_fR[var, n] 
      # FvR[2, var, n] = lambda_fR[2, 1, var, n]*dqdx_fR[var, n] + lambda_fR[2, 2, var, n]*dqdy_fR[var, n] 
    # end
  # end
  # return nothing

  #
  # method-2
  #
  # compute Λ∇q 
  q_grad = sview(eqn.q_grad, :,:,:,:)
  lambda_eL = sview(eqn.lambda, :, :, :, :, eL)
  lambda_eR = sview(eqn.lambda, :, :, :, :, eR)
  Fv_eL = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerElement)
  Fv_eR = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerElement)

  for i = 1 : length(Fv_eL)
    Fv_eL[i] = 0.0
    Fv_eR[i] = 0.0
  end

  for d1 = 1:Tdim
    for d2 = 1:Tdim
      for node = 1:mesh.numNodesPerElement
        for ivar = 1:mesh.numDofPerNode
          Fv_eL[d1, ivar, node] += lambda_eL[d1, d2, ivar, node]*q_grad[ivar, node, eL, d2]  
          Fv_eR[d1, ivar, node] += lambda_eR[d1, d2, ivar, node]*q_grad[ivar, node, eR, d2]  
        end
      end
    end
  end

  # interpolate Λ∇q to face 
  for d = 1 : Tdim
    Fv_eL_d = slice(Fv_eL, d, :, :)
    Fv_eR_d = slice(Fv_eR, d, :, :)
    Fv_fL_d = slice(FvL, d, :, :)
    Fv_fR_d = slice(FvR, d, :, :)
    interiorFaceInterpolate!(sbpface, interface, Fv_eL_d, Fv_eR_d, Fv_fL_d, Fv_fR_d)
  end

  return nothing
end

"""
  compute the interior penalty matrix

  **Input**
   * mesh
   * sbp
   * eqn
   * opts
   * iface: the index of interface

  **Input/Output**
   * pMat: penalty matrix

"""
function cmptIPMat{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                           sbp::AbstractSBP,
                                           eqn::EllipticData{Tsol, Tres, Tdim},
                                           opts,
                                           iface::Int,
                                           pMat::AbstractArray{Tsol, 3})
  interface = mesh.interfaces[iface]
  eL = interface.elementL
  eR = interface.elementR
  fL = interface.faceL
  fR = interface.faceR
  sigma = eqn.params.const_delta
  numDofPerNode = mesh.numDofPerNode
  numNodesPerFace = mesh.numNodesPerFace
  penalty_method = opts["Flux_name"]
  sbpface = mesh.sbpface
  wface = sview(sbpface.wface, :)
  fL = interface.faceL
  fR = interface.faceR
  stencilSize = sbpface.stencilsize
  permL = sview(sbpface.perm, :, fL)
  permR = sview(sbpface.perm, :, fR)
  R = sview(sbpface.interp, :,:)
  RHR = zeros(Tmsh, numDofPerNode, numNodesPerFace, numNodesPerFace)
  nrm = Array(Tmsh, numNodesPerFace, Tdim)
  nrm1 = Array(Tmsh, numNodesPerFace, Tdim)
  area = Array(Tmsh, numNodesPerFace)

  for n=1:mesh.numNodesPerFace
    dxidx = sview(mesh.dxidx_face, :, :, n, iface)
    nrm_xi = sview(mesh.sbpface.normal, :, fL)
    nrm[n,1] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
    nrm[n,2] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]

    area[n] = sqrt(nrm[n,1]*nrm[n,1] + nrm[n,2]*nrm[n,2])
    nrm1[n,1] = nrm[n,1]/area[n]
    nrm1[n,2] = nrm[n,2]/area[n]
  end

  face_area = 0.0
  for n = 1:mesh.numNodesPerFace
    face_area += wface[n]*area[n]
  end
  
  area_sum = sview(eqn.area_sum, :)

  area_weightL = area_sum[eL]/face_area
  area_weightR = area_sum[eR]/face_area

  lambda_eL = sview(eqn.lambda, :,:,:,:,eL)
  lambda_eR = sview(eqn.lambda, :,:,:,:,eR)

  if penalty_method == "SAT0"  ### SATs, Cip = -1
    #
    # Compute SAT matrix (penalty parameter)
    # First compute the eigenvalue of 
    #
    #	|Λxx Λxy| |J 0|
    #	|Λyx Λyy| |0 J|
    #
    eigMaxL = Array(Tmsh, mesh.numDofPerNode)
    eigMaxR = Array(Tmsh, mesh.numDofPerNode)
    eigMaxL[:] = -1.0
    eigMaxR[:] = -1.0

    for dof = 1 : mesh.numDofPerNode
      for n = 1 : mesh.numNodesPerElement 
        # left element
        b = real(lambda_eL[1,1,dof,n] + lambda_eL[2,2,dof,n])
        ac = real(lambda_eL[1,1,dof,n] * lambda_eL[2,2,dof,n] - lambda_eL[1,2,dof,n] * lambda_eL[2,1,dof,n])
        root = 0.5*(b + sqrt(b*b - 4*ac))*mesh.jac[n, eL] 
        eigMaxL[dof] = max(eigMaxL[dof], root)

        # right element
        b = real(lambda_eR[1,1,dof,n] + lambda_eR[2,2,dof,n])
        ac = real(lambda_eR[1,1,dof,n] * lambda_eR[2,2,dof,n] - lambda_eR[1,2,dof,n] * lambda_eR[2,1,dof,n])
        root = 0.5*(b + sqrt(b*b - 4*ac))*mesh.jac[n, eR] 
        eigMaxR[dof] = max(eigMaxR[dof], root)
      end

      for n = 1:mesh.numNodesPerFace
        pMat[dof, n, n] = area_weightL*eigMaxL[dof] + area_weightR*eigMaxR[dof]
        pMat[dof, n, n] *= 0.25*sigma*area[n]*area[n]
      end
    end

  elseif penalty_method == "SAT"
    for i = 1 : length(pMat)
      pMat[i] = 0.0	
    end


    # Compute R Λ/H R^T
    for d1 = 1:Tdim
      for d2 = 1:Tdim
        for i = 1 : length(RHR)
          RHR[i] = 0.0
        end
        for n1 = 1:mesh.numNodesPerFace
          for n2 = 1:mesh.numNodesPerFace
            for k = 1:stencilSize
              RHR[:, n2, n1] += (lambda_eL[d1, d2, :, permL[k]] / eqn.w[permL[k], eL] * area_weightL 
                               + lambda_eR[d1, d2, :, permR[k]] / eqn.w[permR[k], eR] * area_weightR) * R[k, n2] * R[k, n1] 
            end
          end
        end

        # N RHR N
        for n1 = 1:mesh.numNodesPerFace
          for n2 = 1:mesh.numNodesPerFace
            pMat[:, n2, n1] += nrm1[n2, d1]*RHR[:, n2, n1]*nrm1[n1, d2]
          end
        end
      end
    end

    # Left multiply by numFacesPerElem/4*B
    for row = 1:mesh.numNodesPerFace
      for col = 1:mesh.numNodesPerFace
        pMat[:, row, col] *= 0.25*wface[col]*area[col]*area[row]
      end
    end
  else
    error("We should never get here")
  end

  return nothing
end

