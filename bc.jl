abstract AbstractDirichletBC <: BCType
abstract AbstractNeumannBC <: BCType
export isDirichlet, isNeumann

type DirichletAllZero <: AbstractDirichletBC
end
function call{Tmsh, Tsol}(obj::DirichletAllZero,
                          xy::AbstractArray{Tmsh},
                          gD::AbstractArray{Tsol, 1})  
  gD[:] = 0.0
  return nothing
end

type DirichletTrig <: AbstractDirichletBC
end
function call{Tmsh, Tsol}(obj::DirichletTrig, 
                          xy::AbstractArray{Tmsh}, 
                          q::AbstractArray{Tsol, 1})  
  k = 2.0
  q[:] = sin(2*k*pi*xy[1])*sin(2*k*pi*xy[2])
  if abs(q[1]) > 1.0E-10
    error("some problem!")
  end
  return nothing
end
type DirichletExpTrig <: AbstractDirichletBC
end
function call{Tmsh, Tsol}(obj::DirichletExpTrig, 
                          xy::AbstractArray{Tmsh}, 
                          q::AbstractArray{Tsol, 1}) 
  k = 2.0
  q[:] = sin(2*k*pi*xy[1])*sin(2*k*pi*xy[2]) * exp(xy[1] + xy[2])
  if abs(q[1]) > 1.0E-10
    error("some problem!")
  end
  return nothing
end


type DirichletPolynial2nd <: AbstractDirichletBC
end
function call{Tmsh, Tsol}(obj::DirichletPolynial2nd,
                          xy::AbstractArray{Tmsh},
                          gD::AbstractArray{Tsol, 1}) 
  a = 1.0
  b = 1.0
  c = 1.0
  gD[:] =  a*xy[1]*xy[1] + b*xy[1]*xy[2] + c*xy[2]*xy[2]
end

type DirichletHicken2011 <: AbstractDirichletBC
end
function call{Tmsh, Tsol}(obj::DirichletHicken2011,
                          xy::AbstractArray{Tmsh},
                          gD::AbstractArray{Tsol, 1}) 
  ex = exp(xy[2])
  tmp = pi * (exp(xy[1]) -1.0) / (e - 1.0)
  gD[:] = ex * sin(tmp)
end
type NeumannTrig <: AbstractNeumannBC
end
function call{Tmsh, Tsol}(obj::NeumannTrig,
                          xy::AbstractArray{Tmsh, 1},
                          nrm::AbstractArray{Tmsh, 1},
                          gN::AbstractArray{Tsol, 1})  
  k = 2.0
  # ss = sin(2*k*pi*xy[1])*sin(2*k*pi*xy[2])
  cs = cos(2*k*pi*xy[1])*sin(2*k*pi*xy[2])
  sc = sin(2*k*pi*xy[1])*cos(2*k*pi*xy[2])
  # cc = cos(2*k*pi*xy[1])*cos(2*k*pi*xy[2])
  lambdaxx = xy[1]*xy[1] + 1
  lambdaxy = xy[1]*xy[2]
  lambdayy = xy[2]*xy[2] + 1
  q_x = 2*k*pi * cs
  q_y = 2*k*pi * sc
  gN[:] = nrm[1] * (lambdaxx*q_x + lambdaxy*q_y) + 
  nrm[2] * (lambdaxy*q_x + lambdayy*q_y)
  return nothing
end
type NeumannExpTrig <: AbstractNeumannBC
end
function call{Tmsh, Tsol}(obj::NeumannExpTrig,
                          xy::AbstractArray{Tmsh, 1},
                          nrm::AbstractArray{Tmsh, 1},
                          gN::AbstractArray{Tsol, 1})  
  k = 2.0
  ss = sin(2*k*pi*xy[1])*sin(2*k*pi*xy[2])
  cs = cos(2*k*pi*xy[1])*sin(2*k*pi*xy[2])
  sc = sin(2*k*pi*xy[1])*cos(2*k*pi*xy[2])
  cc = cos(2*k*pi*xy[1])*cos(2*k*pi*xy[2])
  ex = exp(xy[1] + xy[2])
  lambdaxx = xy[1]*xy[1] + 1
  lambdaxy = xy[1]*xy[2]
  lambdayy = xy[2]*xy[2] + 1
  q_x = ex * (ss + 2*k*pi * cs)
  q_y = ex * (ss + 2*k*pi * sc)
  gN[:] = nrm[1] * (lambdaxx*q_x + lambdaxy*q_y) + 
  nrm[2] * (lambdaxy*q_x + lambdayy*q_y)
  return nothing
end

global const BCDict = Dict{ASCIIString, BCType}(
  "DirichletAllZero"     => DirichletAllZero(),
  "DirichletPolynial2nd" => DirichletPolynial2nd(),
  "DirichletTrig"        => DirichletTrig(),
  "DirichletExpTrig"     => DirichletExpTrig(),
  "DirichletHicken2011"  => DirichletHicken2011(),
  "NeumannTrig"          => NeumannTrig(),
  "NeumannExpTrig"       => NeumannExpTrig(),
)

"""

  defines what kind of a boundary condition is, Dirichlet, or Neumann.

  **Input**
   * dBC: an object of [`BCType`](@ref)

  **Output**
   true if the boundary condition is of type [`AbstractDirichletBC`]; false if it is not

"""

function isDirichlet(dBC::AbstractDirichletBC)
  return true
end
function isDirichlet(dBC::AbstractNeumannBC)
  return false
end

"""

  similar to [`isDirichlet`](@ref) 

  **Input**
   * dBC: an object of [`BCType`](@ref)

  **Output**
   true if the boundary condition is of type [`AbstractDirichletBC`]; false if it is not

"""

function isNeumann(dBC::AbstractDirichletBC)
  return false
end
function isNeumann(dBC::AbstractNeumannBC)
  return true
end

"""

  get the boundary functions into mesh.bndry_funcs.

  **Inputs**
   * sbp
   * opts
   * opts

  **Inputs/Ouptput**
   * mesh

"""

function getBCFunctors(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EllipticData, opts)
  for i = 1:mesh.numBC
    key = string("BC", i, "_name")
    val = opts[key]
    mesh.bndry_funcs[i] = BCDict[val]
  end
end

"""

  interpolate some scalar variable (usually the solution) from elements to interfaces.

  **Inputs**
   * mesh
   * sbp
   * eqn
   * opts
   * q: the variable to be interpolated

  **Input/Output**
   * q_bndry: the variable on the boundary

"""

function interpolateBoundary{Tsol, Tres}(mesh::AbstractDGMesh,
                                         sbp::AbstractSBP,
                                         eqn::AbstractEllipticData{Tsol, Tres},
                                         opts,
                                         q::AbstractArray{Tsol,3},
                                         q_bndry::AbstractArray{Tsol,3})
  boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, q, q_bndry)
end

function interpolateBoundary{Tsol, Tres}(mesh::AbstractDGMesh,
                                         sbp::AbstractSBP,
                                         eqn::AbstractEllipticData{Tsol, Tres},
                                         opts,
                                         grad::AbstractArray{Tsol,4},
                                         grad_bndry::AbstractArray{Tsol,4})
  # @assert(size(grad, 4) == size(grad_bndry, 4))   # Tdim
  # @assert(size(grad, 1) == size(grad_bndry, 1))   # numDofPerNode

  dim = size(grad, 4)
  for d=1:dim
    g = sview(grad, :, :, :, d)
    g_bndry = sview(grad_bndry, : , : , :, d)
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, g, g_bndry)
  end
end

"""

  compute the BR2 term and assemble into eqn.res
  TODO: Currently it only works for SBP-BR2 and SBP-SIPG. We may need to
        extend it to CDG (compact discontinuous galerkin). 

  **Input**
   * mesh
   * sbp
   * opts

  **Input/Output**
   * eqn
"""

function getBCFluxes{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                             sbp::AbstractSBP,
                                             eqn::EllipticData{Tsol, Tres, Tdim},
                                             opts)
  #
  # The boundary integrals are categorized into 3 classes.
  # The first is numerical and applies to all boundaries (both Dirichlet
  # and Neumann). The second class is `physical` Dirichlet and the third
  # one is `physical` Neumann. These 3 classes are dealed with separatedly
  # in 3 loops.
  #

  #
  # Take into account Direclet and Neumann boundary conditions
  #
  p = opts["order"]
  Cip = opts["Cip"]
  penalty_method = opts["Flux_name"]
  sbpface = mesh.sbpface
  dq = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)  
  sbpface = mesh.sbpface
  nrm = Array(Tmsh, mesh.numNodesPerFace, Tdim)
  nrm1 = Array(Tmsh, mesh.numNodesPerFace, Tdim)
  area = Array(Tmsh, mesh.numNodesPerFace)

  eigMax = Array(Tmsh, mesh.numDofPerNode)
  Fv_face = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)

  pMat = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerFace) 
  R = sview(sbpface.interp, :,:)
  RR = R.'*R
  area_sum = sview(eqn.area_sum, :)

  # perm = zeros(Tmsh, sbp.numnodes, sbpface.stencilsize)
  # Hinv = zeros(Tmsh, sbp.numnodes, sbp.numnodes)
  # Bsqrt = zeros(Tmsh, sbpface.numnodes, sbpface.numnodes)
  # for s = 1:sbpface.stencilsize
    # perm[sbpface.perm[s, 1], s] = 1.0
  # end
  # for i = 1:sbp.numnodes
    # Hinv[i,i] = 1.0/sbp.w[i]
  # end
  # for i = 1:sbpface.numnodes
    # Bsqrt[i,i] = sqrt(sbpface.wface[i])
  # end
  # BsqrtRHinvRtBsqrt = Bsqrt*R.'*perm.'*Hinv*perm*R*Bsqrt
  # sigma = eigmax(BsqrtRHinvRtBsqrt)

  relax_coef = 1.0
  if haskey(opts, "unstable_coef")
    relax_coef = opts["unstable_coef"]
  end

  for iBC = 1:mesh.numBC
    indx0 = mesh.bndry_offsets[iBC]
    indx1 = mesh.bndry_offsets[iBC+1] - 1
    gD = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
    gg = Array(Tsol, mesh.numDofPerNode)
    gN = Array(Tsol, mesh.numDofPerNode)

    bc_func = mesh.bndry_funcs[iBC]

    if isDirichlet(bc_func) # Dirichlet BC, term-5 and term-7
      for f = indx0:indx1
        xflux = sview(eqn.xflux_bndry, :,:,f)
        yflux = sview(eqn.yflux_bndry, :,:,f)
        flux  = sview(eqn.flux_bndry, :,:,f)
        bndry = mesh.bndryfaces[f]
        elem = bndry.element
        face = bndry.face
        perm = sview(sbpface.perm, :, face)
        stencilSize = sbpface.stencilsize

        # Compute geometric info on face
        for n=1:mesh.numNodesPerFace
          dxidx = sview(mesh.dxidx_bndry, :, :, n, f)
          nrm_xi = sview(mesh.sbpface.normal, :, bndry.face)
          nrm[n,1] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
          nrm[n,2] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]
          area[n] = sqrt(nrm[n,1]*nrm[n,1] + nrm[n,2]*nrm[n,2])
          nrm1[n,1] = nrm[n,1]/area[n]
          nrm1[n,2] = nrm[n,2]/area[n]
        end

        q = sview(eqn.q_bndry, :, :, f)
        for n = 1:mesh.numNodesPerFace
          xy = sview(mesh.coords_bndry, :, n, f)
          # bc_func(xy, sview(gD, :, n))
          bc_func(xy, gg)
          # gD[:, n] = gg[:]
          for dof = 1:mesh.numDofPerNode
            dq[dof, n] = q[dof, n] - gg[dof]
          end
        end

        cmptFv_bndry(mesh, eqn, f, Fv_face)

        # term-5, the real penalty term
        cmptBPMat(mesh, sbp, eqn, opts, f, pMat)
        if penalty_method == "SIPG" 
          for n = 1:mesh.numNodesPerFace
            for dof = 1 : mesh.numDofPerNode
              flux[dof, n] += pMat[dof, n, n]*dq[dof, n] * relax_coef
            end
          end
        elseif penalty_method == "BR2" 
          for n = 1:mesh.numNodesPerFace
            for n1 = 1:mesh.numNodesPerFace
              for dof = 1:mesh.numDofPerNode
                flux[dof, n] += pMat[dof, n, n1]*dq[dof, n1] * relax_coef
              end
            end
          end
        else
          error("We should never get here") 
        end

        for n = 1:mesh.numNodesPerFace
          for dof = 1 : mesh.numDofPerNode
            # term-2
            flux[dof, n] += -Fv_face[1, dof, n]*nrm[n,1] - Fv_face[2, dof, n]*nrm[n,2]
            # term-3 and 7
            xflux[dof, n] -= dq[dof,n] * nrm[n,1]
            yflux[dof, n] -= dq[dof,n] * nrm[n,2]
          end
        end
      end
    elseif isNeumann(bc_func)   # Neumann BC, term-8
      nrm_j = Array(Tmsh, 2)
      for f = indx0:indx1
        flux  = sview(eqn.flux_bndry, :,:,f)
        bndry = mesh.bndryfaces[f]
        # Compute geometric info on face
        for n=1:mesh.numNodesPerFace
          dxidx = sview(mesh.dxidx_bndry, :, :, n, f)
          nrm_xi = sview(mesh.sbpface.normal, :, bndry.face)
          nrm[n,1] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
          nrm[n,2] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]
          # area[n] = sqrt(nrm[n,1]*nrm[n,1] + nrm[n,2]*nrm[n,2])
        end

        for j = 1:mesh.numNodesPerFace
          xy = sview(mesh.coords_bndry, :, j, f)
          nrm_j[1] = nrm[j, 1]
          nrm_j[2] = nrm[j, 2]
          bc_func(xy, nrm_j, gN)
          # term-8
          flux[:, j] -= gN
        end
      end
    else
      error("We should never get here!")
    end
  end

  return nothing 
end

"""
  compute the boundary penalty matrix
  **Input**
   * mesh
   * sbp
   * eqn
   * opts
   * iface: the index of boundary face
  **Input/Output**
   * pMat: penalty matrix
"""

function cmptBPMat{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                                   sbp::AbstractSBP,
                                                   eqn::EllipticData{Tsol, Tres, Tdim},
                                                   opts,
                                                   iface::Int,
                                                   pMat::AbstractArray{Tsol, 3})
  bndry = mesh.bndryfaces[iface]
  face = bndry.face
  elem = bndry.element
  penalty_method = opts["Flux_name"]
  sbpface = mesh.sbpface
  stencilSize = sbpface.stencilsize
  perm = sview(sbpface.perm, :, face)
  wface = sview(sbpface.wface, :)
  R = sview(sbpface.interp, :,:)
  RHR = zeros(Tmsh, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerFace)
  nrm = Array(Tmsh, mesh.numNodesPerFace, Tdim)
  nrm1 = Array(Tmsh, mesh.numNodesPerFace, Tdim)
  area = Array(Tmsh, mesh.numNodesPerFace)
  cdelta = eqn.params.const_delta

  for n=1:mesh.numNodesPerFace
    dxidx = sview(mesh.dxidx_bndry, :, :, n, iface)
    nrm_xi = sview(mesh.sbpface.normal, :, face)
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
  
  area_sum = eqn.area_sum[elem]
  lambda_el = sview(eqn.lambda, :,:,:,:, elem)

  # On boundary, the area is weighted twice
  fweight = 0.5*area_sum/face_area

  if penalty_method == "SIPG"
    for dof = 1 : mesh.numDofPerNode
      eigMax = -1.0e10
      # left element
      for n = 1 : mesh.numNodesPerElement 
        b = real(lambda_el[1,1,dof,n]) + real(lambda_el[2,2,dof,n])
        ac = real(lambda_el[1,1,dof,n]) * real(lambda_el[2,2,dof,n]) - real(lambda_el[1,2,dof,n]) * real(lambda_el[2,1,dof,n])
        root = 0.5*(b + sqrt(b*b - 4*ac))*mesh.jac[n, elem]
        eigMax = max(eigMax, root)
      end

      for n = 1:mesh.numNodesPerFace
        # Ref. to Theorem 3 in the journal paper.
        # there are 2 area, one is from penalty parameter δ, 
        # the other is from B_\γ
        pMat[dof, n, n] = eigMax*cdelta*fweight*area[n]*area[n]
      end
    end
  elseif penalty_method == "BR2"
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
              RHR[:, n2, n1] += lambda_el[d1, d2, :, perm[k]] / eqn.w[perm[k], elem] * R[k, n2] * R[k, n1] 
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

    # Left multiply by 1/α*B
    for row = 1:mesh.numNodesPerFace
      for col = 1:mesh.numNodesPerFace
        pMat[:, row, col] *= fweight*wface[col]*area[col] * area[row]
      end
    end
  else
    error("We should never get here")
  end
end

"""
  compute the diffusion flux on a boundary face. There are possibly two ways: 
    1). interpolate the solution from element to boundary face, and then compute
        the flux.
    2). compute the elememt flux and then interpolate it to boundary face.
  The experiments show suble difference between two approaches. But the theoretical
  proof in the journal paper uses the second one.

  **Input**
   * mesh
   * eqn
   * iface: index of boundary face

  **Input/Output**
   * Fv_face: the diffusion flux

"""
function cmptFv_bndry{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                                              eqn::EllipticData{Tsol, Tres, Tdim},
                                              iface::Int, 
                                              Fv_face::AbstractArray{Tsol, 3})
  bndry = mesh.bndryfaces[iface]
  sbpface = mesh.sbpface
  elem = bndry.element
  face = bndry.face

  # prepareation for term-2, ie, Λ∇q. There are 2 ways to compute it:
  # 1) interpolate then compute
  # 2) compute then interpolate

  #
  # method-1
  #
  # dqdx = sview(eqn.q_grad_bndry, :, :, iface, 1)
  # dqdy = sview(eqn.q_grad_bndry, :, :, iface, 2)
  # lambda_face = sview(eqn.lambda_bndry, :,:,:,:, iface)
  # for n = 1:mesh.numNodesPerFace
    # for var = 1:mesh.numDofPerNode
      # Fv_face[1, var, n] = lambda_face[1, 1, var, n]*dqdx[var, n] + lambda_face[1, 2, var, n]*dqdy[var, n] 
      # Fv_face[2, var, n] = lambda_face[2, 1, var, n]*dqdx[var, n] + lambda_face[2, 2, var, n]*dqdy[var, n] 
    # end
  # end

  # return nothing

  #
  # method-2
  #
  # compute Λ∇q 
  dqdx = sview(eqn.q_grad, :, :, elem, 1)
  dqdy = sview(eqn.q_grad, :, :, elem, 2)
  lambda_el = sview(eqn.lambda, :,:,:,:,elem)
  Fv_elem = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerElement)

  for n = 1:mesh.numNodesPerElement
    for var = 1:mesh.numDofPerNode
      Fv_elem[1, var, n] = lambda_el[1, 1, var, n]*dqdx[var, n] + lambda_el[1, 2, var, n]*dqdy[var, n] 
      Fv_elem[2, var, n] = lambda_el[2, 1, var, n]*dqdx[var, n] + lambda_el[2, 2, var, n]*dqdy[var, n] 
    end
  end
  # interpolate Λ∇q to face 
  for dim = 1 : Tdim
    Fv_elem_d = slice(Fv_elem, dim, :, :)
    Fv_face_d = slice(Fv_face, dim, :, :)
    boundaryFaceInterpolate!(sbpface, face, Fv_elem_d, Fv_face_d)
  end

  return nothing
end
