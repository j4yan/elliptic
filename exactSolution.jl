abstract ExactSolutionType

export ExactSolutionType, ExactTrig, calcErrorL2Norm

type ExactTrig <: ExactSolutionType
end
function call{Tmsh, Tsol}(obj::ExactTrig, 
                          xy::AbstractArray{Tmsh}, 
                          q::AbstractArray{Tsol, 1},
                          t=0.0)  # (numDofPerNode))
  k = 2.0
  q[:] = sin(2*k*pi*xy[1])*sin(2*k*pi*xy[2])
  return nothing
end

type ExactTrigUnsteady <: ExactSolutionType
end
function call{Tmsh, Tsol}(obj::ExactTrigUnsteady, 
                          xy::AbstractArray{Tmsh}, 
                          q::AbstractArray{Tsol, 1},
                          t=0.0)  # (numDofPerNode))
  k = 2.0
  c = -10
  q[:] = sin(2*k*pi*xy[1])*sin(2*k*pi*xy[2]) * exp(c*t)
  return nothing
end

type ExactTrigSlightUnsteady <: ExactSolutionType
end
function call{Tmsh, Tsol}(obj::ExactTrigSlightUnsteady, 
                          xy::AbstractArray{Tmsh}, 
                          q::AbstractArray{Tsol, 1},
                          t=0.0)  # (numDofPerNode))
  k = 2.0
  c = -1.0
  q[:] = sin(2*k*pi*xy[1])*sin(2*k*pi*xy[2]) * exp(c*t)
  return nothing
end


type ExactExpTrig <: ExactSolutionType
end
function call{Tmsh, Tsol}(obj::ExactExpTrig, 
                          xy::AbstractArray{Tmsh}, 
                          q::AbstractArray{Tsol, 1},
                          t=0.0)  # (numDofPerNode))
  k = 2.0
  ex = exp(xy[1] + xy[2])
  q[:] = ex * sin(2*k*pi*xy[1]) * sin(2*k*pi*xy[2])
  return nothing
end

type ExactPoly2nd <: ExactSolutionType
end
function call{Tmsh, Tsol}(obj::ExactPoly2nd, 
                          xy::AbstractArray{Tmsh}, 
                          q::AbstractArray{Tsol, 1},
                          t=0.0)
  a = 1.0
  b = 1.0
  c = 1.0
  q[:] =  a*xy[1]*xy[1] + b*xy[1]*xy[2] + c*xy[2]*xy[2]
  return nothing
end

type ExactHicken2011 <: ExactSolutionType
end
function call{Tmsh, Tsol}(obj::ExactHicken2011, 
                          xy::AbstractArray{Tmsh}, 
                          q::AbstractArray{Tsol, 1},
                          t=0.0)  # (numDofPerNode))
  ex = exp(xy[2])
  tmp = pi * (exp(xy[1]) -1.0) / (e - 1.0)
  q[:] = ex * sin(tmp)
  return nothing
end

global const ExactDict = Dict{ASCIIString, ExactSolutionType}(
 "ExactTrig" => ExactTrig(),
 "ExactTrigUnsteady" => ExactTrigUnsteady(),
 "ExactTrigSlightUnsteady" => ExactTrigSlightUnsteady(),
 "ExactPoly2nd" => ExactPoly2nd(),
 "ExactExpTrig" => ExactExpTrig(),
 "ExactHicken2011" => ExactHicken2011(),
)


function calcErrorL2Norm{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh},
                                           sbp::AbstractSBP,
                                           eqn::AbstractEllipticData{Tsol, Tres},
                                           opts)
  t = eqn.params.t
  l2norm::Float64 = 0.
  lInfnorm::Float64 = 0.
  qe = Array(Tsol, mesh.numDofPerNode)
  # exactFunc = ExactDict[opts["exactSolution"]]
  exactFunc = eqn.ExactFunc

  lInfnorm = -1.0
  elem_with_max_dq = 0
  for el = 1 : mesh.numEl
    for n = 1 : mesh.numNodesPerElement
      xy = sview(mesh.coords, :, n, el)
      exactFunc(xy, qe, t)
      q = sview(eqn.q, :, n, el)
      jac = mesh.jac[n, el]
      for v = 1:mesh.numDofPerNode
        dq = real(q[v] - qe[v])
        # dq = Float64(q[v] - qe[v])
        l2norm += dq*dq*sbp.w[n]/jac
        if lInfnorm < abs(dq)
          lInfnorm = abs(dq)
        end
      end
    end
  end

  return sqrt(l2norm), lInfnorm
end

