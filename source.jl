type SrcExpTrigPoly0thDiffn <: SRCType
end
function call(obj::SrcExpTrigPoly0thDiffn, 
              xy::AbstractVector, 
              src::AbstractArray,
              t=0.0)
  k = 2.0
  ss = sin(2*k*pi*xy[1])*sin(2*k*pi*xy[2])
  cs = cos(2*k*pi*xy[1])*sin(2*k*pi*xy[2])
  sc = sin(2*k*pi*xy[1])*cos(2*k*pi*xy[2])
  ex = exp(xy[1] + xy[2])
  src[:] = (2. - 8*k*k*pi*pi) * ss + 4*k*pi*(cs + sc)
  src[:] *= -10.*ex
  return nothing
end

type SrcHicken2011 <: SRCType
end
function call(obj::SrcHicken2011, 
              xy::AbstractVector, 
              src::AbstractArray,
              t=0.0)
  xtmp = pi * (exp(xy[1]) - 1.0) / (e - 1.0) 
  xtmp_1 = pi * exp(xy[1]) / (e - 1.0) 
  xfac = sin(xtmp)
  yfac = exp(xy[2])
  xfac_x = xtmp_1 * cos(xtmp)
  xfac_xx = xfac_x - xtmp_1*xtmp_1*xfac
  yfac_y = yfac
  yfac_yy = yfac
  q = xfac * yfac
  q_x = xfac_x * yfac
  q_y = xfac * yfac_y
  q_xx = xfac_xx * yfac
  q_yy = xfac * yfac_yy
  lambda = pi * exp(xy[1]) / (e - 1) 
  lambda_x = lambda
  src[:] = lambda_x * q_x + lambda * q_xx + lambda * q_yy 
  src[:] = -src[:]
  return nothing
end
#
type SrcExpTrigPoly2ndDiffn <: SRCType
end
function call(obj::SrcExpTrigPoly2ndDiffn, 
              xy::AbstractVector, 
              src::AbstractArray,
              t=0.0)
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
  q_xx = ex * ((1 - 4*k*k*pi*pi) * ss + 4*k*pi * cs)
  q_xy = ex * (ss + 2*k*pi*cs + 2*k*pi*sc + 4*k*k*pi*pi*cc)
  q_yy = ex * ((1 - 4*k*k*pi*pi) * ss + 4*k*pi * sc)
  lambdaxx_x = 2*xy[1]
  lambdaxy_x = xy[2]
  lambdaxy_y = xy[1]
  lambdayy_y = 2*xy[2]
  src[:] = lambdaxx_x * q_x + lambdaxx * q_xx + 
           lambdaxy_x * q_y + lambdaxy * q_xy + 
           lambdaxy_y * q_x + lambdaxy * q_xy + 
           lambdayy_y * q_y + lambdayy * q_yy
  src[:] = -src[:]
  return nothing
end
#
#
# a scalar diffusion coefficient, 10
#
type SrcTrigPoly0thDiffn <: SRCType
end
function call(obj::SrcTrigPoly0thDiffn, 
              xy::AbstractVector, 
              src::AbstractArray,
              t=0.0)
  n = 2.0
  factor = 8.0*n*n*pi*pi
  ss = sin(2*n*pi*xy[1])*sin(2*n*pi*xy[2])
  src[:] = 0.0
  for i = 1:length(src)
    src[i] = factor*ss * 1.0e1
  end
  return nothing
end
#
# a 2nd order diffusion coefficient tensor
#
type SrcTrigPoly2ndDiffn <: SRCType
end
function call(obj::SrcTrigPoly2ndDiffn, 
              xy::AbstractVector, 
              src::AbstractArray,
              t=0.0)
  n = 2.0
  sc = sin(2*n*pi*xy[1])*cos(2*n*pi*xy[2])
  cs = cos(2*n*pi*xy[1])*sin(2*n*pi*xy[2])
  ss = sin(2*n*pi*xy[1])*sin(2*n*pi*xy[2])
  cc = cos(2*n*pi*xy[1])*cos(2*n*pi*xy[2])
  src[:] = 0.0
  for i = 1:length(src)
    src[i] += 6*n*pi*xy[1]*cs
    src[i] += 6*n*pi*xy[2]*sc
    src[i] += 4*n*n*pi*pi*(-2.0 - xy[1]*xy[1] - xy[2]*xy[2])*ss
    src[i] += 8*n*n*pi*pi*xy[1]*xy[2]*cc
    src[i] = -src[i]
  end
  return nothing
end

type SrcTrigPoly2ndDiffnUnsteady <: SRCType
end
function call(obj::SrcTrigPoly2ndDiffnUnsteady, 
              xy::AbstractVector, 
              src::AbstractArray,
              t=0.0)
  n = 2.0
  sc = sin(2*n*pi*xy[1])*cos(2*n*pi*xy[2])
  cs = cos(2*n*pi*xy[1])*sin(2*n*pi*xy[2])
  ss = sin(2*n*pi*xy[1])*sin(2*n*pi*xy[2])
  cc = cos(2*n*pi*xy[1])*cos(2*n*pi*xy[2])
  src[:] = 0.0
  c = -10.0
  for i = 1:length(src)
    src[i] += 6*n*pi*xy[1]*cs
    src[i] += 6*n*pi*xy[2]*sc
    src[i] += 4*n*n*pi*pi*(-2.0 - xy[1]*xy[1] - xy[2]*xy[2])*ss
    src[i] += 8*n*n*pi*pi*xy[1]*xy[2]*cc
    src[i] = exp(c*t)*(c*ss - src[i])
  end
  return nothing
end



#
# a 6th order diffusion coefficient tensor
#
type SrcTrigPoly6thDiffn <: SRCType
end
function call(obj::SrcTrigPoly6thDiffn, 
              xy::AbstractVector, 
              src::AbstractArray,
              t=0.0)
  n = 2.0
  sc = sin(2*n*pi*xy[1])*cos(2*n*pi*xy[2])
  cs = cos(2*n*pi*xy[1])*sin(2*n*pi*xy[2])
  ss = sin(2*n*pi*xy[1])*sin(2*n*pi*xy[2])
  cc = cos(2*n*pi*xy[1])*cos(2*n*pi*xy[2])
  src[:] = 0.0
  for i = 1:length(src)
    src[i] += 2*n*pi*(6*xy[1]^5 + 3*xy[1]^3*xy[2]^2)*cs
    src[i] += 2*n*pi*(6*xy[2]^5 + 3*xy[1]^2*xy[2]^3)*sc
    src[i] += 4*n*n*pi*pi*(-2.0 - xy[1]^6 - xy[2]^6)*ss
    src[i] += 8*n*n*pi*pi*xy[1]^3*xy[2]^3*cc
    src[i] = -src[i]
  end
  return nothing
end

global const SRCDict = Dict{ASCIIString, SRCType}(
 "SrcTrigPoly0thDiffn" => SrcTrigPoly0thDiffn(),
 "SrcExpTrigPoly0thDiffn" => SrcExpTrigPoly0thDiffn(),
 "SrcExpTrigPoly2ndDiffn" => SrcExpTrigPoly2ndDiffn(),
 "SrcTrigPoly2ndDiffnUnsteady" => SrcTrigPoly2ndDiffnUnsteady(),
 "SrcTrigPoly2ndDiffn" => SrcTrigPoly2ndDiffn(),
 "SrcTrigPoly6thDiffn" => SrcTrigPoly6thDiffn(),
 "SrcHicken2011" => SrcHicken2011(),
)


function getSrcFuntors(mesh::AbstractMesh,
                       sbp::AbstractSBP,
                       eqn::EllipticData,
                       opts)
  eqn.src_func = SRCDict[opts["SRC_name"]]
  return nothing
end

function calcSource{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                            sbp::AbstractSBP,
                                            eqn::EllipticData{Tsol, Tres, Tdim},
                                            srcFunc::SRCType,
                                            src::AbstractArray{Tres, 3})

  numElems = mesh.numEl
  numNodesPerElement = mesh.numNodesPerElement
  # xy = Array(Tmsh, Tdim)
  src_tmp = Array(eltype(src), size(src, 1))
  t = eqn.params.t
  for elem = 1:numElems
    # @bp
    for n = 1:numNodesPerElement
      xy = sview(mesh.coords, :, n, elem)
      srcFunc(xy, sview(src, :, n, elem), t)
    end
  end
end
