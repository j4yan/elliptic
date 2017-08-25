"""

  Functor type to compute diffusion coefficients.  

"""
abstract AbstractDiffn

"""
  
  a diffusion coefficent:
    | x^2+1   xy  |
    |  xy   y^2+1 |

"""
type DiffnPoly2nd<: AbstractDiffn
end
function call{Tmsh, Tsol}(obj::DiffnPoly2nd,
                          xy::AbstractArray{Tmsh},
                          lambda::AbstractArray{Tsol, 3})
  # @assert(size(lambda, 1) == Tdim)
  # @assert(size(lambda, 2) == Tdim)
  # the 3rd dimension should be dof per node	
	for dof = 1:size(lambda, 3)
		lambda[1, 1, dof] = xy[1]*xy[1] + 1.0
		lambda[1, 2, dof] = xy[1]*xy[2]
		lambda[2, 2, dof] = xy[2]*xy[2] + 1.0
		lambda[2, 1, dof] = lambda[1, 2, dof]
	end
	return nothing
end

"""

  a diffusion coefficient:
    |x^6+1   (xy)^3|
    |(xy)^3   y^6+1|

"""
type DiffnPoly6th<: AbstractDiffn
end
function call{Tmsh, Tsol}(obj::DiffnPoly6th,
						  xy::AbstractArray{Tmsh},
						  lambda::AbstractArray{Tsol, 3})
	# @assert(size(lambda, 1) == Tdim)
	# @assert(size(lambda, 2) == Tdim)
	# the 3rd dimension should be dof per node	
	for dof = 1:size(lambda, 3)
		lambda[1, 1, dof] = xy[1]^6 + 1.0
		lambda[1, 2, dof] = xy[1]^3*xy[2]^3
		lambda[2, 2, dof] = xy[2]^6 + 1.0
		lambda[2, 1, dof] = lambda[1, 2, dof]
	end
	return nothing
end

"""

  a constant diffusion coefficient

"""
type DiffnPoly0th<: AbstractDiffn
end

function call{Tmsh, Tsol}(obj::DiffnPoly0th,
						  xy::AbstractArray{Tmsh},
						  lambda::AbstractArray{Tsol, 3})
	# @assert(size(lambda, 1) == Tdim)
	# @assert(size(lambda, 2) == Tdim)
	# the 3rd dimension should be dof per node	
	for dof = 1:size(lambda, 3)
		lambda[1, 1, dof] = 1.0e1
		lambda[2, 2, dof] = 1.0e1 
		lambda[1, 2, dof] = 0.0 
		lambda[2, 1, dof] = 0.0
	end
	return nothing
end

"""

  Diffusion coefficient used by Hicken2011 siam jsc.
  This one is used since the adjoint bc is consitent with the objective.

"""
type DiffnHicken2011<: AbstractDiffn
end

function call{Tmsh, Tsol}(obj::DiffnHicken2011,
						  xy::AbstractArray{Tmsh},
						  lambda::AbstractArray{Tsol, 3})
	# @assert(size(lambda, 1) == Tdim)
	# @assert(size(lambda, 2) == Tdim)
	# the 3rd dimension should be dof per node	
	for dof = 1:size(lambda, 3)
    lambda[1, 1, dof] = pi * exp(xy[1]) / (e - 1.0)
    lambda[2, 2, dof] =lambda[1, 1, dof] 
		lambda[1, 2, dof] = 0.0 
		lambda[2, 1, dof] = 0.0
	end
	return nothing
end

"""
  
  Dictionary to look up a diffusion function. Everything a 
  new diffusion function is created, we need to add it 
  to this dictionary.

"""
global const DiffnDict = Dict{ASCIIString, AbstractDiffn}(
	"poly0th" => DiffnPoly0th(),
	"poly2nd" => DiffnPoly2nd(),
	"poly6th" => DiffnPoly6th(),
	"DiffnHicken2011" => DiffnHicken2011(),
)

"""

  get the function to compute diffusion coefficients
  
  **Input**
   * mesh
   * sbp
   * opts

  **Input/Output**
   * eqn

"""
function getDiffnFunc(mesh::AbstractMesh,
                      sbp::AbstractSBP,
                      eqn::EllipticData,
                      opts)
  eqn.diffusion_func = DiffnDict[opts["Diffusion"]]
end

function calcDiffn{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                           sbo::AbstractSBP,
                                           eqn::EllipticData{Tsol, Tres, Tdim},
                                           diffusion_func::AbstractDiffn,
                                           lambda::AbstractArray{Tres, 5} )

	numElems = mesh.numEl
	numNodesPerElement = mesh.numNodesPerElement

	for elem = 1:numElems
		# @bp
		for n = 1:numNodesPerElement
			xy = sview(mesh.coords, :, n, elem)

			diffusion_func(xy, sview(lambda, :, :, :, n, elem))
		end
	end
end

