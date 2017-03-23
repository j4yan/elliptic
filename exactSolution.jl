
abstract ExactSolutionType

export ExactSolutionType, ExactTrig, calcErrorL2Norm

type ExactTrig <: ExactSolutionType
end
function call{Tmsh, Tsol}(obj::ExactTrig, xy::AbstractArray{Tmsh}, q::AbstractArray{Tsol, 1})  # (numDofPerNode))
    k = 2.0
    q[:] = sin(2*k*pi*xy[1])*sin(2*k*pi*xy[2])
    return nothing
end

# type ExactPoly2nd <: ExactSolutionType
# end
# function call{Tmsh, Tsol}(obj::ExactPoly2nd, xy::AbstractArray{Tmsh}, q::AbstractArray{Tsol, 1})
	# q[:] =  -0.25*xy[1]*xy[1] - 0.25*xy[2]*xy[2]
	# return nothing
# end
type ExactPoly2nd <: ExactSolutionType
end
function call{Tmsh, Tsol}(obj::ExactPoly2nd, xy::AbstractArray{Tmsh}, q::AbstractArray{Tsol, 1})
	a = 1.0
	b = 1.0
	c = 1.0
	q[:] =  a*xy[1]*xy[1] + b*xy[1]*xy[2] + c*xy[2]*xy[2]
	return nothing
end


global const ExactDict = Dict{ASCIIString, ExactSolutionType}(
															  "ExactTrig" => ExactTrig(),
	"ExactPoly2nd" => ExactPoly2nd(),
	"ExactTrig" => ExactTrig()
)


@debug function calcErrorL2Norm{Tmsh, Tsol, Tres}(
        mesh::AbstractMesh{Tmsh},
        sbp::AbstractSBP,
        eqn::AbstractEllipticData{Tsol, Tres},
        opts)
    l2norm::Float64 = 0.
    lInfnorm::Float64 = 0.
    qe = Array(Tsol, mesh.numDofPerNode)
    exactFunc = ExactDict[opts["exactSolution"]]

	lInfnorm = -1.0
	elem_with_max_dq = 0
    for e = 1:mesh.numEl
        for n = 1:mesh.numNodesPerElement
            xy = sview(mesh.coords, :, n, e)
            exactFunc(xy, qe)
            q = sview(eqn.q, :, n, e)
			jac = mesh.jac[n, e]
            for v = 1:mesh.numDofPerNode
                # dq = Real(q[v] - qe[v])
                dq = Float64(q[v] - qe[v])
				l2norm += dq*dq*sbp.w[n]/jac
				# l2norm += dq*dq
				if lInfnorm < abs(dq)
					lInfnorm = abs(dq)
				end
			end
		end
    end
	
    return sqrt(l2norm), lInfnorm
end
