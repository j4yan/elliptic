#
# Calculate flux at edge cubature points using face-based form.
#
@debug function calcFaceViscousFlux{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},
													 sbp::AbstractSBP,
													 eqn::EllipticData{Tsol, Tres, Tdim},
													 opts,
													 interfaces::AbstractArray{Interface,1},
													 xflux_face::AbstractArray{Tres, 3},
													 yflux_face::AbstractArray{Tres, 3},
													 flux_face::AbstractArray{Tres, 3})
  # viscous flux jacobian
	Bvis = Array(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim)
	nfaces = length(interfaces)
	penalty_method = opts["Flux_name"]
	p = opts["order"]
	Cip = opts["Cip"]
	dq = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
	factor_shahbazi = Float64(p + 1.0)*Float64(p + Tdim)/(2.0*Tdim)
	sbpface = mesh.sbpface
	numFacesPerElem = 3
	nrm = Array(Tmsh, mesh.numNodesPerFace, Tdim)
	nrm0 = Array(Tmsh, mesh.numNodesPerFace, Tdim)
	area = Array(Tmsh, mesh.numNodesPerFace)
	Bvis_dqdxL = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
	Bvis_dqdxR = Array(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
	eigMaxL = Array(Tmsh, mesh.numDofPerNode)
	eigMaxR = Array(Tmsh, mesh.numDofPerNode)

	#
	# For scalar penalty, we use penalty;
	# for matrix penalty, we use Sat.
	#
	penalty      = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
	penalty_sat0 = Array(Tmsh, mesh.numDofPerNode, mesh.numNodesPerFace)
	penalty_sat  = Array(Tmsh, mesh.numDofPerNode, mesh.numNodesPerFace)
	penalty_shahbazi = Array(Tmsh, mesh.numDofPerNode, mesh.numNodesPerFace)

	Sat = Array(Tmsh, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerFace)
	#
	# |R 0| |Λxx Λxy| |R 0|^T
	# |0 R| |Λyx Λyy| |0 R|
	#
	sbpface = mesh.sbpface
	R = sview(sbpface.interp, :,:)
	RLR_L = Array(Tmsh, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerFace)
	RLR_R = Array(Tmsh, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerFace)
	RLR = Array(Tmsh, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerFace)
	BsqrtRHinvRtBsqrt = Array(Tmsh, mesh.numNodesPerFace, mesh.numNodesPerFace)
	HRBRH = Array(Tmsh, mesh.numNodesPerElement, mesh.numNodesPerElement)
	area_sum = sview(eqn.area_sum, :)
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

	println("rho_min = ", eigmin(BsqrtRHinvRtBsqrt), ", rho_max = ", sigma)

	for f = 1:nfaces    # loop over faces
		face = interfaces[f]
		eL = face.elementL
		eR = face.elementR
		fL = face.faceL
		fR = face.faceR
		permL = sview(sbpface.perm, :, fL)
		permR = sview(sbpface.perm, :, fR)
		stencilSize = sbpface.stencilsize

		#
		# Compute geometric info on face
		#
		for n=1:mesh.numNodesPerFace
			dxidx = sview(mesh.dxidx_face, :, :, n, f)
			# norm vector in reference element
			nrm_xi = sview(sbp.facenormal, :, fL)
			nrm[n,1] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
			nrm[n,2] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]

			area[n] = sqrt(nrm[n,1]*nrm[n,1] + nrm[n,2]*nrm[n,2])
			#
			# norm vector in physical domain without any scale, ie, |nrm|=1
			#
			nrm0[n,1] = nrm[n,1]/area[n]
			nrm0[n,2] = nrm[n,2]/area[n]
		end
		#
		# Compute the size of element and face, and then meas(face)/meas(elem)
		#
		elem_volL = 0.0
		elem_volR = 0.0
		face_area = 0.0
		for n = 1:mesh.numNodesPerElement
			elem_volL += sbp.w[n]/mesh.jac[n, eL]
			elem_volR += sbp.w[n]/mesh.jac[n, eR]
			# elem_volL += eqn.w[n, eL]
			# elem_volR += eqn.w[n, eR]
		end

		for n = 1:mesh.numNodesPerFace
			face_area += sbpface.wface[n]*area[n]
		end

		if penalty_method == "Shahbazi" || penalty_method == "SAT0"
			heL = elem_volL/eqn.area_sum[eL]
			heR = elem_volR/eqn.area_sum[eR]
			he = min(heL, heR)
		end
		area_weightL = area_sum[eL]/face_area
		area_weightR = area_sum[eR]/face_area

		#
		# face variables 
		#
		BvisL = sview(eqn.Bvis_face, :, :, :, 1, :, f) 
		BvisR = sview(eqn.Bvis_face, :, :, :, 2, :, f) 
		dqdxL = sview(eqn.q_grad_face, :, 1, :, f, 1)
		dqdyL = sview(eqn.q_grad_face, :, 1, :, f, 2)
		dqdxR = sview(eqn.q_grad_face, :, 2, :, f, 1)
		dqdyR = sview(eqn.q_grad_face, :, 2, :, f, 2)

		for n = 1:mesh.numNodesPerFace
			for var = 1:mesh.numDofPerNode
				Bvis_dqdxL[1, var, n] = BvisL[1, 1, var, n]*dqdxL[var, n] + BvisL[1, 2, var, n]*dqdyL[var, n] 
				Bvis_dqdxL[2, var, n] = BvisL[2, 1, var, n]*dqdxL[var, n] + BvisL[2, 2, var, n]*dqdyL[var, n] 
				Bvis_dqdxR[1, var, n] = BvisR[1, 1, var, n]*dqdxR[var, n] + BvisR[1, 2, var, n]*dqdyR[var, n] 
				Bvis_dqdxR[2, var, n] = BvisR[2, 1, var, n]*dqdxR[var, n] + BvisR[2, 2, var, n]*dqdyR[var, n] 
			end
		end

		#
		# First compute penalty
		#
		if penalty_method == "SAT0"  ### SATs, Cip = -1
			#
			# Compute SAT matrix (penalty parameter)
			# First compute the eigenvalue of 
			#
			#	|Λxx Λxy| |J 0|
			#	|Λyx Λyy| |0 J|
			#
			eigMaxL[:] = -1.0
			eigMaxR[:] = -1.0
			BvisL = sview(eqn.Bvis, :,:,:,:,eL)
			BvisR = sview(eqn.Bvis, :,:,:,:,eR)

			for dof = 1 : mesh.numDofPerNode
				for n = 1 : mesh.numNodesPerElement 
					# left element
					b = real(BvisL[1,1,dof,n] + BvisL[2,2,dof,n])
					ac = real(BvisL[1,1,dof,n] * BvisL[2,2,dof,n] - BvisL[1,2,dof,n] * BvisL[2,1,dof,n])
					root = 0.5*(b + sqrt(b*b - 4*ac))*mesh.jac[n, eL] 
					eigMaxL[dof] = max(eigMaxL[dof], root)

					# right element
					b = real(BvisR[1,1,dof,n] + BvisR[2,2,dof,n])
					ac = real(BvisR[1,1,dof,n] * BvisR[2,2,dof,n] - BvisR[1,2,dof,n] * BvisR[2,1,dof,n])
					root = 0.5*(b + sqrt(b*b - 4*ac))*mesh.jac[n, eR] 
					eigMaxR[dof] = max(eigMaxR[dof], root)
				end

				for n = 1:mesh.numNodesPerFace
					penalty[dof, n] = area_weightL*eigMaxL[dof] + area_weightR*eigMaxR[dof]
					# penalty[dof, n] = 2.0*he*max(eigMaxL[dof], eigMaxR[dof])
					# penalty[dof, n] *= 0.125*sigma
					penalty[dof, n] *= 0.25*sigma*area[n]*area[n]
				end
			end

		elseif  penalty_method == "Shahbazi"
			#
			# Reference:
			# Shahbazi, Mavriplis, Multigrid algorithms for high order discontinuous 
			# Galerkin discretization of the compressible Navier-Stokes equations, 
			# JCP, volume 228 Issue 21, November, 2009
			#
			Bvis_faceL = sview(eqn.Bvis_face, :,:, :, 1, :, f)
			Bvis_faceR = sview(eqn.Bvis_face, :,:, :, 2, :, f)

			#
			# Compute n^TΛn
			#
			for n = 1:mesh.numNodesPerFace
				for dof = 1:mesh.numDofPerNode
					penalty[dof, n] =  nrm0[n,1]*nrm0[n,1]*(Bvis_faceL[1, 1, dof, n] + Bvis_faceR[1,1,dof, n])
					penalty[dof, n] += nrm0[n,1]*nrm0[n,2]*(Bvis_faceL[1, 2, dof, n] + Bvis_faceR[1,2,dof, n])
					penalty[dof, n] += nrm0[n,2]*nrm0[n,1]*(Bvis_faceL[2, 1, dof, n] + Bvis_faceR[2,1,dof, n])
					penalty[dof, n] += nrm0[n,2]*nrm0[n,2]*(Bvis_faceL[2, 2, dof, n] + Bvis_faceR[2,2,dof, n])
					penalty[dof, n] *= 0.5*factor_shahbazi*area[n]/he
				end
			end

		elseif penalty_method == "SAT"
			#
			BvisL = sview(eqn.Bvis, :,:,:,:,eL)
			BvisR = sview(eqn.Bvis, :,:,:,:,eR)
			Sat[:,:,:] = 0.0	

			# RLR_L and RLR_R are reused for each (d1, d2)
			for d1 = 1:Tdim
				for d2 = 1:Tdim
					# Compute R Λ{d1,d2}/H R^T
					# RLR_L[:,:,:] = 0.0
					# RLR_L[:,:,:] = 0.0
					RLR[:,:,:] = 0.0
					for n1 = 1:mesh.numNodesPerFace
						for n2 = 1:mesh.numNodesPerFace
							for k = 1:stencilSize
								# RLR_L[:, n2, n1] += BvisL[d1, d2, :, permL[k]] / eqn.w[permL[k], eL] * R[k, n2] * R[k, n1] 
								# RLR_R[:, n2, n1] += BvisR[d1, d2, :, permR[k]] / eqn.w[permR[k], eR] * R[k, n2] * R[k, n1] 
								# RLR[:, n2, n1] += (BvisL[d1, d2, :, permL[k]] / eqn.w[permL[k], eL]+  BvisR[d1, d2, :, permL[k]] / eqn.w[permR[k], eL])* R[k, n2] * R[k, n1] 
								RLR[:, n2, n1] += (BvisL[d1, d2, :, permL[k]] / eqn.w[permL[k], eL]*area_weightL+  BvisR[d1, d2, :, permR[k]] / eqn.w[permR[k], eR] *area_weightR)* R[k, n2] * R[k, n1] 
							end
						end
					end

					# N_{d1}i} RLR N_{d2}
					for n1 = 1:mesh.numNodesPerFace
						for n2 = 1:mesh.numNodesPerFace
							# Sat[:, n2, n1] += nrm0[n2, d1]*(RLR_L[:, n2, n1] + RLR_R[:,n2, n1])*nrm0[n1, d2]
							Sat[:, n2, n1] += nrm0[n2, d1]*RLR[:, n2, n1]*nrm0[n1, d2]
						end
					end
				end
			end

			#
			# Left multiply by numFacesPerElem/4*B
			#
			for row = 1:mesh.numNodesPerFace
				for col = 1:mesh.numNodesPerFace
					# Sat[:, row, col] *= 0.25*numFacesPerElem*sbpface.wface[col]*area[col]
					Sat[:, row, col] *= 0.25*sbpface.wface[col]*area[col]
				end
			end
		else
			error("We should never get here")
			for n = 1:mesh.numNodesPerFace
				for dof = 1:mesh.numDofPerNode
					penalty[dof, n] = Cip*p*p/he
				end
			end

		end

		# RHR[:,:,:] = 0.0
		# for n1 = 1:mesh.numNodesPerFace
			# for n2 = 1:mesh.numNodesPerFace
				# for k = 1:stencilSize
					# # RLR_L[:, n2, n1] += BvisL[d1, d2, :, permL[k]] / eqn.w[permL[k], eL] * R[k, n2] * R[k, n1] 
					# # RLR_R[:, n2, n1] += BvisR[d1, d2, :, permR[k]] / eqn.w[permR[k], eR] * R[k, n2] * R[k, n1] 
					# RHR[n2, n1] +=  R[k, n2] * R[k, n1]/sbp.w[permL[k]] 
				# end
			# end
		# end
		# println(eigmax(RHR))

		# # @bp
		# for row = 1:mesh.numNodesPerFace
			# for col = 1:mesh.numNodesPerFace
				# RLR[:, row, col] *= 0.25*numFacesPerElem*sbpface.wface[row]*sbpface.wface[col]
			# end
		# end

		qL = sview(eqn.q_face, :, 1, :, f)
		qR = sview(eqn.q_face, :, 2, :, f)

		for n = 1:mesh.numNodesPerFace
			for dof = 1:mesh.numDofPerNode
				dq[dof, n] = qL[dof, n] - qR[dof, n]
			end
		end

		for n=1:mesh.numNodesPerFace
			# the other is `[[v]]⋅f(u+, u-) = (v+ - v-)g(u+, u-)⋅norm`.
			# vector f(u+, u-) is stored in xflux_face and yflux_face,
			# while scalar g(u+, u-) is stored in flux_face.

			# term-2 vector [[q]]
			xflux_face[:, n, f] = -0.5*dq[:,n]*nrm[n,1]
			yflux_face[:, n, f] = -0.5*dq[:,n]*nrm[n,2]


			# term-3, {{∇q}}⋅norm
			flux_face[:, n, f] = -0.5*((Bvis_dqdxL[1, :, n]+ Bvis_dqdxR[1, :, n])*nrm[n,1] + (Bvis_dqdxL[2, :, n]+ Bvis_dqdxR[2, :, n])*nrm[n,2])

		end
		#
		# term-4, penalty term 
		#            	

		#
		# area comes from dΓ = area*dΓ_{ξ}
		#
		if penalty_method == "Shahbazi" || penalty_method == "SAT0"|| penalty_method == "Hartman"  
			for n=1:mesh.numNodesPerFace
				for dof = 1:mesh.numDofPerNode
					flux_face[dof, n, f] += penalty[dof, n]*dq[dof, n]
				end
			end
		elseif penalty_method == "SAT" 
			for n=1:mesh.numNodesPerFace
				for n1 = 1:mesh.numNodesPerFace
					for dof = 1:mesh.numDofPerNode
						flux_face[dof, n, f] += Sat[dof, n, n1]*dq[dof, n1]*area[n]
					end
				end
			end
		else
			error("We should never get here")	
		end
	end # end of loop over all interfaces

	return nothing
end # end of function calcFaceFlux



function calcViscousFluxJac{Tmsh, Tsol, Tres}(params::ParamType{2, :conservative},
																							q::AbstractArray{Tsol, 1},
																							aux_vars::AbstractArray{Tres, 1},
																							Bvis::AbstractArray{Tsol, 4})
	@assert(size(Bvis, 1) == size(q, 1))
	@assert(size(Bvis, 2) == size(q, 1))
  press = g1etPressure(aux_vars)

  Bvis[1, 1, 1, 1] = 0.0 
  Bvis[1, 2, 1, 1] = 0.0 
  Bvis[1, 3, 1, 1] = 0.0 
  Bvis[1, 4, 1, 1] = 0.0 

	Bvis[2, 1, 1, 1] = -4.0/3.0*q[2]
  Bvis[2, 2, 1, 1] = 4.0/3.0 
  Bvis[2, 3, 1, 1] = 0.0 
  Bvis[2, 4, 1, 1] = 0.0 

	Bvis[3, 1, 1, 1] = -q[3]
  Bvis[3, 2, 1, 1] = 0.0 
  Bvis[3, 3, 1, 1] = 1.0 
  Bvis[3, 4, 1, 1] = 0.0 

	Bvis[4, 1, 1, 1] = -(4.0/3.0*q[2]*q[2] + q[3]*q[3] + )
  Bvis[4, 2, 1, 1] = 0.0 
  Bvis[4, 3, 1, 1] = 1.0 
  Bvis[4, 4, 1, 1] = 0.0 
end

function calcAuxVars{Tmsh, Tsol, Tres}(params::ParamType{2, :conservative},
																			 q::AbstractArray{Tsol, 1},
																			 aux_vars::AbstractArray{Tres, 1})
end
