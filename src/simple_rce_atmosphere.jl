mutable struct SimpleRCEAtmosphere <: AbstractComponent

	# Prognostic variables
	h_m
	G_h_m

	# Diagnostic variables
	T_a
	q_a
	ρ_a
	h_b
	M_u
	p_lcl
	lw_down_s
	lw_up_toa
	Q_rad_lcl

	# Parameters
	RH_s
	p_s
	τ_0
	n_τ
	D_τ
	T_t
	sw_down_s
	ϵ_p

	# Constants
	planet::Thermodynamics.Planet

end

τ = τ_0 (p/p_s)^n 
p = p_s (τ/τ_0)^(1/n)

function surface_air!(atm::SimpleRCEAtmosphere)
	T_a = Thermodynamics.invert_mse(
		atm.planet, atm.h_b, atm.RH_s, 0.0, atm.p_s; T_ini = atm.T_a)
	atm.T_a = T_a
	atm.q_a = atm.RH_s * Thermodynamics.q_sat(atm.planet, T_a, atm.p_s)
	atm.ρ_a = atm.p_s / (atm.planet.R_a * atm.T_a)
end

function p_lcl!(atm::SimpleRCEAtmosphere)
	atm.p_lcl = Thermodynamics.approximate_p_lcl(
		atm.planet, atm.RH_s, atm.p_s)

function T_prof(atm::SimpleRCEAtmosphere, p_grid)
	Tgrid = zeros(size(pgrid))
	for ii in 1:length(pgrid)
		if pgrid[ii] <= p_lcl
			Tgrid[ii] = atm.T_a * (
				(pgrid[ii]/atm.p_s) ^ (atm.planet.R_a / atm.planet.c_pa)
			)
		elseif Tgrid[ii] <=


function rad_lw!(atm::SimpleRCEAtmosphere, lw_up_s; δτ = 0.01)

	# Generate optical depth, pressure, temperature profiles
	τ = τ_0:-δτ:0.0
	p = atm.p_s .* (τ./τ_0).^(1.0/atm.n_τ)
	T = zeros(size(p))
	T[1] = atm.T_a
	for ii in 2:length(p)
		if p[ii] <= p_lcl
			T[ii] = T[0] * (
				(p[ii]/p[0]) ^ (atm.planet.R_a / atm.planet.c_pa)
			)
		elseif T[ii-1] < atm.T_t
			dTdp = atm.R_a / atm.c_pa * Thermodynamics.γ(
				atm.planet, T[ii-1], p[ii-1])
			dp = p[ii] - p[ii-1]
			T[ii] = T[ii-1] + dTdp * dp
		else
			T[ii] = T_t
		end
	end

	# Calculate upward lw flux
	F_up = zeros(size(p))
	F_up[1] = lw_up_s
	for ii in 2:length(τ)
		dFdτ = D_τ * (F_up[ii-1] - atm.planet.σ * T[ii-1]^4.0)
		dτ = τ[ii] - τ[ii-1]
		F_up[ii] = F_up[ii-1] + dFdτ * dτ
	end

	# Calculate downward lw flux
	F_down = zeros(size(p))
	F_down[end] = 0.0
	for ii in (length(τ)-1):1
		dFdτ = -D_τ * (F_down[ii+1] - atm.planet.σ * T[ii+1]^4.0)
		dτ = τ[ii+1] - τ[ii]
		F_down[ii] = F_down[ii+1] + dFdτ * dτ
	end

	# Save longwave down at surface
	atm.lw_down_s = F_down[1]

	# Save upward longwave at TOA
	atm.lw_up_toa = F_up[end]

	# Save radiative cooling rate at LCL 
	ilcl = findall(p .> p_lcl)[1]
	F_net_above = F_up[ilcl] - F_down[ilcl]
	F_net_below = F_up[ilcl-1] - F_down[ilcl-1]
	dm = (p[ilcl-1] - p[ilcl])/atm.planet.g
	atm.Q_rad_lcl = (F_net_below - F_net_above) / dm

end


function M_u!(atm::SimpleRCEAtmosphere)
	T_lcl = (
		atm.T_a * (atm.p_lcl / atm.p_s) ^ 
		(atm.planet.R_a / atm.planet.c_p)
		)
	atm.M_u = (
		atm.Q_rad_lcl * 
		Planet.moist_adiabatic_dTdz(planet, T_lcl, atm.p_lcl) / atm.ϵ_p
		)
end

function h_b!(atm::SimpleRCEAtmosphere, shf, lhf)
	atm.h_b = atm.h_m + (shf + lhf) / atm.M_u
end

function lw_down_s!(atm::SimpleRCEAtmosphere)
	atm.lw_down_s = (
		atm.T_a * (atm.p_lw / atm.planet.p_s) ^ 
		(atm.planet.Ra / atm.planet.c_p)
		)
end

function update_state!(atm::SimpleRCEAtmosphere, shf, lhf)
	m_atm = (atm.p_s - atm.p_t) / atm.planet.g
	atm.G_h_m.set_tendency!((shf + lhf)/m_atm - atm.Q_rad)
	atm.h_m = atm.G_h_m.step!()
end