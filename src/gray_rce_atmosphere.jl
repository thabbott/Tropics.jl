mutable struct GrayRCEAtmosphere <: AbstractComponent

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
	lw_up_s
	sw_down_s
	sw_up_s
	sw_up_toa
	lw_up_toa
	Q_rad_lcl

	# Parameters
	RH_s
	p_s
	τ_0
	n_τ
	D_τ
	sw_down_toa
	ϵ_p

	# Constants
	planet::Thermodynamics.Planet

end

function tropical_atmosphere(Δt; h_m = 310e3, τ_0 = 1.0)
	planet = Thermodynamics.earth()
	G_h_m = ForwardEulerTimeStepper(0.0, Δt)
	RH_s = 0.8
	p_s = 1e5
	n_τ = 4.0
	D_τ = 1.5
	sw_down_toa = 350.
	ϵ_p = 0.25
	T_a = 300.0
	q_a = RH_s * Thermodynamics.q_sat(planet, T_a, p_s)
	ρ_a = p_s / (planet.R_a * T_a)
	h_b = Thermodynamics.mse(planet, T_a, 0.0, q_a)
	return GrayRCEAtmosphere(
		h_m,
		G_h_m,
		T_a, q_a, ρ_a,
		h_b, 0.0, 0.0, 
		0.0, 0.0, 0.0, 
		0.0, 0.0, 0.0, 
		0.0,
		RH_s,
		p_s,
		τ_0,
		n_τ,
		D_τ,
		sw_down_toa,
		ϵ_p,
		planet
	)
end

function surface_air!(atm::GrayRCEAtmosphere)
	T_a = Thermodynamics.invert_mse(
		atm.planet, atm.h_b, atm.RH_s, 0.0, atm.p_s; T_ini = atm.T_a)
	atm.T_a = T_a
	atm.q_a = atm.RH_s * Thermodynamics.q_sat(atm.planet, T_a, atm.p_s)
	atm.ρ_a = atm.p_s / (atm.planet.R_a * atm.T_a)
end

function p_lcl!(atm::GrayRCEAtmosphere)
	atm.p_lcl = 900e2
	# TODO: implement
	# atm.p_lcl = Thermodynamics.approximate_p_lcl(
	# 	atm.planet, atm.RH_s, atm.p_s)
end

function rad_lw!(atm::GrayRCEAtmosphere, T_s; δτ = 0.01)

	# Generate optical depth, pressure, temperature profiles
	τ = atm.τ_0:-δτ:0.0
	p = atm.p_s .* (τ./atm.τ_0).^(1.0/atm.n_τ)
	T = zeros(size(p))
	T[1] = atm.T_a
	for ii in 2:length(p)
		if p[ii] >= atm.p_lcl
			T[ii] = T[1] * (
				(p[ii]/p[1]) ^ (atm.planet.R_a / atm.planet.c_pa)
			)
		else
			dTdp = (
				(T[ii-1] / p[ii-1]) * 
				(atm.planet.R_a / atm.planet.c_pa) * 
				Thermodynamics.γ(atm.planet, T[ii-1], p[ii-1])
			)
			dp = p[ii] - p[ii-1]
			T[ii] = T[ii-1] + dTdp * dp
		end
	end

	# Calculate upward lw flux
	F_up = zeros(size(p))
	F_up[1] = atm.planet.σ * T_s^4.0
	for ii in 2:length(τ)
		dFdτ = atm.D_τ * (F_up[ii-1] - atm.planet.σ * T[ii-1]^4.0)
		dτ = τ[ii] - τ[ii-1]
		F_up[ii] = F_up[ii-1] + dFdτ * dτ
	end

	# Calculate downward lw flux
	F_down = zeros(size(p))
	F_down[end] = 0.0
	for ii in (length(τ)-1):-1:1
		dFdτ = atm.D_τ * (F_down[ii+1] - atm.planet.σ * T[ii+1]^4.0)
		dτ = τ[ii+1] - τ[ii]
		F_down[ii] = F_down[ii+1] + dFdτ * dτ
	end

	# Save longwave fluxes
	atm.lw_down_s = F_down[1]
	atm.lw_up_s = F_up[1]
	atm.lw_up_toa = F_up[end]

	# Save radiative cooling rate at LCL 
	ilcl = findall(p .< atm.p_lcl)[1]
	F_net_above = F_up[ilcl] - F_down[ilcl]
	F_net_below = F_up[ilcl-1] - F_down[ilcl-1]
	dm = (p[ilcl-1] - p[ilcl])/atm.planet.g
	atm.Q_rad_lcl = (F_net_below - F_net_above) / dm

end

function rad_sw!(atm::GrayRCEAtmosphere, α_s)
	atm.sw_down_s = atm.sw_down_toa
	atm.sw_up_s = α_s * atm.sw_down_s
	atm.sw_up_toa = atm.sw_up_s
end

function M_u!(atm::GrayRCEAtmosphere)
	# Calculate mass flux based on BLQE
	T_lcl = (
		atm.T_a * (atm.p_lcl / atm.p_s) ^ 
		(atm.planet.R_a / atm.planet.c_pa)
	)
	dTdz = (
		-atm.planet.g / atm.planet.c_pa *
		Thermodynamics.γ(atm.planet, T_lcl, atm.p_lcl)
	)
	S = atm.planet.c_pa * dTdz + atm.planet.g
	atm.M_u = -atm.Q_rad_lcl / (S * atm.ϵ_p)
end

function S!(atm::GrayRCEAtmosphere)
	# Calculate dry static stability just above LCL
end

function h_b!(atm::GrayRCEAtmosphere, shf, lhf)
	atm.h_b = atm.h_m + (shf + lhf) / atm.M_u
end

function update_state!(atm::GrayRCEAtmosphere, shf, lhf)
	m_atm = atm.p_s / atm.planet.g
	set_tendency!(atm.G_h_m, 
		(
			shf + lhf + 
			atm.lw_up_s - atm.lw_down_s - atm.lw_up_toa +
			atm.sw_down_toa - atm.sw_down_s + atm.sw_up_s - atm.sw_up_toa
		) / m_atm
	)
	atm.h_m = atm.h_m + step(atm.G_h_m)
end