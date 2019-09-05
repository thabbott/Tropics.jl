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
	lw_down_s

	# Parameters
	Q_rad
	G_m
	RH_s
	p_lcl
	p_lw
	p_s
	p_t
	sw_down_s
	ϵ_p

	# Constants
	planet::Planet

end

function surface_air!(atm::SimpleRCEAtmosphere)
	(T_a, q_a) = Planet.invert_mse(
		atm.planet, atm.h_b, atm.RH_s, 0.0, atm.p_s; T_ini = atm.T_a)
	atm.T_a = T_a
	atm.q_a = q_a
	atm.ρ_a = atm.p_s / (atm.planet.R_a * atm.T_a)
end

function M_u!(atm::SimpleRCEAtmosphere)
	T_lcl = (
		atm.T_a * (atm.p_lcl / atm.p_s) ^ 
		(atm.planet.R_a / atm.planet.c_p)
		)
	atm.M_u = (
		atm.Q_rad * 
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