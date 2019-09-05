mutable struct SimpleSurface <: AbstractComponent

	# Prognostic variables and timesteppers
	T_s
	G_T_s::AbstractTimestepper

	# Diagnostic variables
	lw_up_s
	sw_up_s
	shf
	lhf

	# Parameters
	C_s
	α_s
	r_a

	# Planet 
	planet::Planet

end

function lw_up_s!(surf::SimpleSurface)
	surf.lw_up_s = surf.planet.σ * surf.T_s^4.0
end

function sw_up_s!(surf::SimpleSurface, sw_down_s)
	surf.sw_up_s = sw_down_s * surf.α_s
end

function shf!(surf::SimpleSurface, ρ_a, T_a)
	surf.shf = ρ_a * surf.planet.c_p * (T_s - surf.T_a) / r_a
end

function lhf!(surf::SimpleSurface, ρ_a, q_a, p_s)
	q_sat_s = Planet.q_sat(surf.planet, surf.T_s, p_s)
	surf.lhf = ρ_a * surf.planet.L_v * (q_sat_s - q_a) / r_a
end

function update_state!(surf::SimpleSurface, lw_down_s, sw_down_s)

	surf.G_T_s.set_tendency!(
		(sw_down_s + lw_down_s - sw_up_s - lw_up_s - shf - lhf)
		/ surf.C_s)
	surf.T_s = surf.G_T_s.step!()

end