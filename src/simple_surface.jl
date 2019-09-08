mutable struct SimpleSurface <: AbstractComponent

	# Prognostic variables and timesteppers
	T_s
	G_T_s::AbstractTimestepper

	# Diagnostic variables
	shf
	lhf

	# Parameters
	C_s
	α_s
	r_a

	# Planet 
	planet::Thermodynamics.Planet

end

function tropical_ocean(Δt; T_s = 290.0, meters = 1.0)
	planet = Thermodynamics.earth()
	mwe = planet.c_pl * planet.ρ_l
	G_T_s = ForwardEulerTimeStepper(0.0, Δt)
	C_s = meters*mwe
	α_s = 0.2
	r_a = 100.
	return SimpleSurface(
		T_s,
		G_T_s,
		0.0,
		0.0,
		C_s,
		α_s,
		r_a,
		planet
	)
end

function shf!(surf::SimpleSurface, ρ_a, T_a)
	surf.shf = ρ_a * surf.planet.c_pa * (surf.T_s - T_a) / surf.r_a
end

function lhf!(surf::SimpleSurface, ρ_a, q_a, p_s)
	q_sat_s = Thermodynamics.q_sat(surf.planet, surf.T_s, p_s)
	surf.lhf = ρ_a * surf.planet.L_v * (q_sat_s - q_a) / surf.r_a
end

function update_state!(surf::SimpleSurface, lw_down_s, sw_down_s, lw_up_s, sw_up_s)

	set_tendency!(surf.G_T_s,
		(
			sw_down_s + lw_down_s - 
			sw_up_s - lw_up_s - 
			surf.shf - surf.lhf
		)
		/ surf.C_s)
	surf.T_s = surf.T_s + step(surf.G_T_s)

end