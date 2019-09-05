mutable struct SimpleRCESystem <: AbstractCoupler
	surface::SimpleSurface
	atmosphere::SimpleRCEAtmosphere
end

function diagnose!(sys::SimpleRCESystem)

	# Objective function for calculating boundary layer MSE 
	# (coupled to surface fluxes and convective mass flux)
	obj = function(h_b)
		sys.atm.h_b = h_b
		surface_air!(sys.atm)
		M_u!(sys.atm)
		shf!(sys.surf, sys.atm.ρ_a, sys.atm.T_a)
		lhf!(sys.surf, sys.atm.ρ_a, sys.atm.q_a, sys.atm.p_s)
		h_b!(sys.atm, sys.surf.shf, sys.surf.lhf)
		return sys.atm.h_b - h_b
	end

	# Find bracketing interval
	h_b_old = sys.atm.h_b
	fac = 1.0
	h_b_low = h_b_old - fac*sys.atm.planet.c_p
	h_b_high = h_b_old + fac*sys.atm.planet.c_p
	while obj(h_b_low)*obj(h_b_high) >= 0.0 && h_b_low > 0
		fac = 2.0*fac
		h_b_low = h_b_old - fac*sys.atm.planet.c_p
		h_b_high = h_b_old + fac*sys.atm.planet.c_p
	end
	if h_b_low < 0
		@printf('Error: could not find bracketing MSE interval')
	end

	# Find BL MSE
	sys.atm.h_b = find_root(obj, (h_b_low, h_b_high))

	# Calculate other diagnostic quantities
	surface_air!(sys.atm)
	M_u!(sys.atm)
	shf!(sys.surf, sys.atm.ρ_a, sys.atm.T_a)
	lhf!(sys.surf, sys.atm.ρ_a, sys.atm.q_a, sys.atm.p_s)
	lw_down_s!(sys.atm)
	lw_up_s!(sys.surf)
	sw_up_s!(sys.surf, sys.atm.sw_down_s)

end

function update_state!(sys::SimpleRCESystem)
	update_state!(sys.surf, sys.atm.lw_down_s, sys.atm.sw_down_s)
	update_state!(sys.atm, sys.surf.shf, sys.surf.lhf)
end