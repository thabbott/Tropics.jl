mutable struct GrayRCESystem <: AbstractCoupler
	surf::SimpleSurface
	atm::GrayRCEAtmosphere
end

function diagnose!(sys::GrayRCESystem)

	# Objective function for calculating boundary layer MSE 
	# (coupled to surface fluxes and convective mass flux)
	obj = function(h_b)
		sys.atm.h_b = h_b
		surface_air!(sys.atm)
		p_lcl!(sys.atm)
		rad_lw!(sys.atm, sys.surf.T_s)
		M_u!(sys.atm)
		shf!(sys.surf, sys.atm.ρ_a, sys.atm.T_a)
		lhf!(sys.surf, sys.atm.ρ_a, sys.atm.q_a, sys.atm.p_s)
		h_b!(sys.atm, sys.surf.shf, sys.surf.lhf)
		return sys.atm.h_b - h_b
	end

	# Find bracketing interval
	h_b_old = sys.atm.h_b
	fac = 1.0
	h_b_low = h_b_old - fac*sys.atm.planet.c_pa
	h_b_high = h_b_old + fac*sys.atm.planet.c_pa
	while obj(h_b_low)*obj(h_b_high) >= 0.0 && h_b_low > 0
		fac = 2.0*fac
		h_b_low = h_b_old - fac*sys.atm.planet.c_pa
		h_b_high = h_b_old + fac*sys.atm.planet.c_pa
	end
	if h_b_low < 0
		throw(Thermodynamics.InversionError("could not find bracketing MSE interval"))
	end

	# Find BL MSE
	sys.atm.h_b = Roots.find_zero(obj, (h_b_low, h_b_high))

	# Calculate other diagnostic quantities
	surface_air!(sys.atm)
	p_lcl!(sys.atm)
	rad_lw!(sys.atm, sys.surf.T_s)
	rad_sw!(sys.atm, sys.surf.α_s)
	M_u!(sys.atm)
	shf!(sys.surf, sys.atm.ρ_a, sys.atm.T_a)
	lhf!(sys.surf, sys.atm.ρ_a, sys.atm.q_a, sys.atm.p_s)

end

function update_state!(sys::GrayRCESystem)
	update_state!(sys.surf, sys.atm.lw_down_s, sys.atm.sw_down_s,
		sys.atm.lw_up_s, sys.atm.sw_up_s)
	update_state!(sys.atm, sys.surf.shf, sys.surf.lhf)
end