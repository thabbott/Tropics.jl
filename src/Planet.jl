module Planet

struct Planet
	c_p
	g
	L_v
	R_a
	R_v
	e_sat0
	T_0
	σ
end

function e_sat(planet::Planet, T)
	return 	(
				planet.e_sat0 *
				exp(-planet.Lv / planet.Rv * (1.0/T - 1.0/planet.T_0))
			)
end

function q_sat(planet::Planet, T, p)
	return planet.R_a / planet.R_v * (e_sat(T) / p)
end

function invert_mse(planet::Planet, h, RH, z, p; T_ini = nothing)
	if T_ini == nothing
		T_ini = (h - planet.g * z) / planet.c_p
	end
	err = (T) -> 	(	
						planet.c_p * T + planet.g * z +
				  		planet.L_v * RH * q_sat(T, p_s) - h
				 	)
	T_low = T_ini - fac
	T_high = T_ini + fac
	while err(T_low)*err(T_high) >= 0 && T_low > 0.0
		fac = fac*2.0
		T_low = T_ini - fac
		T_high = T_ini + fac
	end
	if T_low < 0
		@printf('Error: could not find bracketing temperature interval')
	end
	T = Roots.find_zero(err, (T_low, T_high))
	q = q_sat(planet, T, p)
	return (T, p)
end

function moist_adiabatic_dTdz(planet::Planet, T, p)

	q = q_sat(planet, T, p)
	ϵ = planet.R_a / planet.R_v
	dry_adiabatic_dTdz = planet.g / planet.c_p
	fac_num = planet.L_v * q_sat / (planet.R_a * T)
	fac_den = planet.L_v^2.0 * q_sat * ϵ / (planet.c_p * planet.R_a * T^2.0)
	return dry_adiabatic_dTdz * (1.0 + fac_num) / (1.0 + fac_den)

end

end # module

