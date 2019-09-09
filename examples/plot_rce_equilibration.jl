@printf("Plotting\n")

# Temperatures
figure()
plot(model_time, T_s, color = "black", label = "\$T_s\$")
plot(model_time, T_a, color = "blue", label = "\$T_a\$")
xlabel("Time (days)")
ylabel("Temperature (K)")
legend()
savefig("examples/rce_equilibration/T_s_T_a.pdf")

# TOA energy balance
figure()
plot(model_time, sw_net_toa, color = "blue", label = "\$F_{sw,toa}\$")
plot(model_time, lw_up_toa, color = "red", label = "\$F_{lw,toa}\$")
plot(model_time, sw_net_toa .- lw_up_toa, 
	color = "black", label = "imbalance")
xlabel("Time (days)")
ylabel("Flux (W m\$^{-2}\$)")
legend()
savefig("examples/rce_equilibration/toa_energy_balance.pdf")

# Surface energy balance
figure()
plot(model_time, sw_net_toa, color = "blue", label = "\$F_{sw,s}\$")
plot(model_time, lw_up_s .- lw_down_s, color = "red",
	label = "\$F_{lw,s}\$")
plot(model_time, shf, color = "magenta", label = "\$F_s\$")
plot(model_time, lhf, color = "green", label = "\$F_l\$")
plot(model_time, sw_net_toa .+ lw_down_s .- lw_up_s .- shf .- lhf,
	color = "black", label = "imbalance")
xlabel("Time (days)")
ylabel("Flux (W m\$^{-2}\$)")
legend()
savefig("examples/rce_equilibration/surface_energy_balance.pdf")