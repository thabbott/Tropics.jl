using Tropics
using Printf
using PyPlot

# ==
# Set parameters for run
# ==
hour = 3600.
day = 24hour
Δt = 1hour
global t = 0.0
T = 120day 

# ==
# Create model
# ==
@printf("Creating model\n")
surf = tropical_ocean(Δt; T_s = 290.0, meters = 0.1)
atm = tropical_atmosphere(Δt; h_m = 310e3, τ_0 = 1.0)
system = Tropics.GrayRCESystem(surf, atm)

# ===
# Simple outputs
# ===
T_s = []
T_a = []
lw_up_toa = []
sw_net_toa = []
model_time = []

# ==
# Run model
# ==
@printf("Running model\n")
while t < T
	diagnose!(system)
	update_state!(system)
	global t += Δt
	push!(T_s, system.surf.T_s)
	push!(T_a, system.atm.T_a)
	push!(model_time, t / day)
	@printf("Day %.2f: E_s = %.1f MJ m-2, E_atm = %.1f MJ m-2, E_tot = %.1f MJ m-2, imbalance = %.1f W m-2\n",
		t / day, system.surf.T_s * system.surf.C_s / 1e6, 
		system.atm.p_s / system.atm.planet.g * system.atm.h_m / 1e6,
		system.surf.T_s * system.surf.C_s / 1e6 +
		system.atm.p_s / system.atm.planet.g * system.atm.h_m / 1e6,
		system.atm.sw_down_toa - system.atm.sw_up_toa - system.atm.lw_up_toa)
end
@printf("Finished\n")

# ==
# Create figures
# ==
@printf("Plotting\n")
figure()
plot(model_time, T_s, color = "black", label = "\$T_s\$")
plot(model_time, T_a, color = "blue", label = "\$T_a\$")
xlabel("Time (days)")
ylabel("Temperature (K)")
legend()
savefig("T_s_T_a.pdf")