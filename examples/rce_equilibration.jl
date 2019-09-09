using Tropics
using Printf
using PyPlot

# ==
# Set parameters for run
# ==
hour = 3600.
day = 24hour
Δt = 6hour
global t = 0.0
T = 500day 

# ==
# Create model
# ==
@printf("Creating model\n")
surf = tropical_ocean(Δt; T_s = 290.0, meters = 1.0)
atm = tropical_atmosphere(Δt; h_m = 310e3, τ_0 = 1.0)
system = Tropics.GrayRCESystem(surf, atm)

# ===
# Simple outputs
# ===
T_s = []
T_a = []
lw_up_toa = []
sw_net_toa = []
lw_down_s = []
lw_up_s = []
shf = []
lhf = []
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
	push!(lw_up_toa, system.atm.lw_up_toa)
	push!(sw_net_toa, system.atm.sw_down_toa - system.atm.sw_up_toa)
	push!(lw_down_s, system.atm.lw_down_s)
	push!(lw_up_s, system.atm.lw_up_s)
	push!(shf, system.surf.shf)
	push!(lhf, system.surf.lhf)
	push!(model_time, t / day)
	@printf(
		"Day %.2f: T_s = %.1f, imbalance @ surface = %.1f W m-2, @ toa = %.1f W m-2\n",
		t / day,
		system.surf.T_s,
		-(system.surf.shf + system.surf.lhf + 
		system.atm.lw_up_s - system.atm.lw_down_s +
		system.atm.sw_up_s - system.atm.sw_down_s),
		system.atm.sw_down_toa - system.atm.sw_up_toa - 
		system.atm.lw_up_toa)
end
@printf("Finished\n")

# ===
# Make plots
# ===
include("plot_rce_equilibration.jl")