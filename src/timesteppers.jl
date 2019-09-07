mutable struct ForwardEulerTimeStepper <: AbstractTimestepper
	G
	Δt
end

function set_tendency!(ts::ForwardEulerTimeStepper, tendency)
	ts.G = tendency
end

function step(ts::ForwardEulerTimeStepper)
	return ts.G * ts.Δt
end