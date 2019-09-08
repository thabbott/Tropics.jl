module Tropics

import Roots
import Thermodynamics
using Printf

export 	tropical_ocean,
		tropical_atmosphere,
		diagnose!,
		update_state!

abstract type AbstractComponent end
abstract type AbstractCoupler end
abstract type AbstractTimestepper end

include("timesteppers.jl")
include("simple_surface.jl")
include("gray_rce_atmosphere.jl")
include("gray_rce_system.jl")

end # module