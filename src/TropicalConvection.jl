module TropicalConvection

import Roots
import Planet
using Printf

abstract type AbstractComponent end
abstract type AbstractCoupler end
abstract type AbstractTimestepper end

include("simple_surface.jl")
include("simple_rce_atmosphere.jl")
include("simple_rce_system.jl")

end # module