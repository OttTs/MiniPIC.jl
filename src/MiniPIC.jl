module MiniPIC

using DataFrames: DataFrame
using ProgressMeter

const BOLTZMANN_CONSTANT = 1.380649e-23
const VACUUM_PERMITTIVITY = 8.8541878128E-12
const ELEMENTARY_CHARGE = 1.60217663E-19
const ELECTRON_MASS = 9.1093837E-31

include("fieldsolver.jl")
include("simulation.jl")

export BOLTZMANN_CONSTANT, VACUUM_PERMITTIVITY
export ELEMENTARY_CHARGE, ELECTRON_MASS
export simulate

end
