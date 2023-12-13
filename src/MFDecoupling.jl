module MFDecoupling

using DifferentialEquations
using LinearAlgebra, DelimitedFiles
using Logging: global_logger
using TerminalLoggers: TerminalLogger
using ProgressMeter
using JLD2
# using ForwardDiff

export read_inputs

global_logger(TerminalLogger());

# ForwardDiff.can_dual(::Type{ComplexF64}) = true

include("IO.jl")
include("rhs_old.jl")
include("rhs.jl")
include("rhs_complex.jl")

end # module MFDecoupling
