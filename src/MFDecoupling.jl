module MFDecoupling

using DifferentialEquations
using LinearAlgebra, DelimitedFiles
using Logging: global_logger
using TerminalLoggers: TerminalLogger
using ProgressMeter
using JLD2

export read_inputs

global_logger(TerminalLogger());

include("IO.jl")
include("rhs_old.jl")
include("rhs.jl")

end # module MFDecoupling
