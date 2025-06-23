module fatscaling

using DataFrames
using Distributions
using LinearAlgebra
using Base.Threads

include("trait_and_rate_functions.jl")
include("dailyforage.jl")
include("smartpath.jl")

export 

#From trait and rate functions
find_velocity,
bite_size_allo,
number_of_chews,
chew_rate_allo,
chew_allo,
gut_capacity_g,
# maxfatstorage,
metabolic_cost,
indperarea,
fatrecoveryrate,
expectedlifetime,
foragingtime,
gut_turnover_time,
dailyfoodintake,
reactionwidth,
reactionheight,
carryingcapacity_costs,

dailyforage,
smartpath


end # module fatscaling
