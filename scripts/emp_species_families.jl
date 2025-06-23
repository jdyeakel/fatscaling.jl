using Revise

using fatscaling

using DataFrames
using Distributions
using LinearAlgebra
using Base.Threads

using Plots
using StatsPlots
using UnicodePlots
using ColorSchemes
using JLD2
using ProgressMeter
using GLM
using LaTeXStrings
using CSV
using RCall

filename = smartpath("data/fat_full_dataset.csv")
fatdata = CSV.read(filename,DataFrame)


# fatdata = fatdata[findall(x->x!="mammal",fatdata.taxa),:] #excludes the Pond 85 data
# fatdata = fatdata[findall(x->x=="measured",fatdata.measure_estimate),:] #excludes extrapolated data


# GET FAMILY INFORMATION FOR EACH GENUS/SPECIES COMBOS

###############################################################################
# 1.  LOAD PACKAGES
###############################################################################
using DataFrames
using HTTP                    # standard library for web requests
using JSON3                   # fast JSON parser
using ProgressMeter           # optional but handy
using Formatting              # neat progress‐bar formatting (optional)

###############################################################################
# 2.  HELPER:  QUERY GBIF AND RETURN FAMILY
###############################################################################
"""
    gbif_family(species::AbstractString; verbose=false) -> Union{String,Missing}

Query GBIF's match endpoint for *species* (underscores or spaces; case-insensitive).
Returns the family name (String) or `missing` if none was found.
"""
function gbif_family(species::AbstractString; verbose=false)
    clean   = replace(species, "_" => " ")
    url     = "https://api.gbif.org/v1/species/match?name=$(HTTP.escapeuri(clean))"
    resp    = HTTP.get(url; readtimeout = 10)   #  timeout so the caller can handle errors
    if resp.status != 200
        verbose && @warn "GBIF returned status $(resp.status) for $species"
        return missing
    end
    data = JSON3.read(resp.body)
    return haskey(data, "family") ? String(data["family"]) : missing
end

###############################################################################
# 3.  BUILD A   species ➜ family   LOOK-UP TABLE WITH CACHING
###############################################################################
species  = unique(fatdata.taxa)          # ~183 distinct names in your CSV
famdict  = Dict{String,Union{String,Missing}}()

@showprogress 1 "Querying GBIF…" for sp in species
    famdict[sp] = gbif_family(sp)
end

###############################################################################
# 4.  MERGE FAMILY BACK INTO THE MAIN DATAFRAME
###############################################################################
fatdata.family = [famdict[row.taxa] for row in eachrow(fatdata)]

###############################################################################
# 5.  OPTIONAL:  WRITE THE NEW TABLE
###############################################################################
CSV.write(smartpath("data/fat_full_dataset_with_family.csv"), fatdata)
