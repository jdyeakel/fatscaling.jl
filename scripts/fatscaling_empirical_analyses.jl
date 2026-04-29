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

using CategoricalArrays, StatsModels, MixedModels


#--------------------------------------------------------

##################################
# OLS
##################################

filename = smartpath("data/fat_full_dataset.csv")
fatdata = CSV.read(filename,DataFrame)

# REMOVE POND 85 DATA (no species info)
fatdata = fatdata[findall(x->x!="mammal",fatdata.taxa),:] #excludes the Pond 85 data

# REMOVE EXTRAPOLATED DATA (isometric by definition)
fatdata = fatdata[findall(x->x=="measured",fatdata.measure_estimate),:] #excludes extrapolated data


##################################
# Fit OLS
##################################

model = lm(@formula(log(fatmass_kg) ~ log(mass_kg)), fatdata)
int = coef(model)[1]
slope = coef(model)[2]



##################################
# Fit OLS; Mass > 50 KG
##################################

model_large = lm(@formula(log(fatmass_kg) ~ log(mass_kg)), fatdata[fatdata.mass_kg .> 0.5, :])


fatdata.log_mass = log.(fatdata.mass_kg)      # once, up top
fatdata.sizecls  = categorical(ifelse.(fatdata.mass_kg .> 0.5, "large", "small"))

# 1.  pooled OLS (for reference)
ols_all = lm(@formula(log(fatmass_kg) ~ 1 + log(mass_kg)), fatdata)

ancova = lm(@formula(log(fatmass_kg) ~ 1 + log_mass +
                                     sizecls + sizecls & log_mass),
            fatdata)


# 3.  test whether the extra 3 terms matter
ft_size = ftest(ancova.model, ols_all.model)     # full vs. reduced
println("\nANCOVA omnibus test:")
display(ft_size)

# 4.  pull out class-specific slopes with SE and CI
β   = coef(ancova)
V   = vcov(ancova)
nms = coefnames(ancova)

function slope_ci(cls; α = 0.05)
    L = zeros(length(nms))
    L[nms .== "log_mass"] .= 1.0                     # correct baseline term
    int_term = "sizecls: $(cls) & log_mass"
    if any(nms .== int_term)                         # interaction may or may not exist
        L[nms .== int_term] .+= 1.0
    end
    b  = dot(L, β)
    se = sqrt(L' * V * L)
    tcrit = quantile(TDist(dof_residual(ancova)), 1 - α/2)
    ci = tcrit * se
    (slope = b, se = se, lower = b - ci, upper = b + ci)
end


sl_small = slope_ci("small")
sl_large = slope_ci("large")

println("\nSize-specific slopes:")
println("  small:  $(round(sl_small.slope; digits=3)) ± $(round(sl_small.se; digits=3))",
        "   (95 % CI $(round(sl_small.lower; digits=3))–$(round(sl_small.upper; digits=3)))")

println("  large:  $(round(sl_large.slope; digits=3)) ± $(round(sl_large.se; digits=3))",
        "   (95 % CI $(round(sl_large.lower; digits=3))–$(round(sl_large.upper; digits=3)))")



#-----------------------------------------------

##################################
# MASS-FAT EMPIRICAL FIGURE
##################################

filename = smartpath("data/fat_full_dataset.csv")
fatdata = CSV.read(filename,DataFrame)

# REMOVE POND 85 DATA (no species info)
fatdata = fatdata[findall(x->x!="mammal",fatdata.taxa),:] #excludes the Pond 85 data

# INCLUDE extrapolated data in the figure, but label appropriately
# fatdata = fatdata[findall(x->x=="measured",fatdata.measure_estimate),:] #excludes extrapolated data

# ID small/large species positions (0.5 kg cutoff)
largepos = findall(x->x>0.5,fatdata[!,:mass_kg])
smallpos = setdiff(collect(1:size(fatdata)[1]),largepos)
estpos = findall(x->x=="estimated",fatdata[!,:measure_estimate])



colors = palette(:tab10)
fatplot = scatter(fatdata[!,:mass_kg][largepos],fatdata[!,:fatmass_kg][largepos],
    xlabel="Body mass, M (kg)",
    ylabel="Fat mass (kg)",
    xscale=:log10,
    yscale=:log10,
    color=:black,
    label="M > 0.5 kg",
    frame=:box)
scatter!(fatplot,fatdata[!,:mass_kg][smallpos],fatdata[!,:fatmass_kg][smallpos],
    label="M < 0.5 kg",
    xscale=:log10,
    yscale=:log10,
    color=:black,
    alpha=0.5)
scatter!(fatplot,fatdata[!,:mass_kg][estpos],fatdata[!,:fatmass_kg][estpos],
    label="Estimated: Mf/M = 44%",
    xscale=:log10,
    yscale=:log10,
    color=:lightblue)
iexp = collect(-2:0.1:5);
massvec = [10^iexp[i] for i=1:length(iexp)]
plot!(fatplot,massvec,exp(coef(model)[1]).*massvec.^(coef(model)[2]),
    width=3,
    color=:orange,
    label = "Best fit: All",
    legend=:topleft,
    foreground_color_legend = nothing)
iexp = collect(-0.5:0.1:5);
massvec = [10^iexp[i] for i=1:length(iexp)]
plot!(fatplot,massvec,exp(coef(model_large)[1]).*massvec.^(coef(model_large)[2]),
    width=3,
    color=:red,
    label = "Best fit: M > 0.5 kg",
    legend=:topleft,
    foreground_color_legend = nothing)

figfile = smartpath("rawfigures/fig_fat.pdf")
Plots.savefig(fatplot, figfile)






#--------------------------------------------------------

##################################
# DATASET VS DATASET COMPARISONS
##################################


#ANCOVA COMPARISONS
filename = smartpath("data/fat_full_dataset.csv")
fatdata = CSV.read(filename,DataFrame)
fatdata = fatdata[findall(x->x!="mammal",fatdata.taxa),:] #excludes the Pond 85 data
fatdata = fatdata[findall(x->x=="measured",fatdata.measure_estimate),:] #excludes extrapolated data
# fatdata = fatdata[findall(x->x>0.5,fatdata[!,:mass_kg]),:]


using CSV, DataFrames, CategoricalArrays, StatsModels, GLM, MixedModels   # add CategoricalArrays

# log-transform once so we never repeat ourselves
fatdata.author   = categorical(fatdata.author)        # now a CategoricalVector
fatdata.log_mass = log.(fatdata.mass_kg)
fatdata.log_fat  = log.(fatdata.fatmass_kg)


# RANDOM INTERCEPT MIXED MODEL
m0 = fit(MixedModel, @formula(log_fat ~ 1 + log_mass + (1 | author)), fatdata) # random intercept *only*

# RANDOM INTERCEPT RANDOM SLOPES MIXED MODEL
mm = fit(MixedModel,
         @formula(log_fat ~ 1 + log_mass + (1 + log_mass | author)),
         fatdata)
print(mm)

#MODEL COMPARISON
MixedModels.likelihoodratiotest(mm, m0)

#Get Delta BIC by comparing BICs in mm and m0
deltabic = bic(mm) - bic(m0) #complex - simple


# FIXED EFFECTS ANCOVA
m_global  = lm(@formula(log_fat ~ 1 + log_mass), fatdata)

m_ancova  = lm(@formula(log_fat ~ 1 + log_mass + author + author & log_mass), fatdata)

ft = ftest(m_ancova.model, m_global.model)   # full, then reduced



##############################################
# TERRESTRIAL VS AQUATIC VS AERIAL COMPARISONS
##############################################

#
# Assess whether habitat category (terrestrial, aquatic, aerial) alters
# the allometric scaling of fat mass with body mass.
#


#ANCOVA COMPARISONS
filename = smartpath("data/fat_full_dataset.csv")
fatdata = CSV.read(filename,DataFrame)
fatdata = fatdata[findall(x->x!="mammal",fatdata.taxa),:] #excludes the Pond 85 data
fatdata = fatdata[findall(x->x=="measured",fatdata.measure_estimate),:] #excludes extrapolated data
# fatdata = fatdata[findall(x->x>0.5,fatdata[!,:mass_kg]),:]


#–– ensure habitat is treated as a categorical predictor ––#
fatdata.habitat = categorical(fatdata.habitat)
fatdata.log_mass = log.(fatdata.mass_kg)
fatdata.log_fat  = log.(fatdata.fatmass_kg)


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
# 1.  Mixed‑effects modelling
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

#
# Determine how many habitat levels are present
hab_lvls = levels(fatdata.habitat)
n_hab    = length(hab_lvls)

if n_hab > 1
    # (a) random intercept only
    m_h0 = fit(
        MixedModel,
        @formula(log_fat ~ 1 + log_mass + (1 | habitat)),
        fatdata;
        REML = false # turn off REML → use (full) maximum-likelihood
    )

    # (b) random intercept + random slope
    m_hm = fit(
        MixedModel,
        @formula(log_fat ~ 1 + log_mass + (1 + log_mass | habitat)),
        fatdata;
        REML = false # turn off REML → use (full) maximum-likelihood
    )

    println(m_hm)                                 # model summary
    MixedModels.likelihoodratiotest(m_hm, m_h0)   # LRT comparison
    deltabic_hab = bic(m_hm) - bic(m_h0)          # ΔBIC (complex – simple)
else
    @info "Only one habitat level detected ($(hab_lvls[1])); skipping mixed‑effects comparison."
end

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
# 2.  Fixed‑effects ANCOVA
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

m_hab_global = lm(
    @formula(log_fat ~ 1 + log_mass),
    fatdata
)

m_hab_ancova = lm(
    @formula(log_fat ~ 1 + log_mass + habitat + habitat & log_mass),
    fatdata
)

ft_hab = ftest(m_hab_ancova.model, m_hab_global.model)   # ANCOVA F‑test

β     = coef(m_hab_ancova)                # vector of estimates
V     = vcov(m_hab_ancova)                # covariance matrix
names = coefnames(m_hab_ancova)           # coefficient names

# helper: build the contrast vector for any habitat
function contrast_vec(hab::AbstractString, names::Vector{String})
    L = zeros(length(names))
    L[names .== "log_mass"] .= 1.0           # baseline slope
    int_term = "habitat: $(hab) & log_mass"  # interaction term
    L[names .== int_term] .+= 1.0
    return L
end

hab_lvls = levels(fatdata.habitat)
out = DataFrame(habitat = String[], slope = Float64[],
                SE = Float64[], lower = Float64[], upper = Float64[])

for h in hab_lvls
    L   = contrast_vec(h, names)
    b   = dot(L, β)                       # estimated slope
    se  = sqrt(L' * V * L)               # SE via delta method
    # ci  = 1.96 * se
    α = 0.05
    tcrit = quantile(TDist(dof_residual(m_hab_ancova)), 1-α/2)
    ci = tcrit * se
    push!(out, (h, b, se, b - ci, b + ci))
end

#Print habitat-specific OLS fits
out




##############################################
# PGLS ANALYSIS
##############################################

# PGLS Analysis
filename = smartpath("data/fat_full_dataset.csv")
fatdata = CSV.read(filename,DataFrame)
fatdata = fatdata[findall(x->x!="mammal",fatdata.taxa),:]
fatdata = fatdata[findall(x->x=="measured",fatdata.measure_estimate),:]
# fatdata = fatdata[findall(x->x>0.5,fatdata[!,:mass_kg]),:]

#Path to dated tree - place in R environ for later loading
tre_path = smartpath("data/dated_tree.tre")   # absolute path as String
@rput tre_path            # now R has an object called tre_path
bysp_clean_path = smartpath("data/bysp_clean.csv")
@rput bysp_clean_path            # now R has an object called bysp_clean_path


R"""
library(rotl); library(ape); library(caper); library(broom); library(phylolm)

# Import cleaned species list that aligns with the tree

# Import dated tree
tree<-read.tree(tre_path)
bysp_clean <- read.csv(bysp_clean_path, row.names = 1)

# replace the space between genus and species with an underscore to match EJ format
rownames(bysp_clean) <- gsub(" ", "_", rownames(bysp_clean))
tree$root.edge = 0



###############################################################################
##  MODELS  ###################################################################
###############################################################################
if (requireNamespace("phylolm", quietly = TRUE)) {
  fit_phy <- phylolm(log_fat ~ log_mass,
                     data = bysp_clean,
                     phy  = tree,
                     model = "BM", method = "ML")
  phy_name <- "PGLS (phylolm)"
} else {
  pic_mass <- pic(bysp_clean$log_mass, tree)
  pic_fat  <- pic(bysp_clean$log_fat,  tree)
  fit_phy  <- lm(pic_fat ~ pic_mass - 1)        # PIC fallback
  phy_name <- "PGLS (PIC)"
}

tidy_phylolm <- function(fit) {
  ct <- summary(fit)$coefficients
  pcol <- if ("p.value" %in% colnames(ct)) "p.value" else "Pr(>|t|)"
  data.frame(
    term      = rownames(ct),
    estimate  = ct[, "Estimate"],
    std.error = ct[, "StdErr"],
    statistic = ct[, "t.value"],
    p.value   = ct[, pcol],
    row.names = NULL
  )
}

#PAGEL'S LAMBDA version
fit_lam <- phylolm(
  log_fat ~ log_mass,
  data  = bysp_clean,
  phy   = tree,
  model = "lambda",
  REML  = FALSE   # optional; ML is the default anyway
)
tidy_lam <- tidy_phylolm(fit_lam)

lambda_hat <- fit_lam$optpar   # this is Pagel's lambda in the lambda model
lambda_hat

fit_ols <- lm(log_fat ~ log_mass, data = bysp_clean)

## --- OLS is always an lm object -------------------------------------------
tidy_ols <- tidy(fit_ols)


tidy_phy <- if (inherits(fit_phy, "phylolm")) tidy_phylolm(fit_phy) else tidy(fit_phy)


## --- combine and print -----------------------------------------------------
out <- rbind(
  cbind(model = "OLS",           tidy_ols),
  cbind(model = phy_name,        tidy_phy)
)
print(out)


#WITH ADDED > large-size class analysis
# 1. Filter to M > 0.5 kg

# PGLS AND OLS for all species (not just larger species)

large_species <- rownames(bysp_clean)[bysp_clean$log_mass > log(0.5)]
bysp_large <- bysp_clean[large_species, , drop = FALSE]

# 2. Prune the tree to match the large-bodied species
tree_large <- drop.tip(tree, setdiff(tree$tip.label, large_species))

# 3. Run models as before (OLS + PGLS)
fit_ols_large <- lm(log_fat ~ log_mass, data = bysp_large)

if (requireNamespace("phylolm", quietly = TRUE)) {
  fit_phy_large <- phylolm(log_fat ~ log_mass,
                           data = bysp_large,
                           phy  = tree_large,
                           model = "BM", method = "ML")
  phy_name_large <- "PGLS >0.5kg (phylolm)"
} else {
  pic_mass_large <- pic(bysp_large$log_mass, tree_large)
  pic_fat_large  <- pic(bysp_large$log_fat,  tree_large)
  fit_phy_large  <- lm(pic_fat_large ~ pic_mass_large - 1)
  phy_name_large <- "PGLS >0.5kg (PIC)"
}

# 4. Tidy and combine results
tidy_ols_large <- tidy(fit_ols_large)

if (inherits(fit_phy_large, "phylolm")) {
  ct   <- summary(fit_phy_large)$coefficients
  pcol <- if ("p.value" %in% colnames(ct)) "p.value" else "Pr(>|t|)"

  tidy_phy_large <- data.frame(
      term      = rownames(ct),
      estimate  = ct[, "Estimate"],
      std.error = ct[, "StdErr"],
      statistic = ct[, "t.value"],
      p.value   = ct[, pcol],
      row.names = NULL
  )
} else {
  tidy_phy_large <- tidy(fit_phy_large)
}

# 5. Append to output table
out_large <- rbind(
  cbind(model = "OLS >0.5kg", tidy_ols_large),
  cbind(model = phy_name_large, tidy_phy_large)
)

out_all <- rbind(out, out_large)  # original results already in `out`
print(out_all)
"""

@rget out_all      # now `out` is back in Julia as a DataFrame




#####################################
# CAPTIVE VS WILD - NAVARRETE DATASET
#####################################
using CSV, DataFrames, CategoricalArrays, StatsModels, GLM, MixedModels   # add 

filename = smartpath("data/Navarrete_fat.csv")
fatdata = CSV.read(filename,DataFrame)


#–– ensure habitat is treated as a categorical predictor ––#
fatdata.habitat = categorical(fatdata.habitat)
fatdata.log_mass = log.(fatdata.mass_kg)
fatdata.log_fat  = log.(fatdata.fatmass_kg)


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
# 1.  Mixed‑effects modelling
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

#
# Determine how many habitat levels are present
hab_lvls = levels(fatdata.habitat)
n_hab    = length(hab_lvls)

if n_hab > 1
    # (a) random intercept only
    m_h0 = fit(
        MixedModel,
        @formula(log_fat ~ 1 + log_mass + (1 | habitat)),
        fatdata;
        REML = false # turn off REML → use (full) maximum-likelihood
    )

    # (b) random intercept + random slope
    m_hm = fit(
        MixedModel,
        @formula(log_fat ~ 1 + log_mass + (1 + log_mass | habitat)),
        fatdata;
        REML = false # turn off REML → use (full) maximum-likelihood
    )

    println(m_hm)                                 # model summary
    MixedModels.likelihoodratiotest(m_hm, m_h0)   # LRT comparison
    deltabic_hab = bic(m_hm) - bic(m_h0)          # ΔBIC (complex – simple)
else
    @info "Only one habitat level detected ($(hab_lvls[1])); skipping mixed‑effects comparison."
end

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
# 2.  Fixed‑effects ANCOVA
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

m_hab_global = lm(
    @formula(log_fat ~ 1 + log_mass),
    fatdata
)

m_hab_ancova = lm(
    @formula(log_fat ~ 1 + log_mass + habitat + habitat & log_mass),
    fatdata
)

ft_hab = ftest(m_hab_ancova.model, m_hab_global.model)   # ANCOVA F‑test



using DataFrames, LinearAlgebra           # GLM already loaded

β     = coef(m_hab_ancova)                # vector of estimates
V     = vcov(m_hab_ancova)                # covariance matrix
names = coefnames(m_hab_ancova)           # coefficient names

# helper: build the contrast vector for any habitat
function contrast_vec(hab::AbstractString, names::Vector{String})
    L = zeros(length(names))
    L[names .== "log_mass"] .= 1.0           # baseline slope
    int_term = "habitat: $(hab) & log_mass"  # interaction term
    L[names .== int_term] .+= 1.0
    return L
end


hab_lvls = levels(fatdata.habitat)
out = DataFrame(habitat = String[], slope = Float64[],
                SE = Float64[], lower = Float64[], upper = Float64[])

for h in hab_lvls
    L   = contrast_vec(h, names)
    b   = dot(L, β)                       # estimated slope
    se  = sqrt(L' * V * L)               # SE via delta method
    # ci  = 1.96 * se
    α = 0.05
    tcrit = quantile(TDist(dof_residual(m_hab_ancova)), 1-α/2)
    ci = tcrit * se
    push!(out, (h, b, se, b - ci, b + ci))
end

out

