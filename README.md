# fatscaling.jl
## Code for manuscript: An origin for the allometric scaling of mammalian fat reflects populations living on the edge  
### Yeakel et al.  

*A framework for exploring how body-size scaling of **foraging gains, costs, and gut capacity** shapes fat storage across mammals.*. 

`fatscaling.jl` implements the core mathematics and Monte-Carlo routines behind Yeakel **et al.**, “Bioenergetic Trade-offs Determine Mass-dependent Fat Reserves” (in prep).  The package provides:  

### Core source files (`/src`)

| File | Purpose |
|------|---------|
| `trait_and_rate_functions.jl` | Mass-scaling relationships and helper utilities (velocity, bite size, chewing mechanics, gut capacity, metabolic rates, etc.). |
| `dailyforage.jl` | Monte-Carlo foraging-bout simulator that returns daily energy **gains**, **costs**, and bite counts. |
| `smartpath.jl` | Portable path helper that builds file locations relative to the package root, with optional index suffixes—keeps I/O reproducible across machines. |
| `fatscaling.jl` | Package entry point; loads the other source files and re-exports the public API. |

<!-- Scripts table will be added once analysis notebooks are finalized -->

### Analysis / figure-generation scripts (`/scripts`)

| Script | What it does |
|--------|--------------|
| `emp_species_families.jl` | Adds taxonomic context to the raw fat-mass CSV: queries the GBIF REST API for each *Genus_species* string, retrieves the matching family, and writes an updated dataset (`fat_full_dataset_with_family.csv`). |
| `empiricalfatscaling.jl` | Reproduces all **empirical** results: OLS & ANCOVA fits, mixed-effects and PGLS models, and the mass–fat scatter-plot figure. Also saves supplementary comparisons (author, habitat, size-class). |
| `fatscaling_analyses.jl` | Runs the **simulation** side of the paper. Monte-Carlo sweeps across body mass, resource richness (μ), patchiness (ζ), and four gut types to generate Figures 3–5 (+ supplementary). Writes large JLD2/CSV data blobs and publishes PDF figures in `rawfigures/`. |

### Mathematica notebooks (`/math`)

| Notebook | Purpose |
|----------|---------|
| `math/gain_costs_allometry_v2.nb` | Symbolic derivation of the allometric gain–cost functions used in the model and the effective slope γ (Eq. 3 in the manuscript). |
| `math/reserveratio.nb` | Analytical exploration of the reserve-ratio metric (ϕ̂ and ϕ^{dem}); contains parameter-sensitivity sweeps and verification plots. |


---

## Installation
```julia
pkg> add https://github.com/jdyeakel/fatscaling.jl
```

OR:  
```bash
# Clone and develop locally (recommended if you’re hacking on the source)
git clone https://github.com/jdyeakel/fatscaling.jl
cd fatscaling.jl

# Launch Julia and activate the project environment
julia --project -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

Prerequisites
	•	Julia ≥ 1.10
	•	Distributions, DataFrames, GLM, Plots, CSV, JLD2, ProgressMeter
These are declared in Project.toml; $ julia --project -e 'using Pkg; Pkg.instantiate()' will grab them.

⸻

```julia
###############################################################################
# Quick start – one-off demo run
###############################################################################
using fatscaling                 # loads the public API

# ----------------------- choose a single parameter set -----------------------
mass      = 35.0          # kg
gut       = "colon"       # one of "caecum", "colon", "non-rumen foregut", "rumen foregut"
teeth     = "all"
μ         = 1e-2          # resource richness (bites m⁻² d⁻¹)
ζ         = 1.0           # patchiness exponent
α         = 4.0           # Gamma shape parameter
edensity  = 18.2          # kJ g⁻¹ dry matter
###############################################################################

# --- mass-scaling traits -----------------------------------------------------
β            = bite_size_allo(mass)                 # bite size (g bite⁻¹)
chewrate     = chew_allo(mass, teeth)               # bites s⁻¹
tchew        = β / chewrate                         # s bite⁻¹
maxgut       = gut_capacity_g(mass, gut)            # g
bcost, fcost = metabolic_cost(mass)                 # kJ s⁻¹ (basal, foraging)
velocity     = find_velocity(mass)                  # m s⁻¹
tmax_bout, _ = foragingtime(mass) .* 3600           # s d⁻¹

# --- competitor load & bite-encounter distribution ---------------------------
_, n      = indperarea(mass)                        # conspecifics sharing patch
width     = reactionwidth(mass)
height    = reactionheight(mass)
na        = n / (width * height)                    # inds m⁻² in reaction plane

m_res     = μ / β                                   # bites m⁻²
m′        = m_res / na                              # competitor-scaled bite density
α′        = α * na^(ζ - 2)                          # patchiness-scaled shape
γ_dist    = Gamma(α′, m′ / α′)                      # Γ(k, θ) encounter distribution

# --- run the foraging simulation --------------------------------------------
gain, cost, n_bites = dailyforage(
    γ_dist, tchew, β, maxgut, velocity,
    bcost, fcost, edensity, tmax_bout
)

println("Daily energy gain  : $(round(gain,  2)) kJ")
println("Daily energetic cost: $(round(cost, 2)) kJ")
println("Daily net energy    : $(round(gain - cost, 2)) kJ")
println("Total bites taken   : $(n_bites)")
```



⸻

Citing

If you use fatscaling.jl in your work, please cite the forthcoming paper (preprint DOI pending).  


⸻

