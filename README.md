# fatscaling.jl

Repository for the manuscript **An origin for the allometric scaling of mammalian fat reflects populations living on the edge**.

Authors: Justin D. Yeakel, Mason Fan, Emily Jane McTavish, Anna Carolina de Almeida, Christopher P. Kempes, Sora Kim, Greg Breed

`fatscaling.jl` contains the data, Julia source code, Julia analysis scripts, and Wolfram Mathematica notebooks used to analyze empirical mammalian fat scaling and a bioenergetic foraging model. The repository is organized so that empirical data feed into `scripts/fatscaling_empirical_analyses.jl`, simulated model outputs feed into the Mathematica notebooks, and selected Mathematica analytical outputs feed back into `scripts/fatscaling_model_analyses.jl`. This README was assembled with help from an AI, with multiple human and silicon audits and sandboxed code checks.  

## Repository Organization

| Path | Description |
|---|---|
| `README.md` | This file. Describes repository contents, data columns, software requirements, and the analysis workflow. |
| `Project.toml` | Julia project file for the local `fatscaling` package. It declares the package name, UUID, author, version, direct package dependencies, and compatibility bounds for the core source code. |
| `Manifest.toml` | Julia manifest generated for the project environment. It records the exact dependency graph for the packages declared in `Project.toml`; it is machine-generated and should not be edited directly. |
| `.gitignore` | Git ignore rules for local/generated files. |
| `src/` | Julia source code defining the reusable `fatscaling` module. |
| `scripts/` | Julia analysis scripts. The two main scripts are `fatscaling_empirical_analyses.jl` and `fatscaling_model_analyses.jl`; `scripts/supplemental/` contains sensitivity-analysis scripts. |
| `data/` | Empirical data, phylogenetic data, source-reference key, and a legacy copy of one model-output CSV. |
| `data/simdata/` | Simulation-derived CSV files used for Julia-Mathematica exchange. |
| `math/` | Wolfram Mathematica notebooks used for analytical derivations and figure generation. |
| `rawfigures/` | Output directory used by the Julia scripts when saving PDF figures. |

Hidden operating-system files such as `.DS_Store` and Git internals are not analytical inputs.

## Julia Source Files

| File | Description |
|---|---|
| `src/fatscaling.jl` | Defines the `fatscaling` module, loads the source files below, and exports the public functions used by the scripts. |
| `src/trait_and_rate_functions.jl` | Allometric helper functions for mammalian traits and rates, including velocity, bite size, chewing rate, gut capacity, metabolic costs, population density in a foraging area, fat recovery rate, expected lifetime, foraging time, gut turnover time, daily food intake, reaction width/height, and carrying-capacity demographic costs. |
| `src/dailyforage.jl` | Monte Carlo daily foraging-bout simulator. Given an encounter-rate distribution and consumer/resource parameters, it simulates travel, chewing, gut filling, energy gains, energy costs, and number of bites taken during one day. |
| `src/smartpath.jl` | Path helper that resolves filenames relative to the package root and can append index suffixes before the file extension. The analysis scripts use this to find data and output locations reproducibly. |

## Main Analysis Scripts

### `scripts/fatscaling_empirical_analyses.jl`

This script reproduces the empirical analyses of mammalian fat mass as a function of body mass.

Main operations:

1. Reads `data/fat_full_dataset.csv`.
2. Removes rows where `taxa == "mammal"` for analyses requiring species information.
3. Removes rows where `measure_estimate != "measured"` for fitted statistical models; estimated rows are retained only in the empirical figure.
4. Fits OLS regressions of `log(fatmass_kg)` on `log(mass_kg)` for all measured records and for records with `mass_kg > 0.5`.
5. Fits an ANCOVA comparing small and large mammals using a 0.5 kg cutoff.
6. Creates the mass-fat scatter plot and saves it to `rawfigures/fig_fat.pdf`.
7. Tests whether fat-mass scaling differs by source dataset using mixed-effects models and fixed-effect ANCOVA with `author`.
8. Tests whether scaling differs by `habitat` using mixed-effects models and fixed-effect ANCOVA.
9. Runs phylogenetic analyses through R via `RCall`, reading `data/bysp_clean.csv` and `data/dated_tree.tre`; R packages used in that block are `rotl`, `ape`, `caper`, `broom`, and `phylolm`.
10. Reads `data/Navarrete_fat.csv` and compares captive/wild fat scaling using mixed-effects models and ANCOVA.

Data inputs:

| Input | Use |
|---|---|
| `data/fat_full_dataset.csv` | Primary empirical body-mass/fat-mass dataset. |
| `data/bysp_clean.csv` | Species-level log-mass/log-fat table aligned with the tree for PGLS. |
| `data/dated_tree.tre` | Dated phylogenetic tree used in PGLS. |
| `data/Navarrete_fat.csv` | Captive/wild comparison dataset. |
| `data/reference_key.xlsx` | Source-reference key for the `author` identifiers in the empirical datasets. It documents sources but is not read directly by the script. |

### `scripts/fatscaling_model_analyses.jl`

This script runs the foraging model, organizes simulation outputs, exchanges data with the Mathematica notebooks, and generates model-based figures.

Main operations:

1. Defines four gut types: `caecum`, `colon`, `non-rumen foregut`, and `rumen foregut`.
2. Defines body masses as `10 .^ collect(1.5:0.1:4.4)`, resource richness values as `10 .^ collect(-3.5:0.02:-1)`, patchiness values `zeta = [1.0, 1.5, 2.0, 2.15]`, energy density `18.2`, Gamma shape parameter `alpha = 4`, and `reps = 1000`.
3. For each gut type, replicate, body mass, richness, and patchiness value, computes allometric traits from `src/trait_and_rate_functions.jl`, builds a Gamma encounter distribution, and calls `dailyforage`.
4. Stores daily gains, costs, encounter counts, gut-capacity energy, gain-minus-gut-capacity, and gain surplus (`gains - costs`).
5. Estimates gain-surplus scaling slopes across richness values.
6. Saves and reloads a large JLD2 object at `../external_data/gainsremainder_acrossgut.jld2`. The script comments note that this file is about 2 GB and may need a user-specific path.
7. Reads `data/simdata/expgainremainderslope.CSV`, which is generated by the Mathematica notebook `math/gain_costs_allometry.nb`, and uses it as an analytical expected-slope curve.
8. Writes `data/simdata/gaincostpoor_data.csv` and the `gaincostpoor_zeta*_data.csv` files for use in `math/gain_costs_allometry.nb`.
9. Computes reserve-ratio quantities and writes `data/simdata/g1c1coords_zeta.csv` for use in `math/reserveratio.nb`.
10. Saves PDF figures to `rawfigures/`, including `fig_gutprediction.pdf`, `fig_gainsremainder_slopev2.pdf`, `fig_kappa.pdf`, `fig_zetastats.pdf`, and `fig_reserveratio_highzeta.pdf`.

## Supplemental Scripts

| File | Description |
|---|---|
| `scripts/supplemental/fatscaling_qualitybiomass.jl` | Sensitivity analysis based on the main model script in which energetic gains are multiplied by a digestible-quality term `theta` that varies with richness. The script writes supplementary figure outputs. |
| `scripts/supplemental/fatscaling_foragetime.jl` | Sensitivity analysis based on the main model script in which the foraging-time allometry is modified through `foragingtime_alt(mass, error)`. The script writes supplementary figure outputs. |

## Mathematica Notebooks and Julia-Mathematica Data Flow

| Notebook | Role in this repository |
|---|---|
| `math/gain_costs_allometry.nb` | Mathematica notebook that receives simulation summaries written by `scripts/fatscaling_model_analyses.jl` as `data/simdata/gaincostpoor_data.csv` and `data/simdata/gaincostpoor_zeta*_data.csv`. It also generates `data/simdata/expgainremainderslope.CSV`, which is read back into the Julia model script as an analytical expected-slope curve. |
| `math/reserveratio.nb` | Mathematica notebook that receives `data/simdata/g1c1coords_zeta.csv` from `scripts/fatscaling_model_analyses.jl` for reserve-ratio figure generation. |
| `math/storage_diff_drift_approx.nb` | Mathematica notebook retained with the analytical materials. It explores drift-dominated vs. diffusion-dominated fat storage assumptions described in the Supplemental Materials, and is used to assemble Figure S7. |

The main bidirectional exchange is:

1. Julia model simulations produce gain/cost CSV files in `data/simdata/`.
2. `math/gain_costs_allometry.nb` uses those CSV files and generates `expgainremainderslope.CSV`.
3. `scripts/fatscaling_model_analyses.jl` reads `expgainremainderslope.CSV` and overlays the analytical expected-slope curve on model results.
4. `scripts/fatscaling_model_analyses.jl` writes `g1c1coords_zeta.csv`.
5. `math/reserveratio.nb` reads `g1c1coords_zeta.csv` for the reserve-ratio analysis/figure.

## Data Files and Column Descriptions

### `data/fat_full_dataset.csv`

Primary empirical dataset with 428 data records. The script excludes rows with `taxa == "mammal"` when species-level identity is required and excludes `measure_estimate == "estimated"` from fitted models.

| Column | Description |
|---|---|
| `mass_kg` | Body mass in kilograms. |
| `fatmass_kg` | Fat mass in kilograms. |
| `taxa` | Taxon identifier, usually a lower-case genus-species string with underscores. The value `mammal` marks rows without species-level taxon information. |
| `author` | Source identifier. Values correspond to `ref id` entries in `data/reference_key.xlsx`. |
| `measure_estimate` | Whether the fat-mass value is `measured` or `estimated`. In the empirical figure, the script labels estimated values as `Mf/M = 44%`; fitted models use measured records only. |
| `habitat` | Habitat category used in habitat comparisons. Values in the file include `terrestrial`, `aquatic`, `aerial`, and `unknown`. |

### `data/Navarrete_fat.csv`

Dataset used for the captive/wild comparison, with 100 data records.

| Column | Description |
|---|---|
| `mass_kg` | Body mass in kilograms. |
| `fatmass_kg` | Fat mass in kilograms. |
| `habitat` | Captivity/status category used by the script as a categorical predictor. Values in the file include `Wild`, `Captive`, and `No data`. |

### `data/bysp_clean.csv`

Species-level table used by the R/PGLS block in `fatscaling_empirical_analyses.jl`, with 158 species records.

| Column | Description |
|---|---|
| unnamed first column | Species name. The R code reads this as row names and replaces spaces with underscores to match tree tip labels. |
| `log_mass` | Natural log of species-level body mass in kilograms. |
| `log_fat` | Natural log of species-level fat mass in kilograms. |

### `data/dated_tree.tre`

Dated phylogenetic tree in Newick format. The empirical script reads it in R with `ape::read.tree()` and uses it with `bysp_clean.csv` for PGLS models.

### `data/reference_key.xlsx`

Reference key for empirical source identifiers. The file contains seven source records.

| Column | Description |
|---|---|
| `ref id` | Short source identifier used in `author` columns, for example `pond_1985`, `pitts_1968`, and `navarrete_2011`. |
| `reference` | Full bibliographic reference for the source identifier. |

### `data/gaincostpoor_data.csv`

Legacy copy of a 30-row simulation gain/cost table with the same column structure and values as `data/simdata/gaincostpoor_data.csv`. The current model script writes and exchanges the version under `data/simdata/`.

| Column | Description |
|---|---|
| `mass` | Body mass in kilograms. |
| `gains` | Mean simulated daily energy gains in kJ. |
| `costs` | Mean simulated daily energy costs in kJ. |

### `data/simdata/gaincostpoor_data.csv`

Simulation summary written by `fatscaling_model_analyses.jl` for caecum gut type, `zeta = 1.0`, and `mu = 10^-3.1`. It contains 30 simulation rows and is exported for `math/gain_costs_allometry.nb`.

| Column | Description |
|---|---|
| `mass` | Body mass in kilograms. |
| `gains` | Mean simulated daily energy gains in kJ. |
| `costs` | Mean simulated daily energy costs in kJ. |

### `data/simdata/gaincostpoor_zeta1_data.csv`

Simulation summary written by `fatscaling_model_analyses.jl` for caecum gut type, `zeta = 1.0`, and `mu = 10^-2.8`. Columns are `mass`, `gains`, and `costs` as described above.

### `data/simdata/gaincostpoor_zeta15_data.csv`

Simulation summary written by `fatscaling_model_analyses.jl` for caecum gut type, `zeta = 1.5`, and `mu = 10^-2.8`. Columns are `mass`, `gains`, and `costs` as described above.

### `data/simdata/gaincostpoor_zeta20_data.csv`

Simulation summary written by `fatscaling_model_analyses.jl` for caecum gut type, `zeta = 2.0`, and `mu = 10^-2.8`. Columns are `mass`, `gains`, and `costs` as described above.

### `data/simdata/gaincostpoor_zeta215_data.csv`

Simulation summary written by `fatscaling_model_analyses.jl` for caecum gut type, `zeta = 2.15`, and `mu = 10^-2.8`. Columns are `mass`, `gains`, and `costs` as described above.

### `data/simdata/g1c1coords_zeta.csv`

Simulation-derived four-row coordinate table written by `fatscaling_model_analyses.jl` and imported into `math/reserveratio.nb`.

| Column | Description |
|---|---|
| `zeta` | Patchiness exponent. |
| `g1` | Gain exponent at the richness value where the gain-surplus exponent is closest to the selected empirical target `gammastar` in the script. |
| `c1` | Cost exponent at the same richness value. |
| `reserveratiodem` | Demographic reserve-ratio value computed in the script for the corresponding `zeta`. |

### `data/simdata/expgainremainderslope.CSV`

Two-column Mathematica-generated CSV with 351 numeric rows and no descriptive text header. With the script's default `CSV.read` call, the first numeric row is treated as column names and the remaining 350 rows are converted to an array for plotting the analytical expected-slope curve.

| Column | Description |
|---|---|
| column 1 | Resource richness `mu` on the original scale. The Julia script plots `log10(column 1)`. |
| column 2 | Analytical expected gain-surplus exponent/slope generated by `math/gain_costs_allometry.nb`. |

## Software and Package Versions

The local project was verified with Julia 1.11.5. The checked-in `Manifest.toml` was generated with Julia 1.11.4 and records `LinearAlgebra` 1.11.0, so Julia 1.11 is the intended Julia series for this repository.

The core package dependencies declared in `Project.toml` are:

| Package | Version/compatibility in this project |
|---|---|
| `DataFrames` | 1.7.0 |
| `Distributions` | 0.25.120 |
| `LinearAlgebra` | 1.11.0 |

The analysis scripts additionally import the following Julia packages. These versions were resolved when running with the `fatscaling` project active and the author's Julia v1.11 environment available on the Julia load path:

| Package | Version |
|---|---|
| `CSV` | 0.10.15 |
| `CategoricalArrays` | 0.10.8 |
| `ColorSchemes` | 3.29.0 |
| `Colors` | 0.12.11 |
| `DataFrames` | 1.7.0 |
| `Distributions` | 0.25.120 |
| `GLM` | 1.9.0 |
| `JLD2` | 0.5.13 |
| `LaTeXStrings` | 1.4.0 |
| `LinearAlgebra` | 1.11.0 |
| `Measures` | 0.3.2 |
| `MixedModels` | 4.35.2 |
| `Plots` | 1.40.13 |
| `ProgressMeter` | 1.10.4 |
| `RCall` | 0.14.8 |
| `Revise` | 3.8.0 |
| `StatsModels` | 0.7.4 |
| `StatsPlots` | 0.15.7 |
| `UnicodePlots` | 3.7.2 |

The empirical PGLS block uses R through `RCall`. The R environment used here was:

| Software/package | Version |
|---|---|
| R | 4.3.1 |
| `rotl` | 3.1.0 |
| `ape` | 5.8.1 |
| `caper` | 1.0.3 |
| `broom` | 1.0.8 |
| `phylolm` | 2.6.5 |

Install these R packages, or equivalent versions listed above, before running the empirical PGLS block on a clean machine.

```bash
Rscript -e 'install.packages(c("rotl", "ape", "caper", "broom", "phylolm"))'
```

The Mathematica notebooks record these front-end versions:

| Notebook | Recorded Wolfram front end |
|---|---|
| `math/gain_costs_allometry.nb` | 14.3 for Mac OS X x86 (64-bit), July 8, 2025 |
| `math/reserveratio.nb` | 14.3 for Mac OS X x86 (64-bit), July 8, 2025 |
| `math/storage_diff_drift_approx.nb` | 14.2 for Mac OS X x86 (64-bit), March 16, 2025 |

`wolframscript -version` reported WolframScript 1.12.0 for Mac OS X ARM (64-bit).

## Running the Analyses

Run commands from the directory that contains this README. If you are one level above the repository directory, enter it first:

```bash
cd fatscaling
```

The author's local shell helper for opening Julia with a selected thread count and project is:

```bash
jt() {
  local threads proj

  if [[ $# -gt 0 && "$1" =~ ^[0-9]+$ ]]; then
    threads=$1
    shift
  else
    threads=${JULIA_NUM_THREADS:-1}
  fi

  export JULIA_NUM_THREADS="$threads"
  echo "JULIA_NUM_THREADS=$JULIA_NUM_THREADS"

  if [[ $# -gt 0 && "$1" != -* ]]; then
    proj="--project=$1"
    shift
  else
    proj="--project"
  fi

  julia $proj "$@"
}
```

For example, from inside `fatscaling/`, `jt 8 .` is equivalent to launching Julia with eight threads and `--project=.`.

To instantiate the checked-in project environment:

```bash
jt 8 . -e 'using Pkg; Pkg.instantiate()'
```

On a clean machine, the additional script-only Julia packages listed above must also be available either in the active project or in the Julia v1.11 stacked environment. To install them into the active project, run:

```bash
jt 8 . -e 'using Pkg; Pkg.add(["CSV", "CategoricalArrays", "ColorSchemes", "Colors", "GLM", "JLD2", "LaTeXStrings", "Measures", "MixedModels", "Plots", "ProgressMeter", "RCall", "Revise", "StatsModels", "StatsPlots", "UnicodePlots"])'
```

To run the empirical analyses:

```bash
jt 8 . scripts/fatscaling_empirical_analyses.jl
```

To run the model analyses:

```bash
jt 8 . scripts/fatscaling_model_analyses.jl
```

The model script writes/loads `../external_data/gainsremainder_acrossgut.jld2` relative to the package root. The script comments indicate that this file is about 2 GB and that users may need to create the sibling `external_data/` directory or adjust the path before running the full model workflow.

## Citation

If you use this repository, please cite the associated manuscript. DOI and public repository citation should be added here after acceptance and public archiving.

## Example Run

The following one-off Julia example runs a single daily foraging simulation using the public `fatscaling` API.

```julia
using fatscaling
using Distributions
using LinearAlgebra
using DataFrames

# ----------------------- choose a single parameter set -----------------------
mass      = 35.0          # kg
gut       = "colon"       # "caecum", "colon", "non-rumen foregut", or "rumen foregut"
teeth     = "all"         # average tooth type determines chew rate
μ         = 1e-2          # resource richness, grams m^-2
ζ         = 1.0           # patchiness exponent
α         = 4.0           # Gamma shape parameter
edensity  = 18.2          # kJ g^-1 dry matter

# --- mass-scaling traits -----------------------------------------------------
β            = bite_size_allo(mass)                 # bite size, g bite^-1
chewrate     = chew_allo(mass, teeth)               # g s^-1
tchew        = β / chewrate                         # s bite^-1
maxgut       = gut_capacity_g(mass, gut)            # g
bcost, fcost = metabolic_cost(mass)                 # kJ s^-1, basal and foraging
velocity     = find_velocity(mass)                  # m s^-1
tmax_bout, _ = foragingtime(mass) .* 3600           # s d^-1

# --- competitor load and bite-encounter distribution -------------------------
_, n      = indperarea(mass)                        # conspecifics sharing patch
width     = reactionwidth(mass)
height    = reactionheight(mass)
na        = n / (width * height)                    # inds m^-2 in reaction plane

m_res     = μ / β                                   # bites m^-2
m′        = m_res / na                              # competitor-scaled bite density
α′        = α * na^(ζ - 2)                          # patchiness-scaled shape
γ_dist    = Gamma(α′, m′ / α′)                      # Gamma encounter distribution

# --- run the foraging simulation --------------------------------------------
gain, cost, n_bites = dailyforage(
    γ_dist, tchew, β, maxgut, velocity,
    bcost, fcost, edensity, tmax_bout
)
net = gain - cost

tbl = DataFrame(
    :energy_gain_kJ => round(gain; digits = 2),
    :energy_cost_kJ => round(cost; digits = 2),
    :net_energy_kJ  => round(net; digits = 2),
    :bites_taken    => n_bites,
)

display(tbl)
```
