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



###########################################################
# Simulations to generate primary results for Figures 3,4,5
###########################################################

#SAME THING BUT ACROSS GUT TYPES

# Define your gut types
gut_types = ["caecum", "colon", "non-rumen foregut", "rumen foregut"]
n_gut = length(gut_types)

# Your pre-existing parameter definitions
massexpvec = collect(1.5:0.1:4.4);
massvec = 10 .^ massexpvec;
l_massvec = length(massvec);

teeth = "all";  # remains constant here

# RESOURCE parameters
alpha = 4;
edensity = 18.2;
# p_bad = 0.05;
# configurations = 20000;
# runs = 200;

# muexpvec = collect(-8.0:0.01:-5.7);
# muexpvec = collect(-3:0.02:0);
muexpvec = collect(-3.5:0.02:-1);
l_muexpvec = length(muexpvec);
zetavec = [1.0, 1.5, 2.0, 2.15];
l_zetavec = length(zetavec);

reps = 1000

# Pre-allocate an array to store slopes for each gut type and each muexp value:
remainder_mu_slope_all = Array{Float64}(undef, n_gut, l_muexpvec);

gains_array = Array{Array{Float64}}(undef,n_gut);
costs_array = Array{Array{Float64}}(undef,n_gut);
encounters_array = Array{Array{Float64}}(undef,n_gut);
cons_maxgut_array = Array{Array{Float64}}(undef,n_gut);
gaindiff_array = Array{Array{Float64}}(undef,n_gut);
gainremainder_array = Array{Array{Float64}}(undef,n_gut);

# Loop over each gut type
for g in 1:n_gut
    current_gut_type = gut_types[g]
    println("Processing gut type: $current_gut_type")
    
    # Pre-allocate arrays for this gut type
    gains = Array{Float64}(undef, reps, l_massvec, l_muexpvec, l_zetavec);
    costs = Array{Float64}(undef, reps, l_massvec, l_muexpvec, l_zetavec);
    encounters = Array{Float64}(undef, reps, l_massvec, l_muexpvec, l_zetavec);
    cons_maxgut = Array{Float64}(undef, l_massvec);
    gaindiff = Array{Float64}(undef, reps, l_massvec, l_muexpvec, l_zetavec);
    gainremainder = Array{Float64}(undef, reps, l_massvec, l_muexpvec, l_zetavec);
    
    # Run simulation over reps, mass, muexp, and zeta
    @showprogress 1 "Computing for gut type $current_gut_type ..." for r in 1:reps
        @threads for m in 1:l_massvec
            for i in 1:l_muexpvec
                for j in 1:l_zetavec
                    mass = massvec[m]
                    muexp = muexpvec[i]
                    zeta = zetavec[j]
                    mu = 10^muexp
                    
                    # Compute consumer and resource properties
                    # Bite size gram/bite
                    beta = bite_size_allo(mass)
                    # Chew rate
                    chewrate = chew_allo(mass, teeth)
                    # Chew time sec per gram
                    t_chewgram = 1 / chewrate
                    # Chew time sec per bite
                    tchew = t_chewgram * beta
                    # Gut capacity grams
                    maxgut = gut_capacity_g(mass, current_gut_type)
                    # Energetic costs kJ per sec
                    bcost_kJps, fcost_kJps = metabolic_cost(mass)
                    # Velocity m/s
                    velocity = find_velocity(mass)
                    # Max bout time sec/day
                    tmax_bout, _ = foragingtime(mass) .* (60 * 60)
                    
                    # Competitor Load (inds)
                    _, n = indperarea(mass)

                    #Reaction plane area
                    width = reactionwidth(mass)
                    height = reactionheight(mass)

                    #Competitor density (per reaction plane area)
                    na = n/(width*height);
                    
                    # Bite density field
                    m_res =  mu * (1 / beta) #* height #* width

                    #Rescaling bite density by Competitor Load
                    mprime = m_res / na
                    alphaprime = alpha * na^(zeta - 2)
                    
                    #Rescaled bite encounter distribution
                    gammadist = Gamma(alphaprime, mprime / alphaprime)
                    
                    # Compute daily foraging outcomes
                    gains_daily, costs_daily, encounters_daily = dailyforage(gammadist, tchew, beta, maxgut, velocity, bcost_kJps, fcost_kJps, edensity, tmax_bout)
                    
                    #Save values
                    gains[r,m,i,j] = gains_daily
                    costs[r,m,i,j] = costs_daily
                    encounters[r,m,i,j] = encounters_daily
                    cons_maxgut[m] = maxgut*edensity

                    gaindiff[r,m,i,j] = gains[r,m,i,j] - cons_maxgut[m]

                    #Calculate gain surplus ***
                    gainremainder[r,m,i,j] = gains[r,m,i,j] - costs[r,m,i,j]
                end
            end
        end
    end

    gains_array[g] = gains
    costs_array[g] = costs
    encounters_array[g] = encounters
    cons_maxgut_array[g] = cons_maxgut
    gaindiff_array[g] = gaindiff
    gainremainder_array[g] = gainremainder
    
end


#EXAMINE SLOPE EXTRACTION
remainder_mu_slope_all = Array{Float64}(undef, n_gut, l_muexpvec);
minmass_array = Array{Array{Float64}}(undef,n_gut)
for g in 1:n_gut
    gains = gains_array[g];
    costs = costs_array[g];
    encounters = encounters_array[g];
    cons_maxgut = cons_maxgut_array[g];
    gaindiff = gaindiff_array[g];
    gainremainder = gainremainder_array[g];

    
    # Calculate slope for each muexp value
    remainder_mu_slope = Array{Float64}(undef, l_muexpvec)
    minmass = Array{Float64}(undef, l_muexpvec)
    for i in 1:l_muexpvec
        # Take mean across reps (and using the first zeta value)
        mean_gainremainder = vec(mean(gainremainder, dims=1)[:, :, i, 1])
        positive_gainremainder = findall(x -> x > 0, mean_gainremainder)
        if !isempty(positive_gainremainder)
            minmass[i] = massvec[positive_gainremainder][1]
        else
            minmass[i] = NaN
        end
        
        # Only perform the linear fit if there are enough positive points
        num_est = 20
        if length(positive_gainremainder) > num_est
            x = log.(massvec[positive_gainremainder])[end-num_est:end]
            y = log.(mean_gainremainder[positive_gainremainder])[end-num_est:end]
            df = DataFrame(x = x, y = y)
            model = lm(@formula(y ~ x), df)
            # Extract coefficients: intercept and slope
            intercept, slope = coef(model)
            remainder_mu_slope[i] = slope
        else
            remainder_mu_slope[i] = NaN
        end
    end
    minmass_array[g] = minmass
    # Save the slope curve for the current gut type
    remainder_mu_slope_all[g, :] = remainder_mu_slope
end

#What is the minimum mu where ANY species can have a positive gain remainder?
minmuany_array = Array{Float64}(undef,n_gut)
for g = 1:n_gut
    minmuany = muexpvec[findall(!isnan,minmass_array[1])[1]]
    minmuany_array[g] = minmuany
end

# 2 GB file

# filename=smartpath("data/simdata/gainsremainder_acrossgut.jld2")

# Pull file from external location
# NOTE: Delete prior to archiving
filename = string(homedir(),"/Dropbox/PostDoc/2024_herbforaging/herbforagingsim/data/simdata/gainsremainder_acrossgut.jld2");

# @save filename gut_types n_gut massexpvec massvec l_massvec teeth alpha edensity muexpvec l_muexpvec zetavec l_zetavec reps gains_array costs_array encounters_array cons_maxgut_array gaindiff_array gainremainder_array remainder_mu_slope_all minmass_array 


###########################################################
# OR Load Data to Create Figures 3,4,5
###########################################################

@load filename gut_types n_gut massexpvec massvec l_massvec teeth alpha edensity muexpvec l_muexpvec zetavec l_zetavec reps gains_array costs_array encounters_array cons_maxgut_array gaindiff_array gainremainder_array remainder_mu_slope_all minmass_array 

# filepath = string(homedir(),"/Dropbox/PostDoc/2024_herbforaging/herbforagingsim/data/expgainremainderslope.CSV")
filename = smartpath("data/expgainremainderslope.CSV")
expgainremainder_poor = CSV.read(filename, DataFrame)
expgainremainder_poor = Array(expgainremainder_poor)


###########################################################
# GENERATE FIGURE 3
###########################################################

gut_types_cap = ["Caecum","Colon","Non-rumen foregut","Rumen foregut"]
gut_type_maxgutslope = [0.860,0.919,0.881,0.897];
pgutslope = plot(
    xlabel = "Empirical gut capacity exponent",
    ylabel = L"Gain surplus exponent, $g_S \propto M^\gamma$",
    xlims = (0.8, 1.0),
    ylims = (0.8, 1.0),
    legend = :topleft,
    frame = :box,
    foreground_color_legend = nothing
              )
colors = palette(:tab10) #, n_gut)  # Get exactly n_gut colors
for g in 1:n_gut
    col = colors[g]
    gutslope = gut_type_maxgutslope[g]
    x = log.(massvec)
    y = log.(vec(mean(costs_array[g],dims=1)[:,:,end,1]))
    df = DataFrame(x = x, y = y)
    model = lm(@formula(y ~ x), df)
    intc, _ = coef(model)
    x = log.(massvec)
    y = log.(vec(mean(gains_array[g],dims=1)[:,:,end,1]))
    df = DataFrame(x = x, y = y)
    model = lm(@formula(y ~ x), df)
    intg, g_slope = coef(model)
    int = exp(intg)/exp(intc)
    c_slope =0.75;
    # exp_slope = c_slope + ((int*(gutslope-c_slope))*massvec[end]^(gutslope-c_slope))/(int*massvec[end]^(gutslope-c_slope) - 1)
    # From Mathematica Solution
    M = 100000; #Asymptotic Mass
    # exp_slope2 = g_slope + (intc*(c_slope - g_slope)*M^c_slope)/(intc*M^c_slope - intg*M^g_slope)
    exp_slope2 = gutslope + (intc*(c_slope - gutslope)*M^c_slope)/(intc*M^c_slope - intg*M^gutslope)
    scatter!(pgutslope, [gutslope], [remainder_mu_slope_all[g, end]],markersize=5,
    color=col,
    label=false,
    m = :rect)
    scatter!(pgutslope, [gutslope], [exp_slope2],markersize=5,color=col,label=gut_types_cap[g])
    plot!(pgutslope,
      [gutslope, gutslope],            # x coordinates (constant)
      [exp_slope2, gutslope],          # y coordinates (from exp_slope2 up to gutslope)
      lw=2,                          # line width
      alpha=0.5,
      color=col,                       # same color if you like
      label=false)                     # no extra legend entry
end
plot!(pgutslope,collect(0.85:0.01:0.95),collect(0.85:0.01:0.95),
    color="black",
    width=2,
    label=false)
scatter!(pgutslope,
    [0, 0],            # x coordinates (constant)
    [0, 0],          # y coordinates (from exp_slope2 up to gutslope)                         # line width
    alpha=1,
    color="white",                       # same color if you like
    markerstrokecolor="white",
    label=" ")    
scatter!(pgutslope,
    [0, 0],            # x coordinates (constant)
    [0, 0],          # y coordinates (from exp_slope2 up to gutslope)                         # line width
    alpha=1,
    color="black",                       # same color if you like
    label="Analytical approx.")                     # no extra legend entry
scatter!(pgutslope,
    [0, 0],            # x coordinates (constant)
    [0, 0],          # y coordinates (from exp_slope2 up to gutslope)                         # line width
    alpha=1,
    m = :rect,
    color="black",                       # same color if you like
    label="Simulation")                     # no extra legend entry
pgutslope
figfile = smartpath("rawfigures/fig_gutprediction.pdf")
Plots.savefig(pgutslope, figfile)

# Plots.scatter(log10.(expgainremainder_poor[:,1]),(expgainremainder_poor[:,2]),xlims=[muexpvec[1],muexpvec[end]])


###########################################################
# GENERATE FIGURE 5
###########################################################

gut_types_cap = ["Caecum","Colon","Non-rumen foregut","Rumen foregut"]
gut_type_maxgutslope = [0.860,0.919,0.881,0.897];
pslope = plot(
    xlabel = L"Richness, $\log_{10}(\mu)$",
    ylabel = L"Gain surplus exponent, $g_S \propto M^\gamma$",
    ylims = (0.8, 1.3),
    legend = :topright,
    frame = :box,
    size=(500,350),
    foreground_color_legend = nothing
              )
colors = palette(:tab10) #, n_gut)  # Get exactly n_gut colors
for g in 1:n_gut
    col = colors[g]
    gutslope = gut_type_maxgutslope[g]
    x = log.(massvec)
    y = log.(vec(mean(costs_array[g],dims=1)[:,:,end,1]))
    df = DataFrame(x = x, y = y)
    model = lm(@formula(y ~ x), df)
    intc, _ = coef(model)
    x = log.(massvec)
    y = log.(vec(mean(gains_array[g],dims=1)[:,:,end,1]))
    df = DataFrame(x = x, y = y)
    model = lm(@formula(y ~ x), df)
    intg, g_slope = coef(model)
    int = exp(intg)/exp(intc)
    c_slope =0.75;
    # exp_slope = c_slope + ((int*(gutslope-c_slope))*massvec[end]^(gutslope-c_slope))/(int*massvec[end]^(gutslope-c_slope) - 1)
    # From Mathematica Solution
    M = 100000; #Asymptotic Mass
    # exp_slope2 = g_slope + (intc*(c_slope - g_slope)*M^c_slope)/(intc*M^c_slope - intg*M^g_slope)
    exp_slope2 = gutslope + (intc*(c_slope - gutslope)*M^c_slope)/(intc*M^c_slope - intg*M^gutslope)
    plot!(muexpvec, remainder_mu_slope_all[g, :], label=gut_types_cap[g], lw=2,color=col)
    scatter!(pslope, [muexpvec[end]], [exp_slope2],markersize=5,color=col,label="")
end
gainpos = log10.(expgainremainder_poor[:,1]) .> -3.2;
Plots.plot!(log10.(expgainremainder_poor[gainpos,1]),(expgainremainder_poor[gainpos,2]),
    xlims=[muexpvec[1],muexpvec[end]+0.1],
    color=:black,
    label="Expected slope",
    width = 2,
    linestyle = :dot
    # markersize=3,
    )
# PLOT EMPIRICAL FAT-MASS EXPONENTS
valid_indices = findall(!isnan, remainder_mu_slope_all[1, :])
lindstedtvalue = 1.19; 
lindstedtcoords = [muexpvec[valid_indices[findmin(abs.(remainder_mu_slope_all[1, valid_indices] .- lindstedtvalue))[2]]], lindstedtvalue]
scatter!(pslope,[lindstedtcoords[1]],[lindstedtcoords[2]],
    markersize=6,
    color="gray",
    label="Previous")
    fullsetvalue = 1.17;
fullsetcoords = [muexpvec[valid_indices[findmin(abs.(remainder_mu_slope_all[1, valid_indices] .- fullsetvalue))[2]]], fullsetvalue]
scatter!(pslope,[fullsetcoords[1]],[fullsetcoords[2]],
    markersize=6,
    color="black",
    label="This paper")
captivevalue = 0.95;
captivecoords = [muexpvec[valid_indices[findmin(abs.(remainder_mu_slope_all[1, valid_indices] .- captivevalue))[2]]], captivevalue]
captivecoords2 = [muexpvec[valid_indices[findmin(abs.(remainder_mu_slope_all[3, valid_indices] .- captivevalue))[2]]], captivevalue]
scatter!(pslope,[captivecoords[1]],[captivecoords[2]],
    markersize=6,
    color="white",
    label="Captive")
plot!(pslope,[captivecoords[1]*0.99,captivecoords2[1]*1.01], [captivecoords[2],captivecoords2[2]],
    width=2,
    color="black",
    label=false)
scatter!(pslope,[captivecoords2[1]],[captivecoords2[2]],
    markersize=6,
    color="white",
    label=false)

display(pslope)

# SAVE FIGURE 4
figfile = smartpath("rawfigures/fig_gainsremainder_slopev2.pdf")
Plots.savefig(pslope, figfile)



#######################################################################
# ORGANIZE DATA FOR FIGURE 4... 
########################################################################
# EXPORT AND THEN IMPORT W/MATHEMATICA FILE: gains_costs_allometry_v2.nb
########################################################################


muvalue = -3.1
mupos = findall(x->x==muvalue,muexpvec)[1]; #1:231
10 .^muexpvec[mupos]
guti = 1
zetai = 1
gainsi = vec(mean(gains_array[guti],dims=1)[:,:,mupos,zetai]);
costsi = vec(mean(costs_array[guti],dims=1)[:,:,mupos,zetai]);
x = log.(massvec);
y = log.(gainsi);
df = DataFrame(x = x, y = y);
model = lm(@formula(y ~ x), df)
plot(log.(massvec),log.(gainsi))
plot!(log.(massvec),log.(costsi))

dfi = DataFrame(mass = massvec,gains = gainsi, costs = costsi)
filename = smartpath("data/simdata/gaincostpoor_data.csv")
CSV.write(filename, dfi)


#For zeta = 1.0
# muvalue = -3.1
muvalue = -2.8
mupos = findall(x->x==muvalue,muexpvec)[1]; #1:231
10 .^muexpvec[mupos]
guti = 1
zetai = 1
gainsi = vec(mean(gains_array[guti],dims=1)[:,:,mupos,zetai]);
costsi = vec(mean(costs_array[guti],dims=1)[:,:,mupos,zetai]);
x = log.(massvec);
y = log.(gainsi);
df = DataFrame(x = x, y = y);
model = lm(@formula(y ~ x), df)
plot(log.(massvec),log.(gainsi))
plot!(log.(massvec),log.(costsi))

dfi = DataFrame(mass = massvec,gains = gainsi, costs = costsi)
filename = smartpath("data/simdata/gaincostpoor_zeta1_data.csv")
CSV.write(filename, dfi)

#For zeta = 1.5
# muvalue = -3.1
muvalue = -2.8
mupos = findall(x->x==muvalue,muexpvec)[1]; #1:231
10 .^muexpvec[mupos]
guti = 1
zetai = 2
gainsi = vec(mean(gains_array[guti],dims=1)[:,:,mupos,zetai]);
costsi = vec(mean(costs_array[guti],dims=1)[:,:,mupos,zetai]);
x = log.(massvec);
y = log.(costsi);
df = DataFrame(x = x, y = y);
model = lm(@formula(y ~ x), df)
plot(log.(massvec),log.(gainsi))
plot!(log.(massvec),log.(costsi))

dfi = DataFrame(mass = massvec,gains = gainsi, costs = costsi)
filename = smartpath("data/simdata/gaincostpoor_zeta15_data.csv")
CSV.write(filename, dfi)


#For zeta = 2.0
# muvalue = -3.1
muvalue = -2.8
mupos = findall(x->x==muvalue,muexpvec)[1]; #1:231
10 .^muexpvec[mupos]
guti = 1
zetai = 3
gainsi = vec(mean(gains_array[guti],dims=1)[:,:,mupos,zetai]);
costsi = vec(mean(costs_array[guti],dims=1)[:,:,mupos,zetai]);
x = log.(massvec);
y = log.(costsi);
df = DataFrame(x = x, y = y);
model = lm(@formula(y ~ x), df)
plot(log.(massvec),log.(gainsi))
plot!(log.(massvec),log.(costsi))

dfi = DataFrame(mass = massvec,gains = gainsi, costs = costsi)
filename = smartpath("data/simdata/gaincostpoor_zeta20_data.csv")
CSV.write(filename, dfi)


#For zeta = 2.15
muvalue = -2.8
mupos = findall(x->x==muvalue,muexpvec)[1]; #1:231
10 .^muexpvec[mupos]
guti = 1
zetai = 4
gainsi = vec(mean(gains_array[guti],dims=1)[:,:,mupos,zetai]);
costsi = vec(mean(costs_array[guti],dims=1)[:,:,mupos,zetai]);
x = log.(massvec);
y = log.(costsi);
df = DataFrame(x = x, y = y);
model = lm(@formula(y ~ x), df)
plot(log.(massvec),log.(gainsi))
plot!(log.(massvec),log.(costsi))

dfi = DataFrame(mass = massvec,gains = gainsi, costs = costsi)
filename = smartpath("data/simdata/gaincostpoor_zeta215_data.csv")
CSV.write(filename, dfi)




# EVALUATE WHERE gamma = 1.19
#For zeta = 1.0
# muvalue = -3.1
muvalue = -3.06
mupos = findall(x->x==muvalue,muexpvec)[1]; #1:231
10 .^muexpvec[mupos]
guti = 1
zetai = 1
gainsi = vec(mean(gains_array[guti],dims=1)[:,:,mupos,zetai]);
costsi = vec(mean(costs_array[guti],dims=1)[:,:,mupos,zetai]);
x = log.(massvec);
y = log.(costsi);
df = DataFrame(x = x, y = y);
model = lm(@formula(y ~ x), df)
plot(log.(massvec),log.(gainsi))
plot!(log.(massvec),log.(costsi))

#For zeta = 1.5
# muvalue = -3.1
muvalue = -3.06
mupos = findall(x->x==muvalue,muexpvec)[1]; #1:231
10 .^muexpvec[mupos]
guti = 1
zetai = 2
gainsi = vec(mean(gains_array[guti],dims=1)[:,:,mupos,zetai]);
costsi = vec(mean(costs_array[guti],dims=1)[:,:,mupos,zetai]);
x = log.(massvec);
y = log.(costsi);
df = DataFrame(x = x, y = y);
model = lm(@formula(y ~ x), df)
plot(log.(massvec),log.(gainsi))
plot!(log.(massvec),log.(costsi))


#For zeta = 2.0
# muvalue = -3.1
muvalue = -3.06
mupos = findall(x->x==muvalue,muexpvec)[1]; #1:231
10 .^muexpvec[mupos]
guti = 1
zetai = 3
gainsi = vec(mean(gains_array[guti],dims=1)[:,:,mupos,zetai]);
costsi = vec(mean(costs_array[guti],dims=1)[:,:,mupos,zetai]);
x = log.(massvec);
y = log.(gainsi);
df = DataFrame(x = x, y = y);
model = lm(@formula(y ~ x), df)
plot(log.(massvec),log.(gainsi))
plot!(log.(massvec),log.(costsi))







enci = vec(mean(encounters_array[guti],dims=1)[:,:,mupos,zetai]);
plot(log.(massvec),log.(enci))

plot(log.(massvec),(gainsi - costsi))

gainremainderi = gainsi .- costsi;
plot(log.(massvec[gainremainderi .> 0]),log.(gainremainderi[gainremainderi .> 0]))

x = log.(massvec);
y = log.(costsi);
df = DataFrame(x = x, y = y);
model = lm(@formula(y ~ x), df);
println(coef(model))

x = log.(massvec);
y = log.(gainsi);
df = DataFrame(x = x, y = y);
model = lm(@formula(y ~ x), df);
println(coef(model))


###################################################################
# DEMOGRAPHIC RESERVE RATIO CALCULATIONS (doesn't take long to run)
###################################################################


gammastar = 1.17;

#Calculate demographic costs:
vals = carryingcapacity_costs.(massvec)
E_kcosts = getindex.(vals,1)
t_lifetime = getindex.(vals,2)
E_kcosts_perdiem = E_kcosts ./ t_lifetime
plot(log.(massvec),log.(E_kcosts_perdiem))
x = log.(massvec)
y = log.(E_kcosts_perdiem)
df = DataFrame(x = x, y = y)
model = lm(@formula(y ~ x), df)
c_k0 = coef(model)[1]

gaincostratio_raw = Array{Float64}(undef,length(gut_types),length(zetavec),length(muexpvec));
gaincostratio = Array{Float64}(undef,length(gut_types),length(zetavec),length(muexpvec));
slopeg = Array{Float64}(undef,length(gut_types),length(zetavec),length(muexpvec));
slopec = Array{Float64}(undef,length(gut_types),length(zetavec),length(muexpvec));
slopegr = Array{Float64}(undef,length(gut_types),length(zetavec),length(muexpvec));
for g=1:length(gut_types)
    for j=1:length(zetavec)
        for i=1:length(muexpvec)
            gainsi = vec(mean(gains_array[g],dims=1)[:,:,i,j]);
            costsi = vec(mean(costs_array[g],dims=1)[:,:,i,j]);

            gremainder = gainsi - costsi;
            
            xc = log.(massvec)
            yc = log.(costsi)
            dfc = DataFrame(xc = xc, yc = yc)
            modelc = lm(@formula(yc ~ xc), dfc)

            xg = log.(massvec)
            yg = log.(gainsi)
            dfg = DataFrame(xg = xg, yg = yg)
            modelg = lm(@formula(yg ~ xg), dfg)
            
            if length(massvec[gremainder .> 0]) > 20
                xgr = log.(massvec[gremainder .> 0][end-20:end])
                ygr = log.(gremainder[gremainder .> 0][end-20:end])
                dfgr = DataFrame(xgr = xgr, ygr = ygr)
                modelgr = lm(@formula(ygr ~ xgr), dfgr)
                gr1 = coef(modelgr)[2]
                slopegr[g,j,i] = gr1
            else
                slopegr[g,j,i] = Inf
            end


            c1 = coef(modelc)[2]
            c0 = coef(modelc)[1]
            
            g1 = coef(modelg)[2]
            

            slopeg[g,j,i] = g1
            slopec[g,j,i] = c1

            #NOT Accounting fo daily demographic costs:
            gaincostratio_raw[g,j,i] = 1 - (c1-g1)/(gammastar-g1) 
            
            #Accounting fo daily demographic costs:
            gaincostratio[g,j,i] = 1 - (c1-g1)/(gammastar-g1) - (exp(c_k0)/exp(c0)) 
        end
    end
end

# DOES NOT CHANGE ACROSS GUT TYPE

#Calculate the gaincostratio exactly where gamma = 1.17

# WITHOUT DEMOGRAPHIC COSTS
gaincostratiostarraw = [gaincostratio_raw[1,j,findmin(abs.(slopegr[1,j,:] .- gammastar))[2]] for j=1:length(zetavec)]

# WITH DEMOGRAPHIC COSTS
gaincostratiostar = [gaincostratio[1,j,findmin(abs.(slopegr[1,j,:] .- gammastar))[2]] for j=1:length(zetavec)]


#what is c1 and g1 at posstar?
# g1
slopeg[1,1,findmin(abs.(slopegr[1,1,:] .- gammastar))[2]]
#c1
slopec[1,1,findmin(abs.(slopegr[1,1,:] .- gammastar))[2]]

#k0/c0
gaincostratiostarraw[1] - gaincostratiostar[1]


#######################################################################
# Gain/Cost Exponent Coordinates for FIGURE 6: Import into Mathematica File: reserveratio.nb
########################################################################

#Create a table of g1,c1 coordinates
g1coords = [slopeg[1,j,findmin(abs.(slopegr[1,j,:] .- gammastar))[2]] for j=1:l_zetavec];
c1coords = [slopec[1,j,findmin(abs.(slopegr[1,j,:] .- gammastar))[2]] for j=1:l_zetavec];

g1c1df = DataFrame(zeta = zetavec, g1 = g1coords, c1 = c1coords, reserveratiodem = gaincostratiostar)   # names can be whatever you like

# EXPORT TO IMPORT INTO MATHEMATIC FILE: reserveratio.nb
filename = smartpath("data/simdata/g1c1coords_zeta.csv")
CSV.write(filename, g1c1df)



#Calculate the costs of growth and growing 2 offspring to 1/2 adult mass
#These are the assumed sunk-costs required of carrying capacity
vals = carryingcapacity_costs.(massvec)
E_kcosts = getindex.(vals,1)
t_lifetime = getindex.(vals,2)
E_kcosts_perdiem = E_kcosts ./ t_lifetime
kappaplot = plot(massvec,E_kcosts_perdiem,
    xscale = :log10,
    yscale = :log10,
    xlabel = "Body mass (kg)",
    ylabel = L"Per-diem demographic costs, $\kappa$ (kJ)",
    label=false,
    frame=:box,
    width=2)
figfile = smartpath("rawfigures/fig_kappa.pdf")
Plots.savefig(kappaplot, figfile)

x = log.(massvec)
y = log.(E_kcosts_perdiem)
df = DataFrame(x = x, y = y)
model = lm(@formula(y ~ x), df)
c_k0 = coef(model)[1]
println(coef(model))

# mu position 22 accords with gamma = 1.19
x = log.(massvec)
y = log.(vec(mean(costs_array[1],dims=1)[:,:,22,1]))
df = DataFrame(x = x, y = y)
model = lm(@formula(y ~ x), df)
println(coef(model))
c_0 = coef(model)[1]

(exp(c_k0)/exp(c_0)) #.* massvec .^ (0.69-0.78)



#######################################################################
# SUPPLEMENTARY FIGURES
########################################################################




#INVESTIGATE HOW ZETA INFLUENCE CV AND D as a function of body size


# Define your gut types
gut_types = ["caecum", "colon", "non-rumen foregut", "rumen foregut"]
n_gut = length(gut_types)

# Your pre-existing parameter definitions
massexpvec = collect(1.5:0.1:4.4);
massvec = 10 .^ massexpvec;
l_massvec = length(massvec);

teeth = "all";  # remains constant here

# RESOURCE parameters
alpha = 4;
edensity = 18.2;
# p_bad = 0.05;
# configurations = 20000;
# runs = 200;

# muexpvec = collect(-8.0:0.01:-5.7);
# muexpvec = collect(-3:0.02:0);
muexpvec = collect(-3.5:0.02:-1);
l_muexpvec = length(muexpvec);
zetavec = [1.0, 1.5, 2.0, 2.15];
l_zetavec = length(zetavec);

reps = 10000

# cv = Array{Float64}(undef,l_massvec)
m_ind = Array{Float64}(undef,l_zetavec,l_massvec)
var_ind = Array{Float64}(undef,l_zetavec,l_massvec)
cv_ind = Array{Float64}(undef,l_zetavec,l_massvec)
dist = Array{Float64}(undef,l_zetavec,l_massvec,reps)

for j=1:l_zetavec
    @threads for m=1:l_massvec
        current_gut_type = "caecum"
        mass = massvec[m]
        muexp = -2.8
        zeta = zetavec[j]
        mu = 10^muexp

        # Compute consumer and resource properties
        beta = bite_size_allo(mass)
        chewrate = chew_allo(mass, teeth)
        t_chewgram = 1 / chewrate
        tchew = t_chewgram * beta
        maxgut = gut_capacity_g(mass, current_gut_type)
        bcost_kJps, fcost_kJps = metabolic_cost(mass)
        velocity = find_velocity(mass)
        tmax_bout, _ = foragingtime(mass) .* (60 * 60)
        ndensity, n = indperarea(mass)
        width = reactionwidth(mass)
        height = reactionheight(mass)
        #adjacent competitor density (inds/m)
        na = n/(width*height);

        m_res =  mu * (1 / beta) #* height #* width
        mprime = m_res / na
        alphaprime = alpha * na^(zeta - 2)

        # \sigma^2 = q^2 / (n^{\zeta} \alpha)
        # var = m_res^2/((na^zeta)*alpha);
        # sigma = sqrt(var)
        # cv[m] = (sigma*na^(1 - (zeta/2)))/m_res

        # var_ind = mprime^2/((na^zeta)*alphaprime);
        # sigma_ind = sqrt(var_ind)

        var = m_res^2/(alpha);
        sigma  = sqrt(var)

        m_ind[j,m] = (alpha*na^(zeta-1))/(m_res*(alpha*na^(zeta-2)-1))

        var_ind[j,m] = ((alpha*na^(zeta-2))^3)/(((alpha*na^(zeta-2)-1)^2)*(alpha*na^(zeta-2) - 2)*(m_res/na)^2);
        
        cv_ind[j,m] = (sigma*na^(1 - (zeta/2)))/m_res
        # cv_ind[m] = (sigma_ind*na^(1 - (zeta/2)))/mprime

        gammadist = Gamma(alphaprime, mprime / alphaprime)

        for r=1:reps
            rg = rand(gammadist);
            distance_to_resource = rand(Exponential(1.0/rg));
            dist[j,m,r] = distance_to_resource
        end


    end
end

zetacolors = ["#4772b2", "#bf8b26", "#708c00", "#2e9999"];

#CV across body mass with increasing zeta
zplot = plot((massvec),cv_ind[1,:],
    ylims = [0,1],
    xlabel="Body mass (kg)",
    ylabel="CV resource density",
    xscale=:log10,
    color=zetacolors[1],
    frame=:box,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    legend=:topleft,
    label=L"\zeta = " * string(1.0),
    width=2)
[plot!(zplot,(massvec),cv_ind[i,:],
    color=zetacolors[i],
    label=L"\zeta = " * string(zetavec[i]),
    width=2
    ) 
    for i=2:4]
zplot

#Mean distances across body mass with increasing zeta
mzplot = scatter((massvec),mean(dist,dims=3)[1,:,:],
    xscale=:log10,
    color=zetacolors[1],
    frame=:box,
    xlabel="Body mass (kg)",
    ylabel="Mean distance (m)",
    label=L"\zeta = " * string(1.0),
    foreground_color_legend = nothing,
    width=2)
[scatter!(mzplot,(massvec),mean(dist,dims=3)[i,:,:],
    color=zetacolors[i],
    label = L"\zeta = " * string(zetavec[i]),
    width=2)
    for i=2:4]
[plot!(mzplot,(massvec),m_ind[i,:],
    color=zetacolors[i],
    markersize=3,
    label=false,
    width =2) 
    for i=1:4]
mzplot


#Std distances across body mass with increasing zeta
pos4 = var_ind[4,:] .> 0;
sdzplot = scatter((massvec),std(dist,dims=3)[1,:,:],
    xscale=:log10,
    yscale=:log10,
    color=zetacolors[1],
    frame=:box,
    xlabel="Body mass (kg)",
    ylabel="Std distance (m)",
    label=L"\zeta = " * string(1.0),
    foreground_color_legend = nothing,
    backgroud_color_legend = nothing,
    legend=:topleft,
    width=2)
[scatter!(sdzplot,(massvec),std(dist,dims=3)[i,:,:],
    color=zetacolors[i],
    label = L"\zeta = " * string(zetavec[i]),
    width=2)
    for i=2:4]
[plot!(sdzplot,
    i == 4 ? massvec[pos4] : massvec,
    i == 4 ? sqrt.(var_ind[i, pos4]) : sqrt.(var_ind[i, :]),
    color = zetacolors[i],
    markersize=3,
    label = false,
    width=2)
    for i = 1:4]
sdzplot

using Measures
zpanelplot = plot(
    zplot, mzplot, sdzplot, 
    layout = (1, 3), 
    size=(1200, 350),
    margin = 6mm,         # Space around entire figure
    plot_margin = 2mm)     # Space *inside* each subplot)

figfile = smartpath("rawfigures/fig_zetastats.pdf")
Plots.savefig(zpanelplot, figfile)




# RESERVE RATIO FOR A DETAILED ZETA > 2 WHEN MU = 10^-2.8
# FOR A SINGLE GUT TYPE

# Your pre-existing parameter definitions
massexpvec = collect(1.5:0.1:4.4);
massvec = 10 .^ massexpvec;
l_massvec = length(massvec);

teeth = "all";  # remains constant here
current_gut_type = "caecum";

# RESOURCE parameters
alpha = 4;
edensity = 18.2;
muexp = -2.8;
mu = 10^muexp
zetavec = collect(2.0:0.01:2.25);
l_zetavec = length(zetavec);

reps = 5000

# Pre-allocate arrays for this gut type
gains = Array{Float64}(undef, reps, l_massvec, l_zetavec);
costs = Array{Float64}(undef, reps, l_massvec, l_zetavec);
encounters = Array{Float64}(undef, reps, l_massvec, l_zetavec);
cons_maxgut = Array{Float64}(undef, l_massvec);
gaindiff = Array{Float64}(undef, reps, l_massvec, l_zetavec);
gainremainder = Array{Float64}(undef, reps, l_massvec, l_zetavec);

# Run simulation over reps, mass, muexp, and zeta
for r in 1:reps
    @threads for m in 1:l_massvec
        for j in 1:l_zetavec
            mass = massvec[m]
            zeta = zetavec[j]

            # Compute consumer and resource properties
            beta = bite_size_allo(mass)
            chewrate = chew_allo(mass, teeth)
            t_chewgram = 1 / chewrate
            tchew = t_chewgram * beta
            maxgut = gut_capacity_g(mass, current_gut_type)
            bcost_kJps, fcost_kJps = metabolic_cost(mass)
            velocity = find_velocity(mass)
            tmax_bout, _ = foragingtime(mass) .* (60 * 60)
            ndensity, n = indperarea(mass)
            width = reactionwidth(mass)
            height = reactionheight(mass)
            #adjacent competitor density (inds/m)
            na = n/(width*height);
            
            m_res =  mu * (1 / beta) #* height #* width
            mprime = m_res / na
            alphaprime = alpha * na^(zeta - 2)
            
            gammadist = Gamma(alphaprime, mprime / alphaprime)
            
            # Compute daily foraging outcomes
            gains_daily, costs_daily, encounters_daily = dailyforage(gammadist, tchew, beta, maxgut, velocity, bcost_kJps, fcost_kJps, edensity, tmax_bout)
            
            gains[r,m,j] = gains_daily
            costs[r,m,j] = costs_daily
            encounters[r,m,j] = encounters_daily
            cons_maxgut[m] = maxgut*edensity

            gaindiff[r,m,j] = gains[r,m,j] - cons_maxgut[m]

            #Calculate remainder gain - cost 
            gainremainder[r,m,j] = gains[r,m,j] - costs[r,m,j]
        end
    end
end

meangains = mean(gains,dims=1)[1,:,:]
meancosts = mean(costs,dims=1)[1,:,:]
meanreserveratio = meangains ./ meancosts

# zetacolors = ["#4772b2", "#bf8b26", "#708c00", "#2e9999"];
using ColorSchemes, Colors
# 1) Define your two end-point colors
c1 = colorant"#708c00"   # RGB(0.44,0.55,0)
c2 = colorant"#2e9999"   # RGB(0.18,0.6,0.6)
c3 = colorant"#8167e6"
# 2) Build a tiny ColorScheme from them
my_scheme = ColorScheme([c1, c2])
my_scheme2 = ColorScheme([c2, c3])
# 3) Sample 16 equally‐spaced colors along that gradient
colors1 = get(my_scheme, range(0, stop=1, length=16))
colors2 = get(my_scheme2, range(0, stop=1, length=10))
colors = [colors1; colors2]

plt = plot(massvec, meanreserveratio[:,1];
    xscale=:log10, 
    xlabel = "Body mass (kg)",
    ylabel = L"Reserve ratio, $\phi=g/c$",
    label=L"\zeta = " * string(zetavec[1]),
    frame=:box,
    color = colors[1],
    width = 2,
    ylims=[0,3],
    foreground_color_legend = nothing,
    backgroud_color_legend = nothing)
maxpos = findmax(meanreserveratio[:,1])
scatter!(plt,[massvec[maxpos[2]]],[maxpos[1]],
    color = colors[1],
    label = false)
scatter!(plt,[massvec[maxpos[2]]],[maxpos[1]],
        color = colors[1],
        markerstrokewidth = 2,
        markersize = 8,
        label = false)
for j in 2:size(meanreserveratio,2)
    maxpos = findmax(meanreserveratio[:,j])
    plot!(plt, massvec, meanreserveratio[:,j],
        label = (j == 16 || j == 26) ? L"\zeta = " * string(zetavec[j]) : false,
        color = colors[j],
        width = 2)
    scatter!(plt,[massvec[maxpos[2]]],[maxpos[1]],
        color = colors[j],
        label = false)
    #Make larger markers for values used in analysis
    if  j == 16
        scatter!(plt,[massvec[maxpos[2]]],[maxpos[1]],
        color = colors[j],
        markerstrokewidth = 2,
        markersize = 8,
        label = false)
        println([[massvec[maxpos[2]]],[maxpos[1]]])
    end
end
plot!(plt, massvec, repeat([1],outer=l_massvec),
    linestyle=:dot,
    label = false,
    color=:black,
    width=2)
display(plt)

masszetacoords = Array{Float64}(undef,l_zetavec,2)
for j = 1:l_zetavec
    maxpos = findmax(meanreserveratio[:,j])
    masszetacoords[j,:] = [maxpos[1],massvec[maxpos[2]]]
end
minfeasiblemasspos = findmin(abs.(masszetacoords[:,1] .- 1))[2]
minfeasiblemass = masszetacoords[minfeasiblemasspos,2]

plt2 = plot(zetavec,masszetacoords[:,2],
    yscale=:log10,
    frame=:box,
    xlabel=L"Patchiness, $\zeta$",
    ylabel="Body mass (kg)",
    legend=false,
    color=:black,
    width=2)
scatter!(plt2,zetavec,masszetacoords[:,2],
    yscale=:log10,
    frame=:box,
    xlabel=L"Patchiness exponent, $\zeta$",
    ylabel="Body mass (kg)",
    legend=false,
    color=colors)
plot!(plt2,zetavec,repeat([minfeasiblemass],outer=l_zetavec),
    width = 2,
    color=:black,
    label=false,
    linestyle=:dot)

    # assuming you already have `plt` and `plt2` defined...
using Measures
combined = plot(
    plt, plt2;
    layout       = (2,1),
    size         = (400, 700),
    margin = 2.5mm,         # Space around entire figure
    plot_margin = 1mm)     # Space *inside* each subplot)
display(combined)

figfile = smartpath("rawfigures/fig_reserveratio_highzeta.pdf")
Plots.savefig(combined, figfile)
















































###################################
# OLD
###################################





# Compressed single pane version
# PLOT
gut_types_cap = ["Caecum", "Colon", "Non-rumen foregut", "Rumen foregut"]
x₀ = -2.5
# colors = cgrad(:roma, 4, categorical = true)
colors = palette(:tab10) #, n_gut)  # Get exactly n_gut colors
#[:red, :blue, :green, :purple]
# 1) build your 2×2 canvas with NO per‐panel axis‐titles:
pratio = plot(
  layout                  = (1,1),
  frame                   = :box,
  foreground_color_legend = nothing,
  legend=:topright
)
zetavec2 = [1,3]; #Int64.(zetavec[[1,3]])
# 2) same double‐loop + your if‐logic for the legend
for (g, gut) in enumerate(gut_types_cap)
  for j in zetavec2
    println(j)
    # j = 1
    # create masks for the two line‐styles
    mask1 = muexpvec .<= x₀
    mask2 = muexpvec .>= x₀
    col = colors[g]
    gammapos = findmin(abs.(slopegr[g,j,:] .- gammastar))[2]
    # — solid up to x₀
    plot!(
      pratio,
      muexpvec[mask1], gaincostratio[g, j, mask1],
      subplot  = 1,
      ylabel   = L"Reserve ratio, $\phi^{\mathrm{dem}}$",
      xlabel   = L"Richness, $\log_{10}(\mu)$",
      ylims = [1.0,1.5],
      linestyle= :solid,
      color     = col, #(g == 2 ? "ζ = $z" : false),
      label    = (j == 1 ? gut_types_cap[g] : false),
      width    = 2
    )
    # — dashed after x₀
    plot!(
      pratio,
      muexpvec[mask2], gaincostratio[g, j, mask2],
      subplot   = 1,
      linestyle = :dash,
      color     = col,
      label     = false,
      width     = 2
    )
    scatter!(
        pratio,
        subplot  = 1,
        [muexpvec[gammapos]],[gaincostratio[g,j,gammapos]],
        color = "black",
        label = ""
    )
  end
end
Plots.annotate!(
      pratio,
      subplot = 1,
      -2.88, 1.36,
      text(L"\zeta=2", 15, halign = :left)
    )
Plots.annotate!(
    pratio,
    subplot = 1,
    -3.4, 1.36,
    text(L"\zeta=1", 15, halign = :left)
    )
pratio

figfile = smartpath("rawfigures/fig_gainsremainder_gcratiov2.pdf");
Plots.savefig(pratio, figfile)











plot!(mzplot,(massvec),mean(dist,dims=3)[2,:,:],
    col=zetacolors[1])
plot!(mzplot,(massvec),mean(dist,dims=3)[3,:,:],
    col=zetacolors[1])
plot!(mzplot,(massvec),mean(dist,dims=3)[4,:,:],
    col=zetacolors[1])
scatter!(mzplot,(massvec),m_ind[1,:],
col=zetacolors[1])
scatter!(mzplot,(massvec),m_ind[2,:])
scatter!(mzplot,(massvec),m_ind[3,:])
scatter!(mzplot,(massvec),m_ind[4,:])


#Mean STD across body mass with increasing zeta
pos4 = var_ind[4,:] .> 0;
sdzplot = plot((massvec),(std(dist,dims=3)[1,:,:]),xscale=:log10,yscale=:log10)
plot!(sdzplot,(massvec),(std(dist,dims=3)[2,:,:]))
plot!(sdzplot,(massvec),(std(dist,dims=3)[3,:,:]))
plot!(sdzplot,(massvec),(std(dist,dims=3)[4,:,:]))
scatter!(sdzplot,(massvec),(sqrt.(var_ind[1,:])))
scatter!(sdzplot,(massvec),(sqrt.(var_ind[2,:])))
scatter!(sdzplot,(massvec),(sqrt.(var_ind[3,:])))
scatter!(sdzplot,(massvec[pos4]),(sqrt.(var_ind[4,pos4])))









#export data table
masszetaCV = DataFrame(mass=massvec,CVzeta1=cv_ind[1,:],CVzeta15=cv_ind[2,:],CVzeta2=cv_ind[3,:],CVzeta215=cv_ind[4,:])
masszetaMean = DataFrame(mass=massvec,Meanzeta1=m_ind[1,:],Meanzeta15=m_ind[2,:],Meanzeta2=m_ind[3,:],Meanzeta215=m_ind[4,:])
masszetaVar = DataFrame(mass=massvec,Varzeta1=var_ind[1,:],Varzeta15=var_ind[2,:],Varzeta2=var_ind[3,:],Varzeta215=var_ind[4,:])






# GainRemainder Plot for HIGH ZETA
zetapos = 4;
gut_types_cap = ["Caecum","Colon","Non-rumen foregut","Rumen foregut"]
gut_type_maxgutslope = [0.860,0.919,0.881,0.897];
pslopehighz = plot(
    xlabel = L"Richness, $\mu$",
    ylabel = L"Forage excess exponent, $g_R \propto M^\gamma$",
    ylims = (0.8, 1.3),
    legend = :topright,
    frame = :box,
    foreground_color_legend = nothing
              )
colors = palette(:tab10) #, n_gut)  # Get exactly n_gut colors
for g in 1:n_gut
    col = colors[g]
    gutslope = gut_type_maxgutslope[g]
    x = log.(massvec)
    y = log.(vec(mean(costs_array[g],dims=1)[:,:,end,zetapos]))
    df = DataFrame(x = x, y = y)
    model = lm(@formula(y ~ x), df)
    intc, _ = coef(model)
    x = log.(massvec)
    y = log.(vec(mean(gains_array[g],dims=1)[:,:,end,zetapos]))
    df = DataFrame(x = x, y = y)
    model = lm(@formula(y ~ x), df)
    intg, g_slope = coef(model)
    int = exp(intg)/exp(intc)
    c_slope =0.75;
    # exp_slope = c_slope + ((int*(gutslope-c_slope))*massvec[end]^(gutslope-c_slope))/(int*massvec[end]^(gutslope-c_slope) - 1)
    # From Mathematica Solution
    M = 100000; #Asymptotic Mass
    exp_slope2 = g_slope + (intc*(c_slope - g_slope)*M^c_slope)/(intc*M^c_slope - intg*M^g_slope)
    plot!(pslopehighz,muexpvec, remainder_mu_slope_all[g, :], label=gut_types_cap[g], lw=2,color=col)
    scatter!(pslopehighz, [muexpvec[end]], [exp_slope2],markersize=5,color=col,label="")
end
#Find the minimum mu where a mass < 100 results in a positive gain remainder
#NOTE: Shouldn't minmass be a function of gut type? - no same across all! 04/07/25
minmu = muexpvec[findall(x->x<50,minmass_array[3])][1]
minmu_yvalues = collect(0.8:0.1:1.3)
minmu_xvalues = repeat([minmu], outer=length(minmu_yvalues))
# Add the line to the plot
Plots.plot!(pslopehighz, minmu_xvalues, minmu_yvalues,
    color=:orange,
    width=2,
    linestyle = :dash,
    label = "")
minmu = muexpvec[findall(x->x<150,minmass_array[3])][1]
minmu_yvalues = collect(0.8:0.1:1.3)
minmu_xvalues = repeat([minmu], outer=length(minmu_yvalues))
# Add the line to the plot
Plots.plot!(pslopehighz, minmu_xvalues, minmu_yvalues,
    color=:red,
    width=2,
    linestyle = :dash,
    label = "")
#because muexpvec is the exponent vector for mu = 10^i, use log10 below!
Plots.plot!(pslopehighz,log10.(expgainremainder_poor[:,1]),(expgainremainder_poor[:,2]),
    xlims=[muexpvec[1],muexpvec[end]+0.1],
    color=:black,
    label="Expected slope",
    width = 2,
    linestyle = :dot
    # markersize=3,
    )
display(pslopehighz)








current_gut_type = gut_types[1]
teeth = "all"
mass = 1000
mu = 10^-3.0
alpha = 4
# zetavec = collect(1.:0.5:3.)
reps = 50000
g_reps = Array{Float64}(undef,length(zetavec),reps)
c_reps = Array{Float64}(undef,length(zetavec),reps)
for i=1:length(zetavec)
    zeta = zetavec[i]
    for r=1:reps
        # Compute consumer and resource properties
        beta = bite_size_allo(mass)
        chewrate = chew_allo(mass, teeth)
        t_chewgram = 1 / chewrate
        tchew = t_chewgram * beta
        maxgut = gut_capacity_g(mass, current_gut_type)
        bcost_kJps, fcost_kJps = metabolic_cost(mass)
        velocity = find_velocity(mass)
        tmax_bout, _ = foragingtime(mass) .* (60 * 60)
        ndensity, n = indperarea(mass)
        width = reactionwidth(mass)
        height = reactionheight(mass)
        #adjacent competitor density (inds/m)
        na = n/(width*height);
                    
        m_res =  mu * (1 / beta) #* height #* width
        mprime = m_res / na
        alphaprime = alpha * na^(zeta - 2)
        
        # configurations = 200000  # note: adjust if needed
        gammadist = Gamma(alphaprime, mprime / alphaprime)
        
        # Compute daily foraging outcomes
        gains_daily, costs_daily, encounters_daily = dailyforage(gammadist, tchew, beta, maxgut, velocity, bcost_kJps, fcost_kJps, edensity, tmax_bout)

        g_reps[i,r] = gains_daily
        c_reps[i,r] = costs_daily
    end
end


exp_g = [mean(g_reps[i,:]) for i=1:length(zetavec)]
std_g = [std(g_reps[i,:]) for i=1:length(zetavec)]
exp_c = [mean(c_reps[i,:]) for i=1:length(zetavec)]
std_c = [std(c_reps[i,:]) for i=1:length(zetavec)]
gain_stats_df = DataFrame(
    zeta = zetavec,
    mean_gain = exp_g,
    std_gain = std_g,
    mean_cost = exp_c,
    std_cost = std_c
)
display(gain_stats_df)




# PLOT
gut_types_cap = ["Caecum", "Colon", "Non-rumen foregut", "Rumen foregut"]
x₀ = -2.5
colors = cgrad(:roma, 4, categorical = true)
#[:red, :blue, :green, :purple]
# 1) build your 2×2 canvas with NO per‐panel axis‐titles:
pratio = plot(
  layout                  = (2,2),
  frame                   = :box,
  foreground_color_legend = nothing,
  legend=:bottomright
)
zetavec2 = zetavec[1:3]
# 2) same double‐loop + your if‐logic for the legend
for (g, gut) in enumerate(gut_types_cap)
  for (j, z) in enumerate(zetavec2)
    # create masks for the two line‐styles
    mask1 = muexpvec .<= x₀
    mask2 = muexpvec .>= x₀
    col = colors[j]
    gammapos = findmin(abs.(slopegr[g,j,:] .- 1.17))[2]
    # — solid up to x₀
    plot!(
      pratio,
      muexpvec[mask1], gaincostratio[g, j, mask1],
      subplot  = g,
      ylabel   = ((g == 1 || g == 3) ? L"Reserve ratio, $\phi$"    : ""),
      xlabel   = ((g == 3 || g == 4) ? L"Richness, $\log(\mu)$" : ""),
      ylims = [1.0,1.5],
      linestyle= :solid,
      color     = col,
      label    = (g == 2 ? "ζ = $z" : false),
      width    = 2
    )
    # — dashed after x₀
    plot!(
      pratio,
      muexpvec[mask2], gaincostratio[g, j, mask2],
      subplot   = g,
      linestyle = :dash,
      color     = col,
      label     = false,
      width     = 2
    )
    scatter!(
        pratio,
        subplot  = g,
        [muexpvec[gammapos]],[gaincostratio[g,j,gammapos]],
        color = col,
        label = ""
    )
    # your annotation
    Plots.annotate!(
      pratio,
      subplot = g,
      -1, 1.45,
      text(string(gut), 8, halign = :right)
    )
  end
end
pratio

figfile = smartpath("rawfigures/fig_gainsremainder_gcratio.pdf")
Plots.savefig(pratio, figfile)

