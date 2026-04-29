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

#NOTE: modified for SA
muexpvec = collect(-3.5:0.02:-1); #collect(-3.5:0.02:-1);
l_muexpvec = length(muexpvec);
#NOTE: modified for SA
zetavec = [1.0]; #[1.0, 1.5, 2.0, 2.15];
l_zetavec = length(zetavec);

#NOTE: modified for SA
reps = 500; #reps = 1000

# Pre-allocate an array to store slopes for each gut type and each muexp value:
remainder_mu_slope_all = Array{Float64}(undef, n_gut, l_muexpvec);

gains_array = Array{Array{Float64}}(undef,n_gut);
costs_array = Array{Array{Float64}}(undef,n_gut);
encounters_array = Array{Array{Float64}}(undef,n_gut);
cons_maxgut_array = Array{Array{Float64}}(undef,n_gut);
gaindiff_array = Array{Array{Float64}}(undef,n_gut);
gainremainder_array = Array{Array{Float64}}(undef,n_gut);

# MODIFICATION FOR QUALITY V BIOMASS SENSITVITY ANALYSIS
# proportion digestible by herbivore
thetavec =  Array{Float64}(undef, l_muexpvec);
# Baseline: set to 1
[thetavec[i] = 1.0 for i in 1:l_muexpvec];
# Include vegetative quality vs. vegetative biomass relationship
# theta0 = 67.968/100;
# theta1 = (0.0246/100); #*5000;
# [thetavec[i] = theta0 - theta1*(10^muexpvec[i]) for i in 1:l_muexpvec];

boutexp = collect(-0.5:0.1:0.5)
gammacurve = Array{Float64}(undef,length(boutexp),l_muexpvec,2)

for s in 1:length(boutexp)

    # Loop over each gut type
    for g in 1:1
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

                        #Proportion digestible
                        theta = thetavec[i];
                        
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
                        tmax_bout, _ = foragingtime_alt(mass,boutexp[s]) .* (60 * 60)
                        
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
                        #NOTE: modified gains for SA
                        gains[r,m,i,j] = theta*gains_daily
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
    for g in 1:1
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
    for g = 1:1
        minmuany = muexpvec[findall(!isnan,minmass_array[1])[1]]
        minmuany_array[g] = minmuany
    end

    gammacurve[s,:,:] = [muexpvec remainder_mu_slope_all[1, :]];

end



#LOAD ANALYTICAL SOLUTION
filename = smartpath("data/expgainremainderslope.CSV")
expgainremainder_poor = CSV.read(filename, DataFrame)
expgainremainder_poor = Array(expgainremainder_poor)


###########################################################
# GENERATE FIGURE 5
###########################################################

# gut_types_cap = ["Caecum","Colon","Non-rumen foregut","Rumen foregut"]
# gut_type_maxgutslope = [0.860,0.919,0.881,0.897];
pslope = plot(
    xlabel = L"Richness, $\log_{10}(\mu)$",
    ylabel = L"Gain surplus exponent, $g_S \propto M^\gamma$",
    ylims = (0.8, 1.3),
    legend = :topright,
    legend_title = "Exp. Diff (%)",
    frame = :box,
    size=(500,350),
    foreground_color_legend = nothing
              )
colors = palette(:matter,11) #, n_gut)  # Get exactly n_gut colors
for s in 1:length(boutexp)
    col = colors[s]
    # gutslope = gut_type_maxgutslope[g]
    # x = log.(massvec)
    # y = log.(vec(mean(costs_array[g],dims=1)[:,:,end,1]))
    # df = DataFrame(x = x, y = y)
    # model = lm(@formula(y ~ x), df)
    # intc, _ = coef(model)
    # x = log.(massvec)
    # y = log.(vec(mean(gains_array[g],dims=1)[:,:,end,1]))
    # df = DataFrame(x = x, y = y)
    # model = lm(@formula(y ~ x), df)
    # intg, g_slope = coef(model)
    # int = exp(intg)/exp(intc)
    # c_slope =0.75;
    # # exp_slope = c_slope + ((int*(gutslope-c_slope))*massvec[end]^(gutslope-c_slope))/(int*massvec[end]^(gutslope-c_slope) - 1)
    # # From Mathematica Solution
    # M = 100000; #Asymptotic Mass
    # # exp_slope2 = g_slope + (intc*(c_slope - g_slope)*M^c_slope)/(intc*M^c_slope - intg*M^g_slope)
    # exp_slope2 = gutslope + (intc*(c_slope - gutslope)*M^c_slope)/(intc*M^c_slope - intg*M^gutslope)
    if isodd(s) || boutexp[s] == 0.0
        plot!(gammacurve[s,:,1], gammacurve[s,:,2], lw=2,color=col,label="$(boutexp[s])")
    else
         plot!(gammacurve[s,:,1], gammacurve[s,:,2], lw=2,color=col,label="")
    end
    # scatter!(pslope, [muexpvec[end]], [exp_slope2],markersize=5,color=col,label="")
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

display(pslope)

# SAVE FIGURE 4
figfile = smartpath("rawfigures/fig_gainsremainder_slopev2_tbouterror.pdf")
Plots.savefig(pslope, figfile)









































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

