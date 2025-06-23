# most but not all of these are currently in use.  

function find_velocity(mass)
    # mass in kg
    # from kramer 2010 "allometric scaling of resource acquisition"
    #Consumer Velocity (meters/second)
    velocity = (0.5 * mass^0.13);
    return velocity
end


function bite_size_allo(mass)
    # from shipley 94 "the scaling of intake rate"
    # mass in kg
    # bite size in g
    #if plant == "browse"
    #    bite_size = (0.057* (mass)^0.63); # [g]
    #elseif plant == "graze"
    #    bite_size = (0.026 * (mass)^0.59); # [g]
    #end

    # Why elephants have trunks Pretorius 2015
    bite_size = (0.002 * (mass)^0.969); # [g]
    return bite_size

end



function number_of_chews(mass)
    # shipley 94
   #chew/g (processed to mean particle size (allo) and swallowed)
   # mass in kg
    
   #uses the correction factor reported in Shipley
   chews_per_gram = 343.71* mass^(-0.83); 
    

   return chews_per_gram

end



function chew_rate_allo(mass, teeth)
   # from "dental functional morphology predicts scaling"
   # mass in kg
   # duration in ms
   if teeth == "bunodont"
       chewing_cycle_duration = (228.0* (mass)^0.246) / 1000; # 2.358 [ms -> s]

   elseif teeth == "acute/obtuse lophs"
       chewing_cycle_duration = (299.2 * (mass)^0.173) / 1000; # 2.476 [ms -> s]

   elseif teeth == "lophs and non-flat"
       chewing_cycle_duration = (320.6 * (mass)^0.154) / 1000; # 2.506 [ms -> s]

   elseif teeth == "lophs and flat"
       chewing_cycle_duration = (262.4* (mass)^0.207) / 1000; # 2.419 [ms -> s]0

   elseif teeth == "all"

        chewing_cycle_duration = (266.07* (mass)^0.196) / 1000; # 2.419 [ms -> s]0

   end

   return 1 / (chewing_cycle_duration )  #[s/chew -> chews/s]

end


function chew_allo(mass, teeth)
    # allometric function for mouth->gut rate
   chew_rate = chew_rate_allo(mass, teeth); # chews/s
   chews_per_gram = number_of_chews(mass)   # chew/g
   chew = chew_rate / chews_per_gram        # chew/s / chew/g -> g/s
   return chew
end



function gut_capacity_g(mass, gut_type)
    # from Muller 2013 - Assessing the Jarmain-Bell Principle
    # dry mass of guts in kg
    # bm in kg

    if gut_type == "caecum"
        capacity = (0.025*(mass)^0.860);

    elseif gut_type == "colon"
        capacity = (0.029 * (mass)^0.919);

    elseif gut_type == "non-rumen foregut"
        capacity = (0.030 * (mass)^0.881);

    elseif gut_type == "rumen foregut"
        capacity = (0.041 * (mass)^0.897);

    end

    return capacity * 1000 #[kg -> g]

end



function metabolic_cost(mass)
    # takes a terrestrial mammal body mass (kg), returns storage masses and metabolic rates
    #Convert to grams
    mass_g = mass*1000;

    b0_bmr = 0.018; #watts g^-0.75
    b0_fmr = 0.047; #watts g^-0.75
    bcost_watts = (b0_bmr*(mass_g^0.75)); #Cost in Watts
    fcost_watts = (b0_fmr*(mass_g^0.75)); #Cost in Watts
    
    #Watt = 0.001 kJ/s
    kJps = 0.001;

    bcost_kJps = bcost_watts*kJps; #kJ per s
    fcost_kJps = fcost_watts*kJps; #kJ per s

    return bcost_kJps, fcost_kJps

end

function indperarea(mass)
    #Enter mass in kg
    #Convert to grams
    # massg = mass*1000;
    # popdensity = (0.0116)*massg^-0.776; #inds/area (from Damuth)
    popdensity = (5.45*10^-5)*mass^-0.776; #inds/area (from Damuth)

    foragebout_hrs, _ = foragingtime(mass); #hrs
    foragebout_s = foragebout_hrs *60*60;
    foragevelocity = find_velocity(mass); #m/s
    areaforaged = pi*((foragevelocity * foragebout_s)/2)^2 #m/s * s * m

    corridorwidth = 2
    corridorareaforaged = foragevelocity * foragebout_s * corridorwidth #corridor 1 m

    #Or use homerange scaling
    #Owen-Smith 1988 pg 96
    HR = 13500*mass^1.25 #m^2

    popinarea = 1. + popdensity*(corridorareaforaged)

    return popdensity, popinarea
end




# function maxfatstorage(mass,edensity_fat)
#     mass_g = mass * 1000; #[kg]->[g]

#     # from Yeakel et al 2018 Dynamics of starvation and recovery
#     #mass at which you die, no fat or muscle
#     fat_mass =  0.02*mass_g^1.19;           #[g]

#     storage_kg = fat_mass/1000;

#     #Joules per gram
#     # joules_per_gram = 20000; #varies between 7000 to 36000 [int] [J/g]
#     # kjoules_per_gram = joules_per_gram / 1000;   # [int/int=int] [kJ/g]
#     kjoules_per_gram = edensity_fat; #kJ/g

#     storage_kj = (fat_mass) * kjoules_per_gram;

#     return storage_kj, storage_kg
# end


function fatrecoveryrate(mass,start_perc,end_perc)
    #From Yeakel et al. 2018
    #mass should be in grams

    mass_g = mass*1000;

    B0 = 0.047; #Wg-3/4
    Emprime = 7000; #J/g
    aprime = B0/Emprime;
    eta = 3/4;

    epsilon_sig = start_perc;
    epsilon_lam = end_perc;
    t_recover = log(((1-(epsilon_sig * epsilon_lam)^(1-eta))/(1-epsilon_lam^(1-eta))))*(mass_g^(1-eta)/(aprime*(1-eta))) #seconds

    fatsynthrate = 1/t_recover; #1/sec

    return fatsynthrate
end

function expectedlifetime(mass)
    #Calder 1984;
    #mass in kg
    #expected lifetime in years
    explifetime = 5.08*mass^0.35;
    return explifetime
end

function foragingtime(mass)
    #Owen Smith 1988 book (data grabbed using PlotDigitizer in /data/)
    #mass in kg
    #foraging time in % of 24 hours

    forageperc = 21.0885*mass^0.124
    forageperc_upper95 = 27*mass^0.17

    #translate to hours
    foragehrs = (forageperc/100)*24
    foragehrs_upper95 = (forageperc_upper95/100)*24

    return foragehrs, foragehrs_upper95 # hours
end

function gut_turnover_time(mass,gut_type,edensity)
    # mass in kg
    # edensity of food in kJ/g

    # maxgut
    maxgut = gut_capacity_g(mass, gut_type) * edensity; #kJ

    #metabolic demand (rate) kJ per day
    #From Calder book pg 126
    dailymetrate = 504*mass^0.76;

    turnovertime = maxgut / dailymetrate; # days
    turnovertime_hrs = turnovertime * 24; # days * 24 hours/day = hrs

    return turnovertime_hrs

end

function dailyfoodintake(mass)
    intake_kJday = 971*mass^0.73
    return intake_kJday
end


#reaction distance
function reactionwidth(mass)
    #mass in kg
    #distance in meters
    #see 'fit' with Pawar data in analysis
    #Note it is not a fit, but anchoring the intercept 

    width = 2 * 5 * mass^(1/3) #meters

    return width

end


# reaction height
function reactionheight(mass)
    #mass in KG
    #distance in meters ~ shoulder height
    #From Larramendi 2016
    ra = 0.1501*mass^0.357
end

function carryingcapacity_costs(mass)

    #birth mass in Kg
    m0 = 0.5581*mass^0.92
    m0g = m0*1000

    B0 = 8.36 #W*kg^-0.75
    B0g = 4.7*10^(-2)

    #cost of tissue synthesis (kJ/kg)
    Em = 5774 * (1/1000) * 1000 #J/g * kJ/J * g/kg

    a = B0/Em;
    ag = B0g/Em

    #lifetime in days
    #time to 1/2 cohort survival - calder book pg 317
    t_lifetime_d = 365*5.44*mass^0.32

    #this uses gram mass
    massg = 1000*mass
    t_lifetime_r =  (-log((1 - 0.99^(1/4))/(1 - (m0g/massg)^(1/4)))*(4*massg^(1/4))/ag)* (1/(60*60*24)) # s * (years/s)

    #proportion of adult mass that ends Juvenile period
    ej = 0.75

    #Costs of growth and replacement in kJ
    E_kcosts = Em*(mass - m0 + 2*ej*mass)

    return E_kcosts, t_lifetime_r, t_lifetime_d
end



# # Bite chew time allometry
# massvec = [10^i for i=collect(0:0.1:5)];
# teeth = "bunodont";
# betavec = bite_size_allo.(massvec);
# chewratevec = chew_allo.(massvec,teeth);
# tchewgramvec = 1 ./ chewratevec;
# tchewvec = tchewgramvec .* betavec
# namespace = smartpath("figures/tchew_allometry.pdf")
# R"""
# pdf($namespace,width=4,height=4)
# plot($massvec,$tchewvec,type='l',lwd=2,xlab='Mass (kg)',ylab='Bite chew time (s)')
# dev.off()
# """
