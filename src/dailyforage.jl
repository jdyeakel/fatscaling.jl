function dailyforage(gammadist,tchew, beta, maxgut, velocity, bcost_kJps, fcost_kJps, edensity, tmax_bout)

    #Initialize
    t_travel = 0.0;
    t_chew = 0.0;
    bites = 0;
    gut = 0.0
    t=0.0;

    #Simulate daily returns and costs
    while t < tmax_bout && gut < maxgut
        
        # tic += 1

        # Draw distance to next resource
        # First draw encounter rate
        rg = rand(gammadist);
        distance_to_resource = rand(Exponential(1.0/rg));
            
        #The forager will move towards the closest resource
        ttravel = distance_to_resource/velocity; # distance / velcity = time
        
        t += ttravel; # time
        t_travel += ttravel; # time
        
        # If the forager successfully reaches the resource
        # Each resource is equal to 1 mouthful
        if (t + ttravel) < tmax_bout
            # number_of_successes += 1;
            
            # Time passes while chewing
            t += tchew; #time
            t_chew += tchew; #time

            # Pass Mouth-Unit to Gut Total (boundary conditions in across-day sim)
            # resgain is kJ per volume (mouth volume)
            gut += beta; #grams/bite
            bites += 1;
            
        end

    end


    #kJ gains
    gains = gut * edensity; #grams * kJ/gram = kJ returns

    #Non-forage time
    #NOTE: changed tmax_bout to t, because you can stop if your gut fills, such that t < tmax_bout, in which case, the remainder is rest.
    #11/15/2024
    t_rest = (24*60*60) - t; #tmax_bout;

    #kJ costs
    costs = fcost_kJps*t_travel + bcost_kJps*(t_chew + t_rest);

    return gains, costs, bites
    
end