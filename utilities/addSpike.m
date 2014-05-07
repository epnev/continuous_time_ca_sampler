function [newSpikeTrain, newCalcium, newLL] = addSpike(oldSpikeTrain,oldCalcium,oldLL,filter,tau,obsCalcium,timeToAdd,Dt,A)

    newSpikeTrain = [oldSpikeTrain timeToAdd]; %possibly inefficient, change if problematic (only likely to be a problem for large numbers of spikes)
    
    %use infinite precision to scale the precomputed FIR approximation to the calcium transient
    wk = A*exp((timeToAdd - Dt*ceil(timeToAdd/Dt))/tau);
    
    newCalcium = oldCalcium;
    tmp = 1 + (floor(timeToAdd):min((length(filter)+floor(timeToAdd)-1),length(newCalcium)-1));
    newCalcium(tmp) = newCalcium(tmp) + wk*filter(1:length(tmp));
    
    %if you really want to, ef*ef' could be precomputed and passed in
    relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
    
    newLL = oldLL - ( wk^2*norm(filter(1:length(tmp)))^2 - 2*relevantResidual*(wk*filter(1:length(tmp))'));