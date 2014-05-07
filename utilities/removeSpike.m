function [newSpikeTrain, newCalcium, newLL] = removeSpike(oldSpikeTrain,oldCalcium,oldLL,filter,tau,obsCalcium,timeToRemove,indx,Dt,A)
    
    newSpikeTrain = oldSpikeTrain;
    newSpikeTrain(indx) = [];
    
    %use infinite precision to scale the precomputed FIR approximation to the calcium transient
    wk = A*exp((timeToRemove - ceil(timeToRemove/Dt)*Dt)/tau);
    %wk = A*exp(-(timeToAdd - Dt*floor(timeToAdd/Dt))/tau);
    
    newCalcium = oldCalcium;
    tmp = 1+ (floor(timeToRemove):min((length(filter)+floor(timeToRemove)-1),length(newCalcium)-1));
    newCalcium(tmp) = newCalcium(tmp) - wk*filter(1:length(tmp));

    %if you really want to, ef*ef' could be precomputed and passed in
    relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
    newLL = oldLL - ( wk^2*norm(filter(1:length(tmp)))^2 + 2*relevantResidual*(wk*filter(1:length(tmp))'));
    
    