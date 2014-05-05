function [newSpikeTrain, newCalcium, newLL] = removeSpike(oldSpikeTrain,oldCalcium,oldLL,filter,tau,obsCalcium,timeToRemove,indx,Dt,A)
    
    newSpikeTrain = oldSpikeTrain;
    newSpikeTrain(indx) = [];
    
    %use infinite precision to scale the precomputed FIR approximation to the calcium transient
    wk = A*exp(-(timeToRemove - floor(timeToRemove/Dt)*Dt)/tau);
    
    newCalcium = oldCalcium;
    tmp = 1+ (floor(timeToRemove):min((length(filter)+floor(timeToRemove)-1),length(newCalcium)-1));
    newCalcium(tmp) = newCalcium(tmp) - wk*filter(1:length(tmp));

    %if you really want to, ef*ef' could be precomputed and passed in
    relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
    newLL = oldLL - ( wk^2*(filter*filter') + 2*relevantResidual*(wk*filter(1:length(tmp))'));
    
    