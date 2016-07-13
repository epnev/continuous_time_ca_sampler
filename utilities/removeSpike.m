function [newSpikeTrain, newCalcium, newLL] = removeSpike(oldSpikeTrain,oldCalcium,oldLL,filter,tau,obsCalcium,timeToRemove,indx,Dt,A) %#codegen
    
    tau_h = tau(1);
    tau_d = tau(2);
    
    ef_h = filter{1,1};
    ef_d = filter{1,2};
    ef_nh = filter{2,1};
    ef_nd = filter{2,2};
    
    newSpikeTrain = oldSpikeTrain;
    newSpikeTrain(indx) = [];
    
    %use infinite precision to scale the precomputed FIR approximation to the calcium transient    
    wk_h = A*exp((timeToRemove - Dt*ceil(timeToRemove/Dt))/tau_h);
    wk_d = A*exp((timeToRemove - Dt*ceil(timeToRemove/Dt))/tau_d);
    
    %%%%%%%%%%%%%%%%%
    %handle ef_h first
    newCalcium = oldCalcium;
    tmp = 1+ (floor(timeToRemove):min((length(ef_h)+floor(timeToRemove)-1),length(newCalcium)-1));
    wef_h = wk_h*ef_h(1:length(tmp));
    newCalcium(tmp) = newCalcium(tmp) - wef_h;

    relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
    relevantResidual(isnan(relevantResidual)) = 0;
    newLL = oldLL - ( wk_h^2*ef_nh(length(tmp)) + 2*relevantResidual*wef_h(:));
    oldCalcium = newCalcium;
    oldLL = newLL;
    %%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%
    %handle ef_d next
    newCalcium = oldCalcium;
    tmp = 1+ (floor(timeToRemove):min((length(ef_d)+floor(timeToRemove)-1),length(newCalcium)-1));
    wef_d = wk_d*ef_d(1:length(tmp));
    newCalcium(tmp) = newCalcium(tmp) - wef_d;

    relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
    relevantResidual(isnan(relevantResidual)) = 0;
    newLL = oldLL - ( wk_d^2*ef_nd(length(tmp)) + 2*relevantResidual*wef_d(:));
    %%%%%%%%%%%%%%%%