function [samples, ci]  = get_next_spikes(curr_spikes,curr_calcium,calciumSignal,ef,tau,calciumNoiseVar, lam, proposalVar, add_move, Dt, A)
    
    %addMoves, dropMoves, and timeMoves give acceptance probabilities for each subclass of move
    %the samples will be a cell array of lists of spike times - the spike times won't be sorted but this shouldn't be a problem.
    
    %noise level here matters for the proposal distribution (how much it 
    %should trust data for proposals vs how much it should mix based on uniform prior)
    %this is accounted for by calciumNoiseVar
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% initialize some parameters
    T = length(calciumSignal); %for all of this, units are bins and spiketrains go from 0 to T where T is number of bins
    ff = ~isnan(calciumSignal); % entries with data
    
%   nsweeps = 1e3; %number of sweeps.
    nsweeps = 1;
    samples = cell(nsweeps);
    
    %% start with initial spiketrain and initial predicted calcium 
    si = curr_spikes; %initial set of spike times has no spikes - this will not be sorted but that shouldn't be a problem
    ci = curr_calcium; %initial calcium is set to baseline 
    
    N = length(si); %number of spikes in spiketrain
        
    %initial logC - compute likelihood initially completely - updates to likelihood will be local
    logC = -norm(ci(ff)-calciumSignal(ff))^2; 
    
    %m = p_spike*Dt*T;
    
    %flag for uniform vs likelihood proposal (if using likelihood proposal, then time shifts are pure Gibbs)
    %this really should be split into four cases, 
    % 1) RW for time shifts with uniform add/drop
    % 2) RW for time shifts with likelihood proposal add/drop
    % 3) Gibbs for time shifts with uniform add/drop
    % 4) Gibbs for time shifts with likelihood proposal add/drop
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% loop over sweeps to generate samples
    addMoves = [0 0]; %first elem is number successful, second is number total
    dropMoves = [0 0];
    timeMoves = [0 0];
    time_move = 0;
    time_add = 0;
    for i = 1:nsweeps
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% loop over spikes, perform spike time move (could parallelize here for non-interacting spikes, i.e. spikes that are far enough away)
        
        %move
        for ni = 1:N %possibly go through in a random order (if you like)
            
            tmpi = si(ni);
            tmpi_ = si(ni)+(proposalVar*randn); %with bouncing off edges
            if tmpi_<0
                tmpi_ = -(tmpi_);
            elseif tmpi_>T
                tmpi_ = T-(tmpi_-T);
            end           

            %set si_ to set of spikes with the move and ci_ to adjusted calcium and update logC_ to adjusted
%           [si_, ci_, logC_] = removeSpike(si,ci,logC,ef,tau,calciumSignal,tmpi,ni,Dt,A);
%           [si_, ci_, logC_] = addSpike(si_,ci_,logC_,ef,tau,calciumSignal,tmpi_,ni,Dt,A);     
            [si_, ci_, logC_] = replaceSpike(si,ci,logC,ef,tau,calciumSignal,tmpi,ni,tmpi_,Dt,A);
            
            %accept or reject
            ratio = exp((logC_-logC)/(2*calciumNoiseVar)*lam(tmpi)/lam(tmpi_));
            if ratio>1 %accept
                si = si_;
                ci = ci_;
                logC = logC_;
                timeMoves = timeMoves + [1 1];
            elseif rand<ratio %accept
                si = si_;
                ci = ci_;
                logC = logC_;
                timeMoves = timeMoves + [1 1];
            else
                %reject - do nothing
                timeMoves = timeMoves + [0 1];
            end

            %visualize
            if (0)
            figure(46)
            subplot(211)
            stem(si,ones(1,length(si)))
            ylim([-.5 2])
            subplot(212)
            plot(calciumSignal)
            hold on
            plot(ci,'r');
            hold off
            end
%             pause
            
        end
        N = length(si);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% loop over add/drop a few times
        %define insertion proposal distribution as the likelihood function
        %define removal proposal distribution as uniform over spikes
        %perhaps better is to choose smarter removals.
        for ii = 1:add_move 
            %% add
            %propose a uniform add
            tmpi = T*Dt*rand;         
            [si_, ci_, logC_] = addSpike(si,ci,logC,ef,tau,calciumSignal,tmpi,length(si_)+1,Dt,A);
        
            %forward probability
            fprob = 1/(T*Dt);
            
            %reverse (remove at that spot) probability
            rprob = 1/(N+1);
            
            %accept or reject
            ratio = exp((1/(2*calciumNoiseVar))*(logC_ - logC))*(rprob/fprob)*lam(tmpi); %posterior times reverse prob/forward prob
            if ratio>1 %accept
                si = si_;
                ci = ci_;
                logC = logC_;
                addMoves = addMoves + [1 1];
            elseif rand<ratio %accept
                si = si_;
                ci = ci_;
                logC = logC_;
                addMoves = addMoves + [1 1];
            else
                %reject - do nothing
                addMoves = addMoves + [0 1];
            end
            N = length(si);

            %% delete
            if N>0                
                %propose a uniform removal
                tmpi = randi(N);
                [si_, ci_, logC_] = removeSpike(si,ci,logC,ef,tau,calciumSignal,si(tmpi),tmpi,Dt,A);

                %reverse probability
                rprob = 1/(T*Dt);

                %compute forward prob
                fprob = 1/N;
                
                %accept or reject
                ratio = exp((1/(2*calciumNoiseVar))*(logC_ - logC))*(rprob/fprob)*(1/lam(si(tmpi))); %posterior times reverse prob/forward prob
               
                if ratio>1 %accept
                    si = si_;
                    ci = ci_;
                    logC = logC_;
                    dropMoves = dropMoves + [1 1];
                elseif rand<ratio %accept
                    si = si_;
                    ci = ci_;
                    logC = logC_;
                    dropMoves = dropMoves + [1 1];
                else
                    %reject - do nothing
                    dropMoves = dropMoves + [0 1];
                end
                N = length(si);
            
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
%        N_sto = [N_sto N];

        samples = si;
        
        %store overall logliklihood as well
%        objective = [objective logC];
        
        if (0)
        figure(48)
        subplot(211)
        plot(N_sto)
        subplot(212)
        plot(objective)
        
        [addMoves(1),addMoves(2), dropMoves(1),dropMoves(2)]
        end  
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%