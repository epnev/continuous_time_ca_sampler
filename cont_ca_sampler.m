function SAMPLES = cont_ca_sampler(Y,params)

% Continuous time sampler
% Y                     data (normalized in [0,1])
% P                     intialization parameters (discrete time constant P.g required)
% params                parameters structure
% params.g              discrete time constant(s) (estimated if not provided)
% params.sn             initializer for noise (estimated if not provided)
% params.b              initializer for baseline (estimated if not provided)
% params.c1             initializer for initial concentration (estimated if not provided)
% params.Nsamples       number of samples after burn in (default 500)
% params.B              number of burn in samples (default 200)
% params.marg           flag for marginalized sampler (default 0)
% params.upd_gam        flag for updating gamma (default 0)
% params.gam_step       number of samples after which gamma is updated (default 50)
% params.std_move       standard deviation of shifting kernel (default 3*Dt)
% params.add_move       number of add moves per iteration (default T/100)
% params.init           initial sample 
% params.f              imaging rate (default 1)
% params.p              order of AR model (p == 1 or p == 2, default 1)

% output struct SAMPLES
% spikes                T x Nsamples matrix with spikes samples
% bp                    Nsamples x 1 vector with samples for spiking prior probability
% Am                    Nsamples x 1 vector with samples for spike amplitude

% If marginalized sampler is used
% Cb                    posterior mean and sd for baseline
% Cin                   posterior mean and sd for initial condition
% else
% Cb                    Nsamples x 1 vector with samples for baseline
% Cin                   Nsamples x 1 vector with samples for initial concentration
% sn                    Nsamples x 1 vector with samples for noise variance

% If gamma is updated
% g                     Nsamples x 1 vector with the gamma updates

% Author: Eftychios A. Pnevmatikakis and Josh Merel

Y = Y(:);
T = length(Y);

% define default parameters
defparams.g = [];
defparams.sn = [];
defparams.b = [];
defparams.c1 = [];
defparams.Nsamples = 100;
defparams.B = 100;
defparams.marg = 0;
defparams.upd_gam = 1; %changed from 1
defparams.gam_step = 1; %changed from 50
defparams.std_move = 3;
defparams.add_move = ceil(T/100);
defparams.init = [];
defparams.f = 1;
defparams.p = 1;

if nargin < 2
    params = defparams;
else
    if ~isfield(params,'g'); params.g = defparams.g; end
    if ~isfield(params,'sn'); params.sn = defparams.sn; end
    if ~isfield(params,'b'); params.b = defparams.b; end
    if ~isfield(params,'c1'); params.c1 = defparams.c1; end
    if ~isfield(params,'Nsamples'); params.Nsamples = defparams.Nsamples; end
    if ~isfield(params,'B'); params.B = defparams.B; end
    if ~isfield(params,'marg'); params.marg = defparams.marg; end
    if ~isfield(params,'upd_gam'); params.upd_gam = defparams.upd_gam; end
    if ~isfield(params,'gam_step'); params.gam_step = defparams.gam_step; end
    if ~isfield(params,'std_move'); params.std_move = defparams.std_move; end
    if ~isfield(params,'add_move'); params.add_move = defparams.add_move; end
    if ~isfield(params,'init'); params.init = defparams.init; end
    if ~isfield(params,'f'); params.f = defparams.f; end
    if ~isfield(params,'p'); params.p = defparams.p; end
end
       
Dt = 1;                                     % length of time bin

marg_flag = params.marg;
gam_flag = params.upd_gam;
gam_step = params.gam_step;
std_move = params.std_move;
add_move = params.add_move;

if gam_flag
    options = optimset('GradObj','On','Display','Off','Algorithm','interior-point','TolX',1e-6);
end

if isempty(params.g)
    p = params.p;
else
    p = length(params.g);                       % order of autoregressive process
end

if isempty(params.init)
   fprintf('Initializing using noise constrained FOOPSI...  ');
   params.init = get_initial_sample(Y,params);
   fprintf('done. \n');
end 
SAM = params.init;

g = SAM.g(:)';
if p == 1
    tau_2 = -Dt/log(g);                         % continuous time constant
    tau_1 = 0;
    tau = [tau_1,tau_2];
    G1 = speye(T);
    G2 = spdiags(ones(T,1)*[-g,1],[-1:0],T,T);
elseif p == 2
    gr = roots([1,-g]);
    p1_continuous = log(min(gr))/Dt; 
    p2_continuous = log(max(gr))/Dt;
    tau_1 = -1/p1_continuous;                   %tau h - smaller (tau_d * tau_r)/(tau_d + tau_r)
    tau_2 = -1/p2_continuous;                   %tau decay - larger
    tau = [tau_1,tau_2];                        %tau_h , tau_d
    G1 = spdiags(ones(T,1)*[-min(gr),1],[-1:0],T,T);
    G2 = spdiags(ones(T,1)*[-max(gr),1],[-1:0],T,T);
else
    error('This order of the AR process is currently not supported');
end

sg = SAM.sg;


SAM = params.init;    
spiketimes_ = SAM.spiketimes_;
lam_ = SAM.lam_;
A_ = SAM.A_;
b_ = SAM.b_;
C_in = SAM.C_in;
    
s_1 = sparse(ceil(spiketimes_/Dt),1,exp((spiketimes_ - Dt*ceil(spiketimes_/Dt))/tau_1),T,1);  
s_2 = sparse(ceil(spiketimes_/Dt),1,exp((spiketimes_ - Dt*ceil(spiketimes_/Dt))/tau_2),T,1);  

prec = 1e-2;     % precision
if p == 1
    h = exp(-(Dt)/tau_2);
    ge = g.^((0:T-1)'); 
    e_support = find(abs(ge)<prec,1);
    ge = sparse(ge);
    ef_d = ge(1:e_support)';
    ef = [{[0,0], ef_d};{[0,0], cumsum(ef_d.^2)}];
else
    if g(2) == 0;
        tau_1 = 1e-3;
    end
    h = exp(-(Dt)/tau_2) - exp(-(Dt)/tau_1);
    ef_d = exp(-(0:T)/tau_2)/h;
    e_support = find(abs(ef_d)<prec,1);
    ef_d = ef_d(1:e_support);

    ef_h = -exp(-(0:T)/tau_1)/h;
    e_support = find(abs(ef_h)<prec,1);
    ef_h = ef_h(1:e_support);
    ef = [{ef_h ef_d};{cumsum(ef_h.^2) cumsum(ef_d.^2)}];
    ge = G2\[1;zeros(T-1,1)];
end
B = params.B;
N = params.Nsamples + B;
Gs = (-G1\s_1(:)+G2\s_2(:))/h;


ss = cell(N,1); 
lam = zeros(N,1);
Am = zeros(N,1);
ns = zeros(N,1);
Gam = zeros(N,p);
if ~marg_flag
    Cb = zeros(N,1);
    Cin = zeros(N,1);
    SG = zeros(N,1);
end

Sp = .1*range(Y)*eye(3);                                    % prior covariance
Ld = inv(Sp);
lb = [0.1*range(Y),0.002*range(Y),0]';      % lower bound for [A,Cb,Cin]

A_ = max(A_,1.1*lb(1));

mu = [A_;b_;C_in];                % prior mean 
Ns = 15;                          % Number of HMC samples

Ym = Y - ones(T,1)*mu(2) - ge*mu(3);

mub = zeros(1+p,1);
Sigb = zeros(1+p,1+p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra tau-related params
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau1_std = .1;
tau2_std = 1;
tauMoves = [0 0];
tau_min = 0;
tau_max = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:N
    if gam_flag
        Gam(i,:) = g;
    end
    sg_ = sg;
    rate = @(t) lambda_rate(t,lam_);
    [spiketimes, ~]  = get_next_spikes(spiketimes_(:)',A_*Gs',Ym',ef,tau,sg_^2, rate, std_move, add_move, Dt, A_);
    spiketimes_ = spiketimes;
    spiketimes(spiketimes<0) = -spiketimes(spiketimes<0);
    spiketimes(spiketimes>T*Dt) = 2*T*Dt - spiketimes(spiketimes>T*Dt); 
    trunc_spikes = ceil(spiketimes/Dt);
    trunc_spikes(trunc_spikes == 0) = 1;
    %s_ = sparse(1,trunc_spikes,exp((spiketimes - Dt*trunc_spikes)/tau(2)) - exp((spiketimes - Dt*trunc_spikes)/tau(1)),1,T);
    s_1 =   sparse(trunc_spikes,1,exp((spiketimes_ - Dt*trunc_spikes)/tau(1)),T,1);  
    s_2 =   sparse(trunc_spikes,1,exp((spiketimes_ - Dt*trunc_spikes)/tau(2)),T,1);  
    Gs = (-G1\s_1(:)+G2\s_2(:))/h;
    ss{i} = spiketimes;
    nsp = length(spiketimes);
    ns(i) = nsp;
    lam(i) = nsp/(T*Dt);
    lam_ = lam(i);
    
    %AM = [G\s_(:),ones(T,1),ge];
    AM = [Gs,ones(T,1),ge];
    L = inv(Ld + AM'*AM/sg^2);
    mu_post = (Ld + AM'*AM/sg^2)\(AM'*Y/sg^2 + Sp\mu);
    if ~marg_flag
        x_in = [A_;b_;C_in];
        if any(x_in < lb)
            x_in = max(x_in,1.1*lb);
        end
        [temp,~] = HMC_exact2(eye(3), -lb, L, mu_post, 1, Ns, x_in);
        Am(i) = temp(1,Ns);
        Cb(i) = temp(2,Ns);
        Cin(i) = temp(3,Ns)';
        A_ = Am(i);
        b_ = Cb(i);
        C_in = Cin(i);

        Ym   = Y - b_ - ge*C_in;
        res   = Ym - A_*Gs;
        sg   = 1./sqrt(gamrnd(1+T/2,1/(0.1 + sum((res.^2)/2))));
        SG(i) = sg;
    else
        repeat = 1;
        while repeat
            A_ = mu_post(1) + sqrt(L(1,1))*randn;
            repeat = (A_<0);
        end                
        Am(i) = A_;
        if i > B
           mub = mub + mu_post(2+(0:p));
           Sigb = Sigb + L(2+(0:p),2+(0:p));
        end
    end
    if gam_flag
        if (i > B) && mod(i-B,gam_step) == 0  % update gamma
            fprintf('updating time constant.. ')
            ss_ = zeros(1,T);
            for rep = i-gam_step+1:i
                trunc_spikes = ceil(ss{rep}/Dt);            
                ss_ = ss_ + Am(rep)*sparse(1,trunc_spikes,exp((ss{rep} - Dt*trunc_spikes)/tau(2)),1,T);
            end
            ss_ = ss_/gam_step;        
            if ~marg_flag
                y_res = Y - mean(Cb(i-gam_step+1:i));
                ss_(1) = ss_(1) + mean(Cin(i-gam_step+1:i));
            else
                y_res = Y - mub(1)/(i-B);
                ss_(1) = ss_(1) + mub(2)/(i-B);
            end
        
            %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
            %update continuous taus
            %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
            
            calciumNoiseVar = sg^2;
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % first update tau(1)
            %%%%%%%%%%%%%%%%%%%%%%%
            if p >= 2
            %initial logC
            logC = -(A_*Gs'-Ym')*(A_*Gs'-Ym')'; 
            
            tau_ = tau;
            
            tau_(1) = tau_(1)+(tau1_std*randn); %with bouncing off min and max
            while tau_(1)>tau(2) || tau_(1)<tau_min
                if tau_(1)<tau_min
                    tau_(1) = tau_min+(tau_min-tau_(1));
                elseif tau_(1)>tau(2)
                    tau_(1) = tau(2)-(tau_(1)-tau(2));
                end
            end 
        
            gr_ = exp(Dt*(-1./tau_));
           
            h_ = exp(-(Dt)/tau_(2)) - exp(-(Dt)/tau_(1));
            G1_ = spdiags(ones(T,1)*[-min(gr_),1],[-1:0],T,T);
            G2_ = spdiags(ones(T,1)*[-max(gr_),1],[-1:0],T,T);
            s_1_ = sparse(ceil(spiketimes_/Dt),1,exp((spiketimes_ - Dt*ceil(spiketimes_/Dt))/tau_(1)),T,1);  
            s_2_ = sparse(ceil(spiketimes_/Dt),1,exp((spiketimes_ - Dt*ceil(spiketimes_/Dt))/tau_(2)),T,1);  
            Gs_ = (-G1_\s_1_(:)+G2_\s_2_(:))/h_;
                        
%             % alternative way of normalizing filters?
%             ef_d = exp(-(0:T)/tau(2));
%             ef_h = -exp(-(0:T)/tau(1));
%             
%             %compute maximum:
%             to = (tau(1)*tau(2))/(tau(2)-tau(1))*log(tau(2)/tau(1)); %time of maximum
%             max_val = exp(-to/tau(2))-exp(-to/tau(1)); %maximum
%             ef = {ef_h/max_val ef_d/max_val};
            
            logC_ = -(A_*Gs_'-Ym')*(A_*Gs_'-Ym')';
                
            %accept or reject
            prior_ratio = 1;
    %         prior_ratio = gampdf(tau_(2),12,1)/gampdf(tau(2),12,1);
            ratio = exp(sum(sum((1./(2*calciumNoiseVar)).*(logC_-logC))))*prior_ratio;
            if ratio>1 %accept
                tau = tau_;
                h = h_; G1 = G1_; G2 = G2_; s_1 = s_1_; s_2 = s_2_; Gs = Gs_;
                tauMoves = tauMoves + [1 1];
            elseif rand<ratio %accept
                tau = tau_;
                h = h_; G1 = G1_; G2 = G2_; s_1 = s_1_; s_2 = s_2_; Gs = Gs_;
                tauMoves = tauMoves + [1 1];
            else
                tauMoves = tauMoves + [0 1];
            end
            end
            %%%%%%%%%%%%%%%%%%%%%%%
            % next update tau(2)
            %%%%%%%%%%%%%%%%%%%%%%%
            
            %initial logC
            logC = -(A_*Gs'-Ym')*(A_*Gs'-Ym')'; 
    
            tau_ = tau;

            tau_(2) = tau_(2)+(tau2_std*randn);
            while tau_(2)>tau_max || tau_(2)<tau_(1)
                if tau_(2)<tau_(1)
                    tau_(2) = tau_(1)+(tau_(1)-tau_(2));
                elseif tau_(2)>tau_max
                    tau_(2) = tau_max-(tau_(2)-tau_max);
                end
            end  
        
            gr_ = exp(Dt*(-1./tau_));
           
            h_ = exp(-(Dt)/tau_(2)) - exp(-(Dt)/tau_(1));
            G1_ = spdiags(ones(T,1)*[-min(gr_),1],[-1:0],T,T);
            G2_ = spdiags(ones(T,1)*[-max(gr_),1],[-1:0],T,T);
            s_1_ = sparse(ceil(spiketimes_/Dt),1,exp((spiketimes_ - Dt*ceil(spiketimes_/Dt))/tau_(1)),T,1);  
            s_2_ = sparse(ceil(spiketimes_/Dt),1,exp((spiketimes_ - Dt*ceil(spiketimes_/Dt))/tau_(2)),T,1);  
            Gs_ = (-G1_\s_1_(:)+G2_\s_2_(:))/h_;
            
            logC_ = -(A_*Gs_'-Ym')*(A_*Gs_'-Ym')';
                
            %accept or reject
            prior_ratio = 1;
    %         prior_ratio = gampdf(tau_(2),12,1)/gampdf(tau(2),12,1);
            ratio = exp(sum(sum((1./(2*calciumNoiseVar)).*(logC_-logC))))*prior_ratio;
            if ratio>1 %accept
                tau = tau_;
                h = h_; G1 = G1_; G2 = G2_; s_1 = s_1_; s_2 = s_2_; Gs = Gs_;
                tauMoves = tauMoves + [1 1];
            elseif rand<ratio %accept
                tau = tau_;
                h = h_; G1 = G1_; G2 = G2_; s_1 = s_1_; s_2 = s_2_; Gs = Gs_;
                tauMoves = tauMoves + [1 1];
            else
                tauMoves = tauMoves + [0 1];
            end            
            
            %convert to g for storage or whatnot 
            g_new = gr_; %is this correct?
            
            %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
%             min_g = @(gam) min_gamma(gam,ss_(:),y_res(:));
%             g_new = fmincon(min_g,g,[],[],[],[],0,1,[],options);
            fprintf('new value %1.5f \n',g_new);
            g = g_new;
            P.g = g;
%             G = spdiags([-g*ones(T,1),ones(T,1)],[-1,0],T,T);
%             ge = P.g.^((0:T-1)'); 
%             ge(ge<prec) = 0;
%             tau(2) = -Dt/log(g);
        end
    end
    if mod(i,100)==0
        fprintf('%i out of total %i samples drawn \n', i, N);
    end 
end
if marg_flag
    mub = mub/(N-B);
    Sigb = Sigb/(N-B)^2;
end

if marg_flag
    SAMPLES.Cb = [mub(1),sqrt(Sigb(1,1))];
    SAMPLES.Cin = [mub(1+(1:p)),sqrt(diag(Sigb(1+(1:p),1+(1:p))))];
else
    SAMPLES.Cb = Cb(B+1:N);
    SAMPLES.Cin = Cin(B+1:N,:);
    SAMPLES.sn2 = SG(B+1:N).^2;
end
SAMPLES.ns = ns(B+1:N);
SAMPLES.ss = ss(B+1:N);
SAMPLES.ld = lam(B+1:N);
SAMPLES.Am = Am(B+1:N)/h;
if gam_flag
    SAMPLES.g = Gam(B+1:N,:);
end
SAMPLES.params = params.init;
