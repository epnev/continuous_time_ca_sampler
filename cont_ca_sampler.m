function SAMPLES = cont_ca_sampler(Y,P,Nsamples,B,params)

% Continuous time sampler
% Y                     data (normalized in [0,1])
% P                     intialization parameters (discrete time constant P.g required)
% Nsamples              number of samples (after burn-in period)
% B                     length of burn-in period
% params                optional additional parameters
% params.marg           flag for marginalized sampler (default 1)
% params.upd_gam        flag for updating gamma (default 0)
% params.gam_step       number of samples after which gamma is updated (default 50)

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

% Author: Eftychios A. Pnevmatikakis

if nargin == 4
    marg_flag = 1;
    gam_flag = 0;
else
    if isfield(params,'marg')
        marg_flag = params.marg;
    else
        marg_flag = 1;
    end
    if isfield(params,'upd_gam')
        gam_flag = params.upd_gam;
    else
        gam_flag = 0;
    end
    if ~isfield(params,'gam_step')
        gam_step = 50;
    else
        gam_step = params.gam_step;
    end
end

if gam_flag
    options = optimset('GradObj','On','Display','Off','Algorithm','interior-point','TolX',1e-6);
end

T = length(Y);
g = P.g;
Dt = 1/P.f;                                     % length of time bin
tau = -Dt/log(g); %Dt/(1-g);                    % continuous time constant
if ~isfield(P,'Cb')
    P.Cb = mean(Y(1:round(T/2)))/2;             % set an arbitrary baseline initializer
end
if ~isfield(P,'sn')
    temp = (xcov(Y,1))/T;
    P.sn = sqrt((g*temp(2) - temp(1))/g);
    if P.sn^2<0
        temp = diff(xcov(Y,1))/T;
        P.sn = sqrt(temp(1));
        if P.sn^2 < 0
            P.sn = std(Y(1:round(T/2)));
        end
    end
end
sg = P.sn;

fprintf('Initializing using noise constrained FOOPSI...  ');
[Z,~,~] = optimal_foopsi_lars(Y(:)',g,P.sn,P.Cb);
fprintf('done. \n');

G = spdiags([-g*ones(T,1),ones(T,1)],[-1,0],T,T);
sp = G*Z(:);                                   % extract spikes
c1 = sp(1);
sp(1) = 0;
s_in = sp>0.45*max(sp);
spiketimes_ = Dt*(find(s_in) + rand(size(find(s_in))) - 0.5);
spiketimes_(spiketimes_ >= T*Dt) = 2*T*Dt - spiketimes_(spiketimes_ >= T*Dt);
lam_ = length(spiketimes_)/(T*Dt);
s_ = sparse(ceil(spiketimes_/Dt),1,exp(-(spiketimes_ - Dt*ceil(spiketimes_/Dt))/tau),T,1);

A_   = max(median(sp(s_in)),0.11);             % initial amplitude value
b_   = P.Cb;                                   % initial baseline value'
C_in = c1;       % initial value sample
ge = P.g.^((0:T-1)'); 
prec = 1e-2;     % precision
ge(ge<prec) = 0;

N = Nsamples + B;

ss = cell(N,1); 
lam = zeros(N,1);
Am = zeros(N,1);
ns = zeros(N,1);
Gam = zeros(N,1);
if ~marg_flag
    Cb = zeros(N,1);
    Cin = zeros(N,1);
    SG = zeros(N,1);
end

Sp = .1*eye(3);          % prior covariance
Ld = inv(Sp);
mu = [A_;b_;C_in];       % prior mean 
lb = [0.1,0.02,0]';      % lower bound for [A,Cb,Cin]
Ns = 15;                 % Number of HMC samples
mu_b = mu(2:3);
Ym = Y - [ones(T,1),ge]*mu_b;

mub = zeros(2,1);
Sigb = zeros(2,2);

for i = 1:N
    if gam_flag
        Gam(i) = g;
    end
    sg_ = sg;
    rate = @(t) lambda_rate(t,lam_);
    [spiketimes, ~]  = get_next_spikes(spiketimes_(:)',(A_*(G\s_(:)))',Ym',ge',tau,sg_^2, rate, 20*Dt, Dt, A_);

    spiketimes_ = spiketimes;
    spiketimes(spiketimes<0) = -spiketimes(spiketimes<0);
    spiketimes(spiketimes>T*Dt) = 2*T*Dt - spiketimes(spiketimes>T*Dt); 
    trunc_spikes = ceil(spiketimes/Dt);
    trunc_spikes(trunc_spikes == 0) = 1;
    s_ = sparse(1,trunc_spikes,exp((spiketimes - Dt*trunc_spikes)/tau),1,T);
    ss{i} = spiketimes;
    nsp = length(spiketimes);
    ns(i) = nsp;
    lam(i) = nsp/(T*Dt);
    lam_ = lam(i);
    
    AM = [G\s_(:),ones(T,1),ge];
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
        Cin(i) = temp(3,Ns);
        A_ = Am(i);
        b_ = Cb(i);
        C_in = Cin(i);

        Ym   = Y - b_ - ge*C_in;
        res   = Ym - A_*filter(1,[1,-g],full(s_'));
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
           mub = mub + mu_post(2:3);
           Sigb = Sigb + L(2:3,2:3);
        end
    end
    if gam_flag
        if (i > B) && mod(i-B,gam_step) == 0  % update gamma
            fprintf('updating time constant.. ')
            ss_ = zeros(1,T);
            for rep = i-gam_step+1:i
                trunc_spikes = ceil(ss{rep}/Dt);            
                ss_ = ss_ + Am(rep)*sparse(1,trunc_spikes,exp((ss{rep} - Dt*trunc_spikes)/tau),1,T);
            end
            ss_ = ss_/gam_step;        
            if ~marg_flag
                y_res = Y - mean(Cb(i-gam_step+1:i));
                ss_(1) = ss_(1) + mean(Cin(i-gam_step+1:i));
            else
                y_res = Y - mub(1)/(i-B);
                ss_(1) = ss_(1) + mub(2)/(i-B);
            end
        
            min_g = @(gam) min_gamma(gam,ss_(:),y_res(:));
            g_new = fmincon(min_g,g,[],[],[],[],0,1,[],options);
            fprintf('new value %1.5f \n',g_new);
            g = g_new;
            P.g = g;
            G = spdiags([-g*ones(T,1),ones(T,1)],[-1,0],T,T);
            ge = P.g.^((0:T-1)'); 
            ge(ge<prec) = 0;
            tau = -Dt/log(g);
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
    SAMPLES.Cin = [mub(2),sqrt(Sigb(2,2))];
else
    SAMPLES.Cb = Cb(B+1:N);
    SAMPLES.Cin = Cin(B+1:N);
    SAMPLES.sn2 = SG(B+1:N).^2;
end
SAMPLES.ns = ns(B+1:N);
SAMPLES.ss = ss(B+1:N);
SAMPLES.ld = lam(B+1:N);
SAMPLES.Am = Am(B+1:N);
if gam_flag
    SAMPLES.g = Gam(B+1:N);
end