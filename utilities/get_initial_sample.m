function SAM = get_initial_sample(Y,P,Z)

if nargin == 2
    if isfield(P,'Cb');
        [Z,~,~] = optimal_foopsi_lars(Y(:)',P.g,P.sn,P.Cb);
    else
        [Z,P.Cb,~] = optimal_foopsi_lars(Y(:)',P.g,P.sn);
    end
end

Dt = 1;
T = length(Y);
p = length(P.g);
G = make_G_matrix(T,P.g);
sp = G*Z(:);                                   % extract spikes
c1 = sp(1:p);
sp(1:p) = 0;
s_in = sp>0.15*max(sp);
spiketimes_ = Dt*(find(s_in) + rand(size(find(s_in))) - 0.5);
spiketimes_(spiketimes_ >= T*Dt) = 2*T*Dt - spiketimes_(spiketimes_ >= T*Dt);
SAM.lam_ = length(spiketimes_)/(T*Dt);
SAM.spiketimes_ = spiketimes_;

SAM.A_   = max(median(sp(s_in)),0.11);             % initial amplitude value
SAM.b_   = P.Cb;                                   % initial baseline value'
SAM.C_in = c1 + 1e-3;                              % initial value sample
SAM.sg = P.sn;