function c_m = make_mean_sample(SAMPLES,P,T,ind)

N = length(SAMPLES.ns);
if nargin == 3
    ind = 1:N;
end
ind(ind>N) = [];
P.f = 1;
g = P.g(:);
p = length(g);
Dt = 1/P.f;                                     % length of time bin
if ~isfield(SAMPLES,'g');
    SAMPLES.g = ones(N,1)*g';
end

if p == 1
    tau_1 = 0;
    tau_2 = -Dt/log(g);                         % continuous time constant
    G1 = speye(T);
    G2 = spdiags(ones(T,1)*[-g,1],[-1:0],T,T);
    ge = P.g.^((0:T-1)');     
elseif p == 2
    gr = roots([1,-g']);
    p1_continuous = log(min(gr))/Dt; 
    p2_continuous = log(max(gr))/Dt;
    tau_1 = -1/p1_continuous;                   %tau h - smaller (tau_d * tau_r)/(tau_d + tau_r)
    tau_2 = -1/p2_continuous;                   %tau decay - larger
    G1 = spdiags(ones(T,1)*[-min(gr),1],[-1:0],T,T);
    G2 = spdiags(ones(T,1)*[-max(gr),1],[-1:0],T,T);
    ge = [-G1\[1;zeros(T-1,1)],G2\[1;zeros(T-1,1)]];
else
    error('This order of the AR process is currently not supported');
end



if length(SAMPLES.Cb) == 2
    marg = 1;       % marginalized sampler
else
    marg = 0;       % full sampler
end

C_rec = zeros(length(ind),T);
for i = 1:length(ind)
    rep = ind(i);
    %trunc_spikes = ceil(SAMPLES.ss{rep}/Dt);
    s_1 =   sparse(ceil(SAMPLES.ss{rep}/Dt),1,exp((SAMPLES.ss{rep} - Dt*ceil(SAMPLES.ss{rep}/Dt))/tau_1),T,1);  
    s_2 =   sparse(ceil(SAMPLES.ss{rep}/Dt),1,exp((SAMPLES.ss{rep} - Dt*ceil(SAMPLES.ss{rep}/Dt))/tau_2),T,1);  
    Gs = (-G1\s_1(:)+ G2\s_2(:));
    if marg
        %C_rec(rep,:) = SAMPLES.Cb(1) + SAMPLES.Am(rep)*filter(1,[1,-SAMPLES.g(rep,:)],full(s_)+[SAMPLES.Cin(:,1)',zeros(1,T-p)]);
        C_rec(i,:) = SAMPLES.Cb(1) + SAMPLES.Am(rep)*Gs + (ge*SAMPLES.Cin(:,1));
    else
        %C_rec(rep,:) = SAMPLES.Cb(rep) + SAMPLES.Am(rep)*filter(1,[1,-SAMPLES.g(rep,:)],full(s_)+[SAMPLES.Cin(rep,:),zeros(1,T-p)]);
        C_rec(i,:) = SAMPLES.Cb(rep) + SAMPLES.Am(rep)*Gs + (ge*SAMPLES.Cin(rep,:)');
    end
end

if length(ind)>1
    c_m = mean(C_rec);
else
    c_m = C_rec;
end