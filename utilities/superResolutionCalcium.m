%dimensions
T = 1000; %start with large T and sample in this space

%model parameters
gamma = .85;
%I think I got this next part right
gamma_continuous = log(gamma)/1; %divide by bin size (I am treating bin sizes as unit length)
tau = -1/gamma_continuous;
b = 1;
A = 1;
p_spike = .01;
c_noise = .2;
% c_noise = 1;
% c_noise = 0; %noiseless calcium signal

%start with univariate calcium signal
s = rand(1,T)<p_spike;

c = b*ones(1,T);
y = zeros(1,T);
for t = 1:T
    if t>1
        c(t) = (1-gamma)*b + gamma*c(t-1) + A*s(t);
    else
        c(t) = (1-gamma)*b + A*s(t);
    end
    y(t) = c(t) + c_noise*randn;
end

%%
figure(45)
subplot(211)
plot(s)
ylim([-.5 2])
subplot(212)
plot(c);hold on; plot(y,'r--');hold off



%% compute exponential filter
T2 = 1000; %this might need to be bigger if sampling rate is high or decay is slow
%note - b is the steady state
c2 = b*ones(1,T2);
for t = 1:T
    if t==1
        c2(t) = (1-gamma)*b + gamma*c2(1) + A*1;
    else
        c2(t) = (1-gamma)*b + gamma*c2(t-1);
    end
end

e_support = find((c2-b)<1e-3,1);
ef = c2(1:e_support)-b;
figure;plot(ef)


%% try sampling function
proposal_move_var = 3; %tune to get a good acceptance rate
[samples, addMoves, dropMoves, timeMoves, N_sto]  = sampleSpikes(y,ef,tau,b,c_noise,p_spike,proposal_move_var,50);

%%
figure;
plot(s,'r')
hold on;
plot(mean(samples_cell2mat(samples,T)))
hold off


