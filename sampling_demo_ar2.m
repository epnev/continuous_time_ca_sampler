clear all;
addpath utilities
dt = 5e-3;
T = 7000;
ld = 0.1;   % rate spikes per second

s = rand(1,round(T/dt)) < ld*dt;
%s = zeros(1,round(T/dt));
%s(0.95*round(T/dt)) = 1;
tau_rise = 0.3;
tau_decay = 1.5;
hmax = tau_decay/(tau_decay+tau_rise)*(tau_rise/(tau_decay+tau_rise))^(tau_rise/tau_decay);
[g,h1] = tau_c2d(tau_rise,tau_decay,dt);

b = hmax/4;
cin = [.2*b,.15*b];
c = [cin,filter(h1,[1,-g],s(3:end),filtic(h1,[1,-g],cin))] + b;
c_true = c(round(1/dt):round(1/dt):round(T/dt));
sg = hmax/4;        % noise level
y = c_true + sg*randn(1,length(c_true));

figure;plot(dt:dt:T,c); hold all;stem(1:T,y); drawnow;
        legend('True Calcium','Observed Values');
%%  constrained foopsi
[g2,h2] = tau_c2d(tau_rise,tau_decay,1);
%P = arpfit(y,2);  % uncomment to estimate parameters
P.f = 1;
P.g = g2;
P.sn = sg;
params.marg = 0;
y_r = y(:);
[ca_foopsi,cb] = optimal_foopsi_lars(y(:)',g2,sg);  %% constrained foopsi
G = make_G_matrix(T,g2);
spikes_foopsi = G*ca_foopsi(:);
spiketimes{1} =  find(s)*dt;
spikeRaster = samples_cell2mat(spiketimes,T,1);
f = find(spikeRaster);
spikes = zeros(sum(spikeRaster),2);
count = 0;
for cnt = 1:length(f)
    spikes(count+(1:spikeRaster(f(cnt))),1) = f(cnt);
    spikes(count+(1:spikeRaster(f(cnt))),2) = 0.95 + 0.025*max(spikes_foopsi)*(1:spikeRaster(f(cnt)));
    count = count + spikeRaster(f(cnt));
end

figure;

stem(spikes_foopsi); hold all; 
    scatter(spikes(:,1),spikes(:,2)-0.95+max(spikes_foopsi),15,'magenta','filled');
    axis([1,T,0,max(spikes(:,2))-0.95+max(spikes_foopsi)]);
    title('Foopsi Spikes','FontWeight','bold','Fontsize',14); xlabel('Timestep','FontWeight','bold','Fontsize',16);
    legend('Foopsi Spikes','Ground Truth');
    
%% MCMC    
SAMPLES = cont_ca_sampler(y_r,P,500,500,params);    %% MCMC        
plot_continuous_samples(SAMPLES,P,y_r(:));