function plot_continuous_samples(SAMPLES,P,Y)

T = length(Y);
N = length(SAMPLES.ns);
show_gamma = 1;
P.f = 1;
Dt = 1/P.f;                                     % length of time bin
if ~isfield(SAMPLES,'g');
    show_gamma = 0;
    SAMPLES.g = P.g*ones(N,1);
end
g = SAMPLES.g;
tau = -Dt./log(g);

if length(SAMPLES.Cb) == 2
    marg = 1;       % marginalized sampler
else
    marg = 0;       % full sampler
end

C_rec = zeros(N,T);
for rep = 1:N
    trunc_spikes = ceil(SAMPLES.ss{rep}/Dt);
    s_ = sparse(1,trunc_spikes,exp((SAMPLES.ss{rep} - Dt*trunc_spikes)/tau(rep)),1,T);
    if marg
        C_rec(rep,:) = SAMPLES.Cb(1) + SAMPLES.Am(rep)*filter(1,[1,-g(rep)],full(s_)+[SAMPLES.Cin(1),zeros(1,T-1)]);
    else
        C_rec(rep,:) = SAMPLES.Cb(rep) + SAMPLES.Am(rep)*filter(1,[1,-g(rep)],full(s_)+[SAMPLES.Cin(rep),zeros(1,T-1)]);
    end
end
Nc = 60;

if marg
    rows = 4;
else
    rows = 5;
end

figure;
    set(gcf, 'PaperUnits', 'inches','Units', 'inches')           
    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'PaperPosition',[0,0, 14, 15])
    set(gcf, 'Position',[2,2, 14, 15])
    subplot(rows,4,[1:4]);plot(Dt*(1:T),Y); hold all; plot(Dt*(1:T),mean(C_rec,1),'linewidth',2); 
        title('Calcium traces','fontweight','bold','fontsize',14)
        legend('Raw data','Mean sample');
    subplot(rows,4,[5:8]); imagesc((1:T)*Dt,1:N,samples_cell2mat(SAMPLES.ss,T)); 
        title('Spike raster plot','fontweight','bold','fontsize',14)
    subplot(rows,4,4+5); plot(1:N,SAMPLES.ns); title('# of spikes','fontweight','bold','fontsize',14)
    subplot(rows,4,4+6); plot(-Nc:Nc,xcov(SAMPLES.ns,Nc,'coeff')); set(gca,'XLim',[-Nc,Nc]);
        title('Autocorrelation','fontweight','bold','fontsize',14)
    subplot(rows,4,4+7); plot(1:N,SAMPLES.ld); title('Firing Rate','fontweight','bold','fontsize',14)
    if ~show_gamma
        subplot(rows,4,4+8); plot(-Nc:Nc,xcov(SAMPLES.ld,Nc,'coeff')); set(gca,'XLim',[-Nc,Nc])
        title('Autocorrelation','fontweight','bold','fontsize',14)
    else
        subplot(rows,4,4+8);  plot(1:N,g); title('Discrete time constant','fontweight','bold','fontsize',14)
    end
    
    subplot(rows,4,4+9); plot(1:N,SAMPLES.Am); title('Spike Amplitude','fontweight','bold','fontsize',14)
    subplot(rows,4,4+10); plot(-Nc:Nc,xcov(SAMPLES.Am,Nc,'coeff')); set(gca,'XLim',[-Nc,Nc])
        title('Autocorrelation','fontweight','bold','fontsize',14)
    if marg
        xx = SAMPLES.Cb(1) + linspace(-4*SAMPLES.Cb(2),4*SAMPLES.Cb(2));
        subplot(4,4,15); plot(xx,normpdf(xx,SAMPLES.Cb(1),SAMPLES.Cb(2)));
            set(gca,'XLim',[xx(1),xx(end)])
            title('Marg. post. of baseline','fontweight','bold','fontsize',14)
        xx = SAMPLES.Cin(1) + linspace(-4*SAMPLES.Cin(2),4*SAMPLES.Cin(2));
        subplot(4,4,16); plot(xx,normpdf(xx,SAMPLES.Cin(1),SAMPLES.Cin(2)));
            set(gca,'XLim',[xx(1),xx(end)])
            title('Marg. post. of initial con','fontweight','bold','fontsize',14)
    else
        subplot(5,4,4+11); plot(1:N,SAMPLES.Cb); title('Baseline','fontweight','bold','fontsize',14)
        subplot(5,4,4+12); plot(-Nc:Nc,xcov(SAMPLES.Cb,Nc,'coeff')); set(gca,'XLim',[-Nc,Nc])
            title('Autocorrelation','fontweight','bold','fontsize',14)
        subplot(5,4,4+13); plot(1:N,SAMPLES.Cin); title('Initial Concentration','fontweight','bold','fontsize',14)
        subplot(5,4,4+14); plot(-Nc:Nc,xcov(SAMPLES.Cin,Nc,'coeff')); set(gca,'XLim',[-Nc,Nc])
            title('Autocorrelation','fontweight','bold','fontsize',14)
        subplot(5,4,4+15); plot(1:N,SAMPLES.sn2); title('Noise variance','fontweight','bold','fontsize',14)
        subplot(5,4,4+16); plot(-Nc:Nc,xcov(SAMPLES.sn2,Nc,'coeff')); set(gca,'XLim',[-Nc,Nc])
            title('Autocorrelation','fontweight','bold','fontsize',14)
    end
    drawnow;