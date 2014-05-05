function plot_continuous_samples_l(SAMPLES,P,Y,traceData)

T = length(Y);
N = length(SAMPLES.ns);

Dt = 1/P.f;                                     % length of time bin
if ~isfield(SAMPLES,'g');
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
    rows = 3;
else
    rows = 5;
end
Mat = samples_cell2mat(SAMPLES.ss,T);
M = zeros(5,T);
for i = 1:T
    M(:,i) = hist(Mat(:,i),0:4);
end
    
figure;
    set(gcf, 'PaperUnits', 'inches','Units', 'inches')           
    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'PaperPosition',[0,0, 14, 15])
    set(gcf, 'Position',[2,2, 14, 15])
    subplot(rows,4,[1:4]);plot(Dt*(1:T),Y); hold all; plot(Dt*(1:T),mean(C_rec,1),'linewidth',2); 
        title('Calcium traces  ','fontweight','bold','fontsize',14)
        legend('Raw data ','Mean sample'); set(gca,'XTick',[],'Ytick',[])
    subplot(rows,4,[5:8]); imagesc((1:T)*Dt,1:N,Mat);  set(gca,'Xtick',[])
        title('Spike raster plot','fontweight','bold','fontsize',14);
        ylabel('Sample # ','fontweight','bold','fontsize',14);
    subplot(rows,4,[9:12]); imagesc(M/1000); axis xy; set(gca,'Ytick',[0.5:4.5],'Yticklabel',[-1:4]); colormap('gray'); %hold all; plot(traceData.spikeFrames+1); set(gca,'YLim',[1,5])
set(gca,'YLim',[1.5,5]);
hold all; scatter(1:T,traceData.spikeFrames+1.5);
 ylabel('# of Spikes ','fontweight','bold','fontsize',14);
 xlabel('Timestep ','fontweight','bold','fontsize',14);
 title('Spike Histogram ','fontweight','bold','fontsize',14);
    drawnow;