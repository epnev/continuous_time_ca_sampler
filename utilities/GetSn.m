function [ sn, thr ] = GetSn(Y)
%GETSN Summary of this function goes here
% input:
% Y - calcium trace
% output:
% Sn - noise standard diviation in trace
% thr - residual noise ( = norm( noise trace) )
L=length(Y);
[psd_Y,ff]=pwelch(Y,round(L/8),[],1000,1);
%[psd_Y,ff]=pwelch(Y,512,[],1000,1);

semilogy(ff,psd_Y);
%max(psd_Y)
xlim([min(ff) max(ff)])

range_min=0.1;
range_max=0.5; % slightly lower then where data was cut-off
ind=ff<range_max;
ind(ff<range_min)=0;
sn=sqrt(median(psd_Y(ind)/2));
thr=sn*sqrt(L);

end

