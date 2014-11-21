function [Z,Cb,SNR] = optimal_foopsi_lars(y,g,sn,Cb)

% Y     vector of data in column format with 0 <= min(Y) < max(Y) <= 1
% g     discrete time constant
% sn    noise standard deviation
% b     baseline (optional)

% outputs
% Z     calcium trace
% b     baseline

T = length(y);
G = spdiags(ones(T,1)*[-g(end:-1:1),1],-length(g):0,T,T);
if T > 5e3;
   use_cvx = 1;
end
if nargin == 3
    if use_cvx
        cvx_begin quiet
            variable Z(T)
            variable Cb
            minimize(sum(G*Z))
            subject to
                G*Z>=0;
                norm(y'-Z-Cb)<=sqrt(T)*sn;
                Cb>=0;
        cvx_end
        Z = Z';
    else
        Ginv = [full(G\speye(T)),ones(T,1)];
        [~, ~, spikes, ~, ~] = lars_regression_noise(y', Ginv, 1, sn^2*T);
        Cb = spikes(end);
        Z = filter(1,[1,-g],spikes(1:T));
    end
elseif nargin == 4
    if use_cvx
        cvx_begin quiet
            variable Z(T)
            minimize(sum(G*Z))
            subject to
                G*Z>=0;
                norm(y'-Z-Cb)<=sqrt(T)*sn;
        cvx_end
        Z = Z';
    else
        Ginv = full(G\speye(T));
        %Ginv = toeplitz(g.^(0:T-1),[1,zeros(1,T-1)]);
        [~, ~, spikes, ~, ~] = lars_regression_noise(y'-Cb, Ginv, 1, sn^2*T);
        Z = filter(1,[1,-g],spikes(1:T));
    end
else
    error('wrong number of inputs provided');
end


SNR = 10*log10(norm(Z)^2/norm(y(:)-Z(:)-Cb)^2);