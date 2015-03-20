function P = arpfit(Y,p,ph)

% estimate a global time constant of order p and noise power for each pixel

if ndims(Y) == 3
    Y = reshape(Y,size(Y,1)*size(Y,2),size(Y,3));
end

[d,T] = size(Y);

if nargin == 2
    fprintf('Estimating time constant through autocorrelation function.. \n');
    tt1 = tic;
    lags = 15;
    XC = zeros(d,2*(lags+p)+1);

    for j = 1:d
        XC(j,:) = xcov(Y(j,:),lags+p,'unbiased');
    end

    g = XC(:,lags+p:-1:1);
    
    A = zeros(d*lags,p);
    for i = 1:d
        A((i-1)*lags + (1:lags),:) = toeplitz(g(i,p:p+lags-1),g(i,p:-1:1));
    end
    gv = g(:,p+1:end)';
    %ph = pinv(A)*gv(:);
    ph = A\gv(:);
    disp(ph);
    fprintf('Done after %2.2f seconds. \n',toc(tt1));
end
mg = sum(ph);

sn = zeros(d,1);

cov_flag = exist('XC','var');

for i = 1:d
    if cov_flag
        xx = XC(i,lags+(1:2*p+1))';
    else
        xx = xcov(Y(i,:)',p,'unbiased');
    end
    s_est = zeros(p,1);
    for k = 1:p
        s_est(k) = (ph'*xx(p+1+k-(1:p)) - xx(p+1+k))/ph(k);
    end
    
    sn(i) = sqrt(mean(s_est));
end

G = spdiags(ones(T,1)*[-ph(end:-1:1)',1],[-length(ph):0],T,T); %foopsi_matrix(T,P.g);
Sp = Y*G';
Sp = Sp(:,p+1:end);
Cb = max(quantile(Sp,0.3,2),0)/(1-mg);

P.sn = sn;
P.Cb = Cb;
P.g = ph(:);

ind = find(abs(imag(P.sn)>0));
if ~isempty(ind);
    fprintf('Correcting complex number estimates. \n');
    %ind2 = setdiff(1:d,ind);
    P.sn(ind) = std(Y(ind,:),[],2); %mean(P.sn(ind2));
    P.Cb(ind) = quantile(Y(ind,:),0.1,2); %mean(P.Cb(ind2));
    fprintf('done. \n');
end