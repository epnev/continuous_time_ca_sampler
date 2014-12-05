function [f,grad] = min_gamma(g,s,y)

T = length(s);
prec = 1e-4;

if length(g) == 1
    K = min(ceil(log(prec)/log(g)),T);
    vec_f = g.^(0:K-1);
    G_f = toeplitz(sparse(1:K,1,vec_f(:),T,1),sparse(1,1,1,1,T));
    vec_g = (0:K-1).*(g.^(-1:K-2));
    G_g = toeplitz(sparse(1:K,1,vec_g(:),T,1),sparse(1,T));

    f = 0.5*norm(y-G_f*s)^2;
    grad = -(y-G_f*s)'*(G_g*s);
else
    G = spdiags(ones(T,1)*[-flipud(g(:))',1],-length(g):0,T,T);
    Gi = G\[1;zeros(T-1,1)];
    K = min(find(Gi<prec,1,'first'),T);
    Gs = G\s;
    f = 0.5*norm(y-Gs)^2;
    g1 = g(1); g2 = g(2);
    d = sqrt(g1^2+4*g2);
    v1 = 2.^(1:K).*(g1*((g1-d).^(1:K) - (g1+d).^(1:K)) + d*((g1-d).^(1:K)+(g1+d).^(1:K)).*(1:K))/d^3;
    v2 = 2.^(1-(1:K)).*((g1-d).^(1:K) - (g1+d).^(1:K) + d*((g1-d).^(0:K-1) + (g1+d).^(0:K-1)).*(1:K))/d^3;
    G_g1 = toeplitz(sparse(1:K,1,v1(:),T,1),sparse(1,T));
    G_g2 = toeplitz(sparse(1:K,1,v2(:),T,1),sparse(1,T));
    grad = -((y-Gs)'*[G_g1*s,G_g2*s])';
end
    
    