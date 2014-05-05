function [f,grad] = min_gamma(g,s,y)

T = length(s);
prec = 1e-4;
K = min(ceil(log(prec)/log(g)),T);
vec_f = g.^(0:K-1);
G_f = toeplitz(sparse(1:K,1,vec_f(:),T,1),sparse(1,1,1,1,T));
vec_g = (0:K-1).*(g.^(-1:K-2));
G_g = toeplitz(sparse(1:K,1,vec_g(:),T,1),sparse(1,T));

f = 0.5*norm(y-G_f*s)^2;
grad = -(y-G_f*s)'*(G_g*s);