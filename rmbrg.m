function [x,R,err,N] = rmbrg(f,a,b,tol,K)
% Romberg.m from Int.zip. tol=0 forces all K rows for the N=4,8 task.
if nargin < 5, K = 10; end
if nargin < 4, tol = 1e-6; end
h = b-a; N = 1; err = inf;
R = zeros(K);
R(1,1) = h/2*(feval(f,a)+feval(f,b));
for k = 2:K
    h = h/2; N = N*2;
    R(k,1) = R(k-1,1)/2+h*sum(feval(f,a+(1:2:N-1)*h));
    factor = 1;
    for n = 2:k
        factor = factor*4;
        R(k,n) = (factor*R(k,n-1)-R(k-1,n-1))/(factor-1);
    end
    err = abs(R(k,k-1)-R(k-1,k-1))/(factor-1);
    if err < tol, break; end
end
if K == 1, k = 1; end
x = R(k,k);
R = R(1:k,1:k);
end
