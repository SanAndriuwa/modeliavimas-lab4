function I = Gauss_Legendre(f,a,b,N,varargin)
% Supplied Gauss_Legendre.m: polynomial roots, moment weights, interval map.
[t,w] = Gausslp(N);
x = ((b-a)*t+a+b)/2;
fx = feval(f,x,varargin{:});
I = w*fx(:)*(b-a)/2;
end

function [t,w] = Gausslp(N)
if N < 1, error('N must be positive.'); end
t = roots(Lgndrp(N))';
A = zeros(N); b = zeros(1,N);
A(1,:) = ones(1,N); b(1) = 2;
for n = 2:N
    A(n,:) = A(n-1,:).*t;
    if mod(n,2) ~= 0, b(n) = 2/n; end
end
w = b/A';
end

function p = Lgndrp(N)
% The Legendre recurrence from the supplied function.
if N == 0
    p = 1;
elseif N == 1
    p = [1 0];
else
    l1 = (2*N-1)*[Lgndrp(N-1) 0];
    l2 = (N-1)*[0 0 Lgndrp(N-2)];
    p = (l1-l2)/N;
end
end
