function [I,t,w] = Gauss_Hermite(g,N,varargin)
% Corrected supplied Gauss_Hermite.m.
% Integrates exp(-x^2)*g(x) on the whole real line; no finite-interval map.
[t,w] = Gausshp(N);
values = feval(g,t,varargin{:});
I = w*values(:);
end

function [t,w] = Gausshp(N)
if N < 1, error('N must be positive.'); end
t = roots(Hermitp(N))';
A = zeros(N); b = zeros(1,N);
A(1,:) = ones(1,N); b(1) = sqrt(pi);
for n = 2:N
    A(n,:) = A(n-1,:).*t;
    if mod(n,2) == 1, b(n) = (n-2)/2*b(n-2); end
end
w = b/A';
end

function p = Hermitp(N)
% The Hermite recurrence from the supplied function.
if N == 0
    p = 1;
else
    p = [2 0];
    for n = 2:N
        p = 2*[p 0]-[0 0 polyder(p)];
    end
end
end
