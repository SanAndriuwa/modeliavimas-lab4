function INTf = smpsns(f,a,b,N,varargin)
% Simpsons.m from Int.zip; same Simpson weights and even-N rule.
if nargin < 4, N = 100; end
if abs(b-a) < 1e-12 || N <= 0, INTf = 0; return; end
if mod(N,2) ~= 0, N = N+1; end
h = (b-a)/N; x = a+(0:N)*h;
fx = feval(f,x,varargin{:});
fx(fx == inf) = realmax; fx(fx == -inf) = -realmax;
kodd = 2:2:N; keven = 3:2:N-1;
INTf = h/3*(fx(1)+fx(N+1)+4*sum(fx(kodd))+2*sum(fx(keven)));
end
