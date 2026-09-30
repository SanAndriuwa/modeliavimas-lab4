% Numerical differentiation and integration: all PDF variants.
clear; clc; close all;
format long g;
functions = {@(x) x.^3-2*x, @(x) sin(x), @(x) exp(x)};
points = [1 pi/3 0];
trueFirst = [1 0.5 1]; trueSecond = [6 -sqrt(3)/2 1];
steps = [0.1 0.01];
c13 = difapx(1,[-1 1]); c15 = difapx(1,[-2 2]);
c23 = difapx(2,[-1 1]); c25 = difapx(2,[-2 2]);
for variant = 1:3
    fprintf('\nDIFFERENTIATION VARIANT %d\n',variant);
    f = functions{variant}; x0 = points(variant);
    errors = zeros(2,4);
    for i = 1:2
        h = steps(i);
        y3 = f(x0+(-1:1)*h); y5 = f(x0+(-2:2)*h);
        estimates = [c13*y3'/h c15*y5'/h ...
            c23*y3'/h^2 c25*y5'/h^2];
        errors(i,:) = abs(estimates-[trueFirst(variant)*[1 1] ...
            trueSecond(variant)*[1 1]]);
    end
    disp(array2table([steps' errors],'VariableNames', ...
        {'h','First3','First5','Second3','Second5'}));
end

% Use the rounded data from the PDF exactly, including exp(-x) in set 3.
xData = [0.8 0.9 1 1.1 1.2; ...
    0.8472 0.9472 1.0472 1.1472 1.2472; -0.2 -0.1 0 0.1 0.2];
yData = [-1.0880 -1.0710 -1 -0.8690 -0.6720; ...
    0.7494 0.8118 0.8660 0.9116 0.9481; ...
    1.2214 1.1052 1 0.9048 0.8187];
for variant = 1:3
    fprintf('\nTABULATED DATA VARIANT %d\n',variant);
    x = xData(variant,:); y = yData(variant,:);
    h = x(2)-x(1); x0 = x(3);
    finite = [c13*y(2:4)'/h c15*y'/h c23*y(2:4)'/h^2 c25*y'/h^2];
    p2 = lagran(x(2:4),y(2:4)); p4 = lagran(x,y);
    derivatives = [polyval(polyder(p2),x0) polyval(polyder(p4),x0) ...
        polyval(polyder(polyder(p2)),x0) polyval(polyder(polyder(p4)),x0)];
    disp(array2table([finite; derivatives], ...
        'VariableNames',{'First3','First5','Second3','Second5'}, ...
        'RowNames',{'Finite differences','Lagrange derivatives'}));
end

integrands = {@(x) x.^3-2*x, @(x) sin(x), @(x) exp(-x)};
upper = [2 pi/2 1]; exact = [0 1 1-exp(-1)];
Ns = [4 8];
for variant = 1:3
    fprintf('\nFINITE INTEGRAL VARIANT %d\n',variant);
    f = integrands{variant}; errors = zeros(2,3);
    for i = 1:2
        N = Ns(i);
        estimates = [trapezoid(f,0,upper(variant),N) ...
            simpson(f,0,upper(variant),N) romberg(f,0,upper(variant),N)];
        errors(i,:) = abs(estimates-exact(variant));
    end
    disp(array2table([Ns' errors],'VariableNames', ...
        {'N','TrapezoidError','SimpsonError','RombergError'}));
end

% Separate numerical error on [a,b] from the omitted interval and tail.
a = 0.001; b = 100; N = 200; tolerance = 1e-4; gaussN = 20;
hasSymbolic = exist('syms','file') == 2;
quadratureErrors = zeros(2,5);
if hasSymbolic
    syms z positive
end
for variant = 1:2
    fprintf('\nIMPROPER INTEGRAL VARIANT %d\n',variant);
    if variant == 1
        f = @(x) sin(x)./x;
        infiniteExact = pi/2; symbolicText = 'pi/2';
        finiteExact = integral(f,a,b,'AbsTol',1e-12,'RelTol',1e-12);
        if hasSymbolic
            fSymbolic = sin(z)/z;
            finiteExact = double(sinint(sym(b))-sinint(sym(a)));
        end
    else
        f = @(x) exp(-x.^2);
        infiniteExact = sqrt(pi)/2; symbolicText = 'sqrt(pi)/2';
        finiteExact = sqrt(pi)/2*(erf(b)-erf(a));
        if hasSymbolic
            fSymbolic = exp(-z^2);
        end
    end
    if hasSymbolic
        infiniteSymbolic = int(fSymbolic,z,0,inf);
        infiniteExact = double(infiniteSymbolic);
        disp('Symbolic integral on [0,infinity):'); disp(infiniteSymbolic);
    else
        fprintf('Analytical integral on [0,infinity): %s\n',symbolicText);
        fprintf('Symbolic Math Toolbox is unavailable: int() was not executed.\n');
    end
    fprintf('Reference on [a,100]: %.12f\n',finiteExact);
    fprintf('Omitted interval and tail: %.6g\n',infiniteExact-finiteExact);
    estimates = zeros(5,1);
    estimates(1) = simpson(f,a,b,N);
    estimates(2) = adaptive_simpson(f,a,b,tolerance);
    estimates(3) = quad(f,a,b,tolerance);
    estimates(4) = quadl(f,a,b,tolerance);
    if variant == 1
        estimates(5) = gauss_legendre(f,a,b,gaussN);
        gaussName = 'Gauss-Legendre';
    else
        % Hermite integrates exp(-x^2)*g(x) on the whole real line.
        % For g=1 the half-line integral is half the sum of weights.
        estimates(5) = gauss_hermite_half(gaussN)-integral(f,0,a)-integral(f,b,inf);
        gaussName = 'Gauss-Hermite';
    end
    Method = {'Simpson'; 'Adaptive Simpson'; 'quad'; 'quadl'; gaussName};
    quadratureErrors(variant,:) = abs(estimates-finiteExact)';
    disp(table(Method,estimates,abs(estimates-finiteExact), ...
        abs(estimates-infiniteExact),'VariableNames', ...
        {'Method','Integral','QuadratureError','TotalErrorToInfinite'}));
    figure;
    plot(linspace(a,10,1000),f(linspace(a,10,1000)));
    xlabel('x'); ylabel('f(x)'); grid on;
    title(sprintf('Improper integral %d: integrand on [a,10]',variant));
end

function I = trapezoid(f,a,b,N)
h = (b-a)/N; y = f(linspace(a,b,N+1));
I = h*((y(1)+y(end))/2+sum(y(2:end-1)));
end

function I = simpson(f,a,b,N)
h = (b-a)/N; y = f(linspace(a,b,N+1));
I = h/3*(y(1)+y(end)+4*sum(y(2:2:end-1))+2*sum(y(3:2:end-2)));
end

function I = romberg(f,a,b,N)
% The final row uses exactly N intervals; no early stopping.
levels = round(log2(N))+1; R = zeros(levels);
for row = 1:levels
    R(row,1) = trapezoid(f,a,b,2^(row-1));
    for col = 2:row
        R(row,col) = (4^(col-1)*R(row,col-1)-R(row-1,col-1))/(4^(col-1)-1);
    end
end
I = R(end,end);
end

function I = adaptive_simpson(f,a,b,tolerance)
% Begin with short intervals so oscillations cannot be missed on [a,100].
segments = max(1,ceil(b-a));
edges = linspace(a,b,segments+1);
I = 0;
for i = 1:segments
    left = edges(i); right = edges(i+1);
    whole = simpson(f,left,right,2);
    I = I+refine(f,left,right,whole,tolerance/segments,0);
end
end

function I = refine(f,a,b,whole,tolerance,depth)
mid = (a+b)/2;
left = simpson(f,a,mid,2); right = simpson(f,mid,b,2);
difference = left+right-whole;
if abs(difference) <= 15*tolerance
    I = left+right+difference/15;
elseif depth >= 30
    error('Adaptive Simpson reached the recursion limit.');
else
    I = refine(f,a,mid,left,tolerance/2,depth+1) ...
        +refine(f,mid,b,right,tolerance/2,depth+1);
end
end

function I = gauss_legendre(f,a,b,N)
% Nodes and weights from the symmetric Legendre recurrence matrix.
j = 1:N-1; offDiagonal = j./sqrt(4*j.^2-1);
J = diag(offDiagonal,1)+diag(offDiagonal,-1);
[vectors,D] = eig(J);
nodes = diag(D)'; weights = 2*vectors(1,:).^2;
x = (a+b)/2+(b-a)/2*nodes;
I = (b-a)/2*sum(weights.*f(x));
end

function I = gauss_hermite_half(N)
offDiagonal = sqrt((1:N-1)/2);
[vectors,~] = eig(diag(offDiagonal,1)+diag(offDiagonal,-1));
weights = sqrt(pi)*vectors(1,:).^2;
I = sum(weights)/2;
end
