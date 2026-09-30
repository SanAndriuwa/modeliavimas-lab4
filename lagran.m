function [l,L]=lagran(x,y)

%Input : x = [x0 x1 ... xN], y = [y0 y1 ... yN]
%Output: l = N eiles Lagranzo daugianaris
% L = Lagranzo koeficientu daugianaris
N = length(x)-1; %daugianario eile
l = 0;
L = zeros(N+1,N+1);
for m = 1:N + 1
    P = 1;
    for k = 1:N + 1
        if k ~= m, P = conv(P,[1 -x(k)])/(x(m)-x(k)); end
        %MATLAB dirba su daugianariu koeficientu vektoriumi (mazejancia
        %tvarka), o dvieju daugianariu daugyba atitinka koeficientu
        %vektoriu konvoliucijai...
    end
    L(m,:) = P; %Lagranzo koeficientu daugianaris
    l = l + y(m)*P; %Lagranzo daugianaris (3 lytis)
end
