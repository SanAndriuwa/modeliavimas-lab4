% Skaitinis diferencijavimas ir integravimas.
clear; clc; close all;
format long g;
variantas = 1; % 1, 2 arba 3: nurodo destytojas.
switch variantas
    case 1
        f = @(x) x.^3-2*x; x0 = 1;
        tikros = [1 6];
        x = [0.8 0.9 1 1.1 1.2];
        y = [-1.0880 -1.0710 -1 -0.8690 -0.6720];
        g = @(x) x.^3-2*x; b = 2; tikslus = 0;
    case 2
        f = @(x) sin(x); x0 = pi/3;
        tikros = [0.5 -sqrt(3)/2];
        x = [0.8472 0.9472 1.0472 1.1472 1.2472];
        y = [0.7494 0.8118 0.8660 0.9116 0.9481];
        g = @(x) sin(x); b = pi/2; tikslus = 1;
    case 3
        f = @(x) exp(x); x0 = 0;
        tikros = [1 1];
        % Lenteles duomenys atitinka exp(-x), ne exp(x).
        x = [-0.2 -0.1 0 0.1 0.2];
        y = [1.2214 1.1052 1 0.9048 0.8187];
        g = @(x) exp(-x); b = 1; tikslus = 1-exp(-1);
    otherwise
        error('Variantas turi buti 1, 2 arba 3.');
end
fprintf('Variantas: %d\n',variantas);
% f - diferencijavimo funkcija; g - integravimo funkcija.
h_reiksmes = [0.1 0.01];
paklaidos = zeros(2,4);
for i = 1:2
    h = h_reiksmes(i);
    yf = f(x0+(-2:2)*h);
    d = isvestines(yf,h);
    paklaidos(i,:) = abs(d-tikros([1 1 2 2]));
end
disp('Analitines funkcijos isvestiniu absoliucios paklaidos:');
disp(array2table([h_reiksmes' paklaidos],'VariableNames', ...
    {'h','First3','First5','Second3','Second5'}));

% Lenteles duomenu diferencijavimas.
d = isvestines(y,x(2)-x(1));
p2 = lagran(x(2:4),y(2:4));
p4 = lagran(x,y);
dp = [polyval(polyder(p2),x(3)) polyval(polyder(p4),x(3)) ...
    polyval(polyder(polyder(p2)),x(3)) polyval(polyder(polyder(p4)),x(3))];
disp('Lenteles duomenu isvestines:');
disp(array2table([d; dp],'VariableNames', ...
    {'First3','First5','Second3','Second5'}, ...
    'RowNames',{'Finite differences','Lagrange'}));

% Baigtinis integralas: tas pats destytojo paskirtas variantas.
Ns = [4 8]; paklaidos = zeros(2,3);
for i = 1:2
    N = Ns(i);
    I = [trpzds(g,0,b,N) smpsns(g,0,b,N) ...
        rmbrg(g,0,b,0,log2(N)+1)];
    paklaidos(i,:) = abs(I-tikslus);
end
disp('Baigtinio integralo absoliucios paklaidos:');
disp(array2table([Ns' paklaidos],'VariableNames', ...
    {'N','Trapezoid','Simpson','Romberg'}));

% Netiesioginis integralas sin(x)/x nuo 0 iki begalybes.
if exist('syms','file') == 2
    syms z
    I_simbolinis = int(sin(z)/z,z,0,inf);
    disp('Simbolinis integralas:'); disp(I_simbolinis);
    I_tikslus = double(I_simbolinis);
else
    I_tikslus = pi/2;
    disp('Nera Symbolic Math Toolbox: int nevykdomas; zinoma reiksme pi/2.');
end
f = @(x) sin(x)./x;
a = 0.001; b = 100; N = 200; tol = 1e-4; IGL = 20;
I = [smpsns(f,a,b,N); adapt_smpsn(f,a,b,tol); ...
    quad(f,a,b,tol); quadl(f,a,b,tol); Gauss_Legendre(f,a,b,IGL)];
% Atskiriame skaitine paklaida nuo nukirstos integralo dalies.
I_baigtinis = integral(f,a,b,'AbsTol',1e-12,'RelTol',1e-12);
metodai = {'Simpson'; 'Adaptive Simpson'; 'quad'; 'quadl'; 'Gauss-Legendre'};
disp(table(metodai,I,abs(I-I_tikslus),abs(I-I_baigtinis), ...
    'VariableNames',{'Method','Integral','ErrorToInfinite','ErrorOnFiniteInterval'}));

function d = isvestines(y,h)
% Keturios uzduoties formules: 3 ir 5 tasku isvestines.
d = [(y(4)-y(2))/(2*h), ...
    (-y(5)+8*y(4)-8*y(2)+y(1))/(12*h), ...
    (y(4)-2*y(3)+y(2))/h^2, ...
    (-y(5)+16*y(4)-30*y(3)+16*y(2)-y(1))/(12*h^2)];
end
