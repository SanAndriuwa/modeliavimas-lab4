function c = difapx(order,points)
% Finite-difference weights for the integer offsets points(1):points(2).
offsets = points(1):points(2);
count = length(offsets);
if count < order+1
    error('More nodes are needed for this derivative.');
end
A = zeros(count);
for power = 0:count-1
    A(power+1,:) = offsets.^power/factorial(power);
end
rhs = zeros(count,1); rhs(order+1) = 1;
c = (A\rhs)';
end
