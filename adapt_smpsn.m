function [INTf,nodes,err] = adapt_smpsn(f,a,b,tol,varargin)
% Adaptive_Simpsons.m from Int.zip, with short initial intervals.
% Initial splitting prevents missing the oscillations of sin(x)/x.
segments = max(1,ceil(abs(b-a)));
edges = linspace(a,b,segments+1);
INTf = 0; err = 0; nodes = a;
for i = 1:segments
    left = edges(i); right = edges(i+1);
    whole = smpsns(f,left,right,2,varargin{:});
    [value,partNodes,partError] = adap_smpsn(f,left,right,whole,tol/segments,0,varargin{:});
    INTf = INTf+value;
    err = err+partError;
    nodes = [nodes partNodes(2:end)]; %#ok<AGROW>
end
end

function [INTf,nodes,err] = adap_smpsn(f,a,b,INTf,tol,depth,varargin)
% The supplied recursive Simpson comparison and tolerance splitting.
c = (a+b)/2;
INTf1 = smpsns(f,a,c,2,varargin{:});
INTf2 = smpsns(f,c,b,2,varargin{:});
INTf12 = INTf1+INTf2;
err = abs(INTf12-INTf)/15;
if err <= tol
    INTf = INTf12+(INTf12-INTf)/15;
    nodes = [a c b];
elseif depth >= 30
    error('Adaptive Simpson reached the recursion limit.');
else
    [INTf1,nodes1,err1] = adap_smpsn(f,a,c,INTf1,tol/2,depth+1,varargin{:});
    [INTf2,nodes2,err2] = adap_smpsn(f,c,b,INTf2,tol/2,depth+1,varargin{:});
    INTf = INTf1+INTf2;
    nodes = [nodes1 nodes2(2:end)];
    err = err1+err2;
end
end
