function [F, J, S, y1, y2] = nonlinfun1(x,a,c,m)
%NONLINFUN generates a nonlinear problem
%   [F,J]=nonlinfun(x,a,c) assigns parameters a, c, and returns F(x) and J(x).
%
%   [F,Jf,S,y]=nonlinfun(x,a,c) assigns parameters a, c, and returns F(x)
%   and the function Jf such that Jf(x)=J(x); S and y are the scalar and
%   the vector which define the jacobian
%   J=S*I+2*y*y'.
%
%   [F,...]=nonlinfun(x,a,c,m), with m<=n(=length(x)), reduces the size of F(x)

%   F. Pes and G. Rodriguez
%   University of Cagliari, Italy

% Last revised April 4, 2025

n = size(x,1);
if nargin<4, m = n; end
if m>n, error('m too big.'), end

S = norm((x-c)./a)^2 - 1;
y2 = (1./a.^2).*(x-c);
x = x(1:m);
c = c(1:m);
y1 = x-c;
F = S*(x-c);
if nargout>2
	J = @(x) S*x(1:m)+2*(y2'*x)*y1;
else
	J = S*eye(m,n) + 2*y1*y2';
end
