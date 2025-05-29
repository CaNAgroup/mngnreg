function J = jack(funz,x,dx)
%JACK approximate the Jacobian matrix by finite differences
%   J = jack(funz,x,dx)
%   where
%      funz - function to compute the Jacobian
%      x    - point at which the Jacobian is evaluated
%      dx   - step
%   and
%      J    - Jacobian matrix of funz

%   G.P. Deidda, P. Díaz de Alba, C. Fenu, and G. Rodriguez
%   University of Cagliari, Italy
%
%   Last revised June 10, 2019

x=x(:);
n = length(x);

if nargin < 3 || isempty(dx) , dx = 1e-3; end

y = funz(x);

m = length(y);

J = zeros(m,n);

for j = 1:n
	x_dx = x;
	x_dx(j) = x(j) + dx;
	y_dx = funz(x_dx);
	J(:,j) = (y_dx-y)./dx;
	%[y y_dx y_dx-y J(:,j)], pause
end


