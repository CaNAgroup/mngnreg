function [F, J] = funOver(x)
%FUNOVER generates a simple overdetermined nonlinear problem
%   [F, J] = funOver(x)
%   where
%      x - input vector
%   and
%      F - nonlinear function
%      J - Jacobian matrix

%   F. Pes and G. Rodriguez
%   University of Cagliari, Italy

% Last revised April 4, 2025

F = [ x(1).^2 + x(2).^2 - 2*x(2) - 3;...
	x(1).^2 + x(2).^2 - 4*x(1) - 2*x(2) + 1;...
	x(1) - 1 ];

if nargout>1
	J = [ 2*x(1), 2*x(2)-2;...
		2*x(1)-4, 2*x(2)-2;...
		1, 0 ];
end
