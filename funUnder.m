function [F, J, H] = funUnder(x)
%FUNUNDER generates a simple underdetermined nonlinear problem
%   [F, J, H] = funUnder(x)
%   where
%      x - input vector
%   and
%      F - nonlinear function
%      J - Jacobian matrix
%      H - Hessian matrix

%   F. Pes and G. Rodriguez
%   University of Cagliari, Italy

% Last revised April 4, 2025

p = 1/9; q = 1/9;
% p = 1; q = 1/10;

F = ((p*(x(1)-1).^2 + q*(x(2)-1).^2 - 1).^2);

if nargout>1
	J = [4*p*(x(1)-1)*(p*(x(1)-1).^2+q*(x(2)-1).^2 -1), 4*q*(x(2)-1)*(p*(x(1)-1).^2+q*(x(2)-1).^2 -1)];
end

if nargout>2
	H = [12*p^2*(x(1)-1).^2+4*p*q*(x(2)-1).^2-4*p 8*p*q*(x(1)-1).*(x(2)-1); ...
		8*p*q*(x(1)-1).*(x(2)-1) 4*p*q*(x(1)-1).^2+12*q^2*(x(2)-1).^2-4*q];
end
