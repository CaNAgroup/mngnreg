function [x, k, rho, fail, X, Res, Alphas, Betas] = tiklgn( fun, b, L, x0, lam, opts)
%TIKLGN Gauss-Newton method regularized by Tikhonov in general form
%   [x, k, rho, fail, X, Res, alpha, beta] = tiklgn( fun, b, L, x0, lam, opts)
%   computes the regularized solution to the nonlinear least-square
%   problem min||fun(x)-b||, i.e.,
%   min { ||fun(x)-b||^2 + lam^2||L(x-xbar)||^2 }.
%   where
%      fun   - function and Jacobian: [f,J]=fun(x)
%      b     - data vector
%      L     - regularization matrix
%      x0    - initial point
%      lam   - regularization parameter
%      opts  - options (struct, every options has a default)
%              niter: max number of iterations
%              damped: 1 damped method, 0 undamped
%              dampos: 1 positive solution, 0 unconstrained
%              findiff: finite differences to approximate the Jacobian
%              alphamin: minimum value for alpha and beta
%              mnflag: minimal-norm solution
%              xbar: a priori estimate of solution (only if mnflag)
%              tau: stopping tolerance
%              eta1: tolerance for residual increase
%              eta2: tolerance for residual increase
%              kres: number of residuals to determine eta
%              curegsvd: 1 "cure" gsvd bugs, 0 no
%              factgsvd: small norm in gsvd
%   and
%      x      - solution
%      k      - number of iterations performed
%      rho    - solution residual
%      fail   - failure (>1) or success (0 or 1)
%      X      - solution at each iteration
%      Res    - residual at each iteration
%      Alphas - step-lengths for damped Gauss-Newton method
%      Betas  - step lengths for projection vector

%   F. Pes and G. Rodriguez
%   University of Cagliari, Italy

% Last revised April 4, 2025

if nargin<5 || isempty(lam), lam = 0; end
if lam<0, error('Illegal regularization parameter lambda.'), end

nam = fieldnames(opts);

% maximum number of iterations
if ~any(strcmpi('niter',nam)), opts.niter = 100; end
% damped method
if ~any(strcmpi('damped',nam)), opts.damped = 1; end
%  damping parameter is used to keep the solution positive
if ~any(strcmpi('dampos',nam)), opts.dampos = 0; end
% finite differences to approximate the Jacobian
if ~any(strcmpi('findiff',nam)), opts.findiff = 0; end
% minimum value for alpha and beta
if ~any(strcmpi('alphamin',nam)), opts.alphamin = 1e-9; end
% minimal norm solution
if ~any(strcmpi('mnflag',nam)), opts.mnflag = 4; end
% a priori estimate of solution (only if mnflag)
if ~any(strcmpi('xbar',nam)), opts.xbar = zeros(size(x0)); end
% stopping tolerance
if ~any(strcmpi('tau',nam)), opts.tau = 1e-8; end
% tolerance for residual increase
if ~any(strcmpi('eta1',nam)), opts.eta1 = 8; end
% tolerance for residual increase
if ~any(strcmpi('eta2',nam)), opts.eta2 = 1/8; end
% number of residuals to determine eta
if ~any(strcmpi('kres',nam)), opts.kres = 5; end
% cure gsvd bugs
if ~any(strcmpi('curegsvd',nam)), opts.curegsvd = 1; end
% small norm in gsvd
if ~any(strcmpi('factgsvd',nam)), opts.factgsvd = 1e-6; end

if opts.mnflag==0 && ~(isequal(opts.xbar,zeros(size(x0))))
	warning('xbar is useful only if mnflag>0.')
end

if ((opts.mnflag==5) || (opts.mnflag==6)) && (opts.damped==1)
	opts.damped=0;
	warning('using undamped method')
end

eta1 = opts.eta1;
eta2 = opts.eta2;

n = size(x0,1);
m = size(b,1);
p = size(L,1);
d = n-p;
mmn = min(m,n);

if (n~=size(L,2)) || (d<0), error('Incorrect dimensions for L.'), end

X = zeros(n,opts.niter+1);
Res = zeros(opts.niter+1,1);
Alphas = zeros(opts.niter,1);
Betas = zeros(opts.niter,1);

normLx0 = norm(L*x0);
X(:,1) = x0;
rGN = fun(X(:,1))-b;
rho = norm(rGN);
Res(1) = rho;

k = 0;
flag = 1;
beta = 1;
flagalpha = 1;
while flag
	k = k+1;

	if opts.findiff
		fGN = fun(X(:,k));
		dx = 1e-6;
		jacGN = jack(fun,X(:,k),dx);
	else
		[fGN,jacGN] = fun(X(:,k));
	end
	rGN = fGN-b;

	smalljac = 0;
	if opts.curegsvd && (norm(jacGN,inf) < opts.factgsvd)
		warning('small Jacobian!')
		smalljac = 1;
		[U,V,WI,C,S] = gsvd(jacGN/opts.factgsvd,L);
	else
		[U,V,WI,C,S] = gsvd(jacGN,L);
	end
	WI = WI';	% WI=inv(W), with the formalism of the paper
	%W = inv(WI);	% check if it's necessary
	[LWI, UWI, PWI] = lu(WI);

	if mmn == 1
		c = C(1,n);
	else
		c = diag(C(1:mmn,n-mmn+1:n));
	end
	r = sum(c>0);
	c = c(mmn-r+1:mmn-d);		% excludes the ones
	%c = c(mmn-r+1:mmn);		% includes the ones - useless
	%W1 = W(:,1:n-r);
	%W2 = W(:,n-r+1:n);

	if p == 1
		s = S(1,1);
	else
		s = diag(S);
	end
	s = s(p-r+d+1:p);

	U = U(:,1:mmn);		% mmn = min(m,n)

	if smalljac
		gk = U'*(rGN/opts.factgsvd);
	else
		gk = U'*rGN;
	end

	if ~opts.mnflag
		y = (c.*gk(mmn-r+1:mmn-d))./(c.^2+lam^2*s.^2);
	end
	zk = WI*(X(:,k)-opts.xbar);
	if r<n
		%tk = W1*zk(1:n-r);
		rhs = [zk(1:n-r);zeros(r,1)];
		sol = LWI \ (PWI*rhs);
		tk = UWI \ sol;
	else
		tk = 0;
	end

	% damped
	if opts.damped
		alpha = 2;
		flagd = 1;
		while flagd
			alpha = alpha/2;
			if opts.mnflag
				y = (c.*gk(mmn-r+1:mmn-d) + ...
					alpha*lam^2*s.^2.*zk(n-r+1:p)) ./ ...
					(c.^2+alpha^2*lam^2*s.^2);
			end
			%step = -W2*[y;gk(mmn-d+1:mmn)];
			rhs = [zeros(n-r,1);[y;gk(mmn-d+1:mmn)]];
			sol = LWI \ (PWI*rhs);
			step = -UWI \ sol;
			if opts.mnflag == 2
				step = step - tk;
			end
			x1 = X(:,k) + alpha * step;
			rho1 = norm(fun(x1)-b);
			alphas = norm(alpha*step);
			% Armijo-Goldstein principle
			if smalljac
				njs2 = (norm(jacGN*step)/opts.factgsvd)^2;
			else
				njs2 = norm(jacGN*step)^2;
			end
			flagd = (rho^2 - rho1^2) < (1/2)*alpha*njs2;
			if opts.dampos, flagd = flagd || any(x1<0); end
			flagalpha = alphas > opts.alphamin;
			flagd = flagd && flagalpha;
		end
		Alphas(k) = alpha;
		% x1 = X(:,k) + alpha*step; % already computed
	else
		if opts.mnflag
			y = (c.*gk(mmn-r+1:mmn-d) + ...
				lam^2*s.^2.*zk(n-r+1:p)) ./ ...
				(c.^2+lam^2*s.^2);
		end
		%step = -W2*[y;gk(mmn-d+1:mmn)];
		rhs = [zeros(n-r,1);[y;gk(mmn-d+1:mmn)]];
		sol = LWI \ (PWI*rhs);
		step = -UWI \ sol;
		if opts.mnflag == 2
			step = step - tk;
		end
		x1 = X(:,k) + step;
	end

	switch opts.mnflag	% minimal-norm and step lenght beta
		case 0		% no minimal-norm
			x2 = x1;
		case 1		% beta=1
			x2 = x1 - tk;
		case 2		% 1 parameter, beta=alpha
			x2 = x1;
		case 3		% 2 parameters, constant eta
			rho1 = norm(fun(x1)-b);
			if beta<1, beta = beta*2; end
			x2 = x1 - beta*tk;
			rho2 = norm(fun(x2)-b);
			while (rho2>(1+eta1)*(rho1+eps)) && (beta>opts.alphamin)
				beta = beta/2;
				x2 = x1 - beta*tk;
				rho2 = norm(fun(x2)-b);
			end
		case 4		% 2 parameters, adaptive eta
			rho1 = norm(fun(x1)-b);
			if beta<1, beta = beta*2; end
			x2 = x1 - beta*tk;
			rho2 = norm(fun(x2)-b);
			if k>=opts.kres
				pol = polyfit([1:opts.kres]',...
					log10(Res(k-opts.kres+1:k)),1);
				if (pol(1)>-1e-2)
					if eta2<1e4, eta2 = 2*eta2; end
				elseif (pol(1)<-1/2)
					eta2 = eta2/2;
				end
			end
			while (rho2>rho1+(rho1+eps)^eta2) && (beta>opts.alphamin)
				beta = beta/2;
				x2 = x1 - beta*tk;
				rho2 = norm(fun(x2)-b);
			end
		case 5		% CKB1
			alpha = 0.5^k;
			x2 = x1 - alpha*tk;
		case 6		% CKB2
			alpha = 0.5^(2^(k-1));
			x2 = x1 - alpha*tk;
		otherwise
			error('Wrong mnflag value.')
	end
	Betas(k) = beta;
	X(:,k+1) = x2;
	normx = norm(X(:,k+1));
	normLx = norm(L*X(:,k+1));

	rho = norm(fun(x2)-b);
	Res(k+1) = rho;

	diff = norm(X(:,k+1)-X(:,k)); % norm(x_k+1 - x_k)

	flagdiff = (diff > opts.tau*max(1,normx)); % tolerance = 1e-8
	flagiter = (k < opts.niter);
	flagnorm = (normLx < 1e8*normLx0);
	%[flagdiff flagiter flagnorm flagalpha]
	flag = flagdiff && flagiter && flagnorm && flagalpha;
end

fail = flagdiff * ((~flagalpha)+10*(~flagiter)+100*(~flagnorm));
x = X(:,k+1);
X = X(:,1:k+1);
Res = Res(1:k+1);
Alphas = Alphas(1:k);
Betas = Betas(1:k);

