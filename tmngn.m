function [x, k, rho, fail, X, Res, Alphas, Betas, ells] = tmngn( fun, b, x0, ell, opts)
%TMNGN Gauss-Newton method regularized by TSVD
%   [x, k, rho, fail, X, Res, alpha, beta, ells] = tmngn( fun, b, x0, ell, opts)
%   computes the minimal-norm solution to the nonlinear least-square
%   problem min||fun(x)-b||,
%   where
%      fun   - nonlinear function and Jacobian matrix: [f,J]=fun(x)
%      b     - data vector
%      x0    - initial point
%      ell   - either rank or truncation parameter
%      opts  - options (struct, every options has a default)
%              niter: max number of iterations
%              damped: 1 damped method, 0 undamped
%              dampos: 1 positive solution, 0 unconstrained
%              findiff: finite differences to approximate the Jacobian
%              rankflag: to estimate the rank of the Jacobian
%              tol: tolerance for Jacobian rank
%              alphamin: minimum value for alpha and beta
%              mnflag: minimal-norm solution
%              xbar: a priori estimate of solution (only if mnflag)
%              tau: stopping tolerance
%              eta1: tolerance for residual increase
%              eta2: tolerance for residual increase
%              kres: number of residuals to determine eta2
%              flagdiffrep: number of convergence conditions to fulfill
%   and
%      x      - solution
%      k      - number of iterations performed
%      rho    - solution residual
%      fail   - failure (>1) or success (0 or 1)
%      X      - solution at each iteration
%      Res    - residual at each iteration
%      Alphas - step lengths for damped Gauss-Newton method
%      Betas  - step lengths for projection vector
%      ells   - numerical rank at each iteration

%   F. Pes and G. Rodriguez
%   University of Cagliari, Italy

% Last revised April 4, 2025

if nargin<4 || isempty(ell), ell = inf; end

nam = fieldnames(opts);

% maximum number of iterations
if ~any(strcmpi('niter',nam)), opts.niter = 100; end
% damped method
if ~any(strcmpi('damped',nam)), opts.damped = 1; end
%  damping parameter is used to keep the solution positive
if ~any(strcmpi('dampos',nam)), opts.dampos = 0; end
% finite differences to approximate the Jacobian
if ~any(strcmpi('findiff',nam)), opts.findiff = 0; end
% rankflag
if ~any(strcmpi('rankflag',nam))
	if isfinite(ell)
		opts.rankflag = 0;
	else
		opts.rankflag = 1;
	end
end
% tolerance for Jacobian rank
if ~any(strcmpi('tol',nam)), opts.tol = sqrt(eps); end
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
% number of residuals to determine eta2
if ~any(strcmpi('kres',nam)), opts.kres = 5; end
% number of convergence conditions to fulfill
if ~any(strcmpi('flagdiffrep',nam)), opts.flagdiffrep = 1; end

if opts.mnflag==0 && ~(isequal(opts.xbar,zeros(size(x0))))
	warning('xbar is useful only if mnflag>0.')
end

if ((opts.mnflag==5) || (opts.mnflag==6)) && (opts.damped==1)
	opts.damped=0;
	warning('using undamped method')
end

if isfinite(ell) && opts.rankflag
	opts.rankflag = 0;
	warning('opts.rankflag set to 0.')
end

eta1 = opts.eta1;
eta2 = opts.eta2;

n = size(x0,1);
m = size(b,1);
if ~isfinite(ell)		% max value for ell
	ell = min(m,n);
end

X = zeros(n,opts.niter+1);
Res = zeros(opts.niter+1,1);
Alphas = zeros(opts.niter,1);
Betas = zeros(opts.niter,1);
ells = zeros(opts.niter,1);

normx0 = norm(x0);
X(:,1) = x0;
rGN = fun(X(:,1))-b;
rho = norm(rGN);
Res(1) = rho;

k = 0;
flag = 1;
beta = 1;
flagalpha = 1;
kflag = 0;
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

	[U,S,V] = svd(jacGN);
	if min(m,n)==1
		d = S(1,1);
	else
		d = diag(S);
	end

	if (min(m,n)>1) 	% estimating rank of Jacobian
		switch opts.rankflag
		case 0
			;	% ell has already been set
		case 1
			ratio = d(1:end-1)./d(2:end);
			vv = find((ratio>1e2) & (d(1:end-1)>opts.tol));
			if isempty(vv)
				ell = min(m,n);
			else
				[~, ii] = max(ratio(vv));
				ell = vv(ii);
			end
		case 2
			vv = find(d<opts.tol);
			if isempty(vv)
				ell = min(m,n);
			else
				ell = vv(1)-1;
			end
		otherwise
			error('wrong rankflag.')
		end
	end
	ells(k) = ell;

	if ell<n
		V2 = V(:,ell+1:n);
	else
		V2 = 0;
	end
	U = U(:,1:ell);
	V = V(:,1:ell);
	d = d(1:ell);

	% step = - pinv(jacGN)*rGN;
	y = U'*rGN;
	y = -y./d;
	step = V*y;
	if opts.mnflag
		tk = V2*(V2'*(X(:,k)-opts.xbar));
	end

	if opts.mnflag == 2
		step = step - tk;
	end

	% damped
	if opts.damped
		alpha = 2;
		flagd = 1;
		while flagd
			alpha = alpha/2;
			x1 = X(:,k) + alpha * step;
			rho1 = norm(fun(x1)-b);
			alphas = norm(alpha*step);
			% Armijo-Goldstein principle
			njs2 = norm(jacGN*step)^2;
			flagd = (rho^2 - rho1^2) < (1/2)*alpha*njs2;
			if opts.dampos, flagd = flagd || any(x1<0); end
			flagalpha = alphas > opts.alphamin;
			flagd = flagd && flagalpha;
		end
		Alphas(k) = alpha;
		% x1 = X(:,k) + alpha*step; % already computed
	else
		x1 = X(:,k) + step;
	end

	switch opts.mnflag	% minimal-norm and step length beta
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

	rho = norm(fun(x2)-b);
	Res(k+1) = rho;

	diff = norm(X(:,k+1)-X(:,k)); % norm(x_k+1 - x_k)

	flagdiff = (diff > opts.tau*max(1,normx)); % tolerance = 1e-8
	if flagdiff
		kflag = 0;
	else
		kflag = kflag+1;
		if kflag < opts.flagdiffrep, flagdiff = 1; end
	end
	flagiter = (k < opts.niter);
	flagnorm = (normx < 1e8*normx0);
	%[flagdiff flagiter flagnorm flagalpha]
	flag = flagdiff && flagiter && flagnorm && flagalpha;

end

fail = flagdiff * ((~flagalpha)+10*(~flagiter)+100*(~flagnorm));
x = X(:,k+1);
X = X(:,1:k+1);
Res = Res(1:k+1);
Alphas = Alphas(1:k);
Betas = Betas(1:k);
ells = ells(1:k);

