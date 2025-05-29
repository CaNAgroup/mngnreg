%DRIVEROVER test program for tmlngn.m

%   F. Pes and G. Rodriguez
%   University of Cagliari, Italy

% Last revised April 4, 2025

fun = @(x) funOver(x);
b = [0; 0; 0];
x0 = [2;1.2]; % initial point
n = length(x0);
der = 1;
L = full(get_l(n,der));
ell = inf;
niter = 40; % max number of iterations
rankflag = 1;
opts = struct( 'niter', niter, ...
	'rankflag', rankflag );

opts.mnflag = 0;
[xLGN, k1, rho1, fail1, X1, Res1] = tmlngn( fun, b, L, x0, ell, opts);
opts.mnflag = 1;
[xMLNGN, k2, rho2, fail2, X2, Res2] = tmlngn( fun, b, L, x0, ell, opts);

nu = 41;
nv = 41;
u = linspace(0,2,nu)';
v = linspace(1,3.5,nv)';
ff = @(x) norm(funOver(x))^2/2;
F = zeros(nv,nu);
resGN = zeros(size(X1,2),1);
resMNGN = zeros(size(X2,2),1);
for i = 1:nu
	for j = 1:nv
		F(j,i) = ff([u(i);v(j)]);
	end
end
for i = 1:size(X1,2)
	resGN(i) = ff([X1(1,i);X1(2,i)]);
end
for i = 1:size(X2,2)
	resMNGN(i) = ff([X2(1,i);X2(2,i)]);
end

figure(1)
surf(u,v,F,'FaceColor','interp','EdgeColor','none')
colormap jet
hold on
plot3(X1(1,:),X1(2,:),resGN,'k.-','MarkerSize',15)
plot3(X2(1,:),X2(2,:),resMNGN,'r.-','MarkerSize',15)
hold off
legend('f','GN','MNGN','Location','Best')
view([-80 40])

figure(2)
minf = min(F(:));
maxf = max(F(:));
contour(u,v,F,linspace(minf,maxf,50))
hold on
plot(X1(1,:),X1(2,:),'k.-','MarkerSize',15)
plot(X2(1,:),X2(2,:),'r.-','MarkerSize',15)
hold off
colormap jet
colorbar
legend('contour','GN','MNGN')
grid
