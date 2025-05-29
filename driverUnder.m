%DRIVERUNDER test program for tmngn.m

%   F. Pes and G. Rodriguez
%   University of Cagliari, Italy

% Last revised April 4, 2025

fun = @(x) funUnder(x);
b = -1;
x0 = [5;3]; % initial point
n = length(x0);
ell = inf;
niter = 40; % max number of iterations
rankflag = 1;
opts = struct( 'niter', niter, ...
	'rankflag', rankflag );

opts.mnflag = 0;
[xGN, k1, rho1, fail1, X1, Res1] = tmngn( fun, b, x0, ell, opts);
opts.mnflag = 1;
[xMNGN, k2, rho2, fail2, X2, Res2] = tmngn( fun, b, x0, ell, opts);

nu = 41;
nv = 41;
u = linspace(-6,8,nu)';
v = linspace(-6,8,nv)';
ff = @(x) norm(funUnder(x)-b)^2/2;
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
zlim([0,4])
clim([0,4])
view([110 -40])
colormap jet
hold on
plot3(X1(1,:),X1(2,:),resGN,'k.-','MarkerSize',15)
plot3(X2(1,:),X2(2,:),resMNGN,'r.-','MarkerSize',15)
hold off
legend('f','GN','MNGN','Location','best')
xlim([-6,6])
ylim([-6,6])

figure(2)
minf = min(F(:));
maxf = max(F(:));
contour(u,v,F,logspace(log10(minf),log10(maxf),50))
clim([0 150])
hold on
plot(X1(1,:),X1(2,:),'k.-','MarkerSize',15)
plot(X2(1,:),X2(2,:),'r.-','MarkerSize',15)
hold off
colormap jet
colorbar
legend('contour','GN','MNGN')
xlim([-6 8])
ylim([-6 8])
axis equal
grid
