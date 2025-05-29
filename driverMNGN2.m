%DRIVERMNGN2 test program for tmngn.m

%   F. Pes and G. Rodriguez
%   University of Cagliari, Italy

% Last revised April 4, 2025

n = 3; m = 2;
a = ones(n,1); % sphere or ellipsoid
c = [2;zeros(n-1,1)]; % center of the sphere
fun = @(x) nonlinfun1(x,a,c,m);
b = zeros(m,1);
x0 = [0;3;3]; % initial point

ell = inf;
niter = 60; % max number of iterations
opts4 = struct( 'niter', niter );
opts5 = struct( 'niter', niter );

opts4.mnflag = 4;	% MNGN2
[sol4, k4, rho4, fail4, X4, Res4, al4, be4] = tmngn( fun, b, x0, ell, opts4);

opts5.mnflag = 5;	% CKB
[sol5, k5, rho5, fail5, X5, Res5, al5, be5] = tmngn( fun, b, x0, ell, opts5);


N = 51;
[X1, Y1, Z1] = sphere(N);
X1 = X1+2;

xv = linspace(-1,3,N);
[Y2, Z2] = meshgrid(xv,xv);
X2 = 2*ones(N);

xv = linspace(-1,3,N);
[X3, Z3] = meshgrid(xv,xv);
Y3 = zeros(N);

figure(1)
plot3(X4(1,:),X4(2,:),X4(3,:),'.-b','markersize',12)
hold on
plot3(X5(1,:),X5(2,:),X5(3,:),'.--r','markersize',12)
plot3(1,0,0,'ok','markersize',8)
surfl(X1,Y1,Z1)
surfl(X2,Y2,Z2)
surfl(X3,Y3,Z3)
hold off
%title('nonlinfun1 - c=[2;0;0]')
legend('MNGN2','CKB','min-norm sol')
axis equal
shading interp
alpha(0.6)
grid

fprintf('\n\n')
fprintf('locus of the solutions: sphere + line intersection of the 2 planes\n')
fprintf('min-norm sol: %.2g %.2g %.2g\n', [1 0 0])
fprintf('MNGN2       : %.2g %.2g %.2g\n', sol4)
fprintf('CKB         : %.2g %.2g %.2g\n', sol5)

