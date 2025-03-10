% script to construct a viscous flow model of a landslide with friction BC
% this solution is a steady-state solution (static BEM)
% 
% AUTHORS
% Rishav Mallick, JPL, Caltech
% March 2025

clear

addpath functions/
import('geometry.*')
% assume landlside properties
rhog = 1; % landslide density - ρg
theta = 40;% slop in degrees
% load landslide geometry from file
eM = geometry.LDhs(1,0.25);% dummy parameters
load('mesh/rcv.mat','rcv');
Lx = 5; % horizontal dimension of landslide

% plot landslide geometry
figure(1),clf
plotpatch2d(rcv)
quiver(rcv.xc(:,1),rcv.xc(:,2),rcv.nv(:,1),rcv.nv(:,2),0.2,'k')
quiver(rcv.xc(:,1),rcv.xc(:,2),rcv.dv(:,1),rcv.dv(:,2),0.2,'r')
plot(rcv.xc(:,1),rcv.xc(:,2),'kx')
axis tight equal

%% compute potential and potential-gradient kernels
[Ku,KDu] = compute_disp_stress_kernels_force(rcv,rcv.xc);
Kdudx = KDu(:,:,1);
Kdudz = KDu(:,:,2);

%% first solve pressure BC problem Δp = 0
% at the top boundary p = 0
% at the bottom boundary dp/dn = -ρg(nx.sinθ - nz.cosθ)
% the pressure flux is outward facing for the mesh (see Figure 1)
% for bottom boundary nv = [0,-1] so dp/dn =  -ρgcosθ
% for left boundary   dp/dn =  0
% BC labels are stored in rcv.BClabel 
% 0 - velocity BC (left boundary)
% 1 - traction BC (top boundary)
% 2 - mixed BC (bottom boundary)
BCvec = zeros(rcv.N,1);
index = rcv.BClabel==2;
BCvec(index) = -rhog*(rcv.nv(index,1).*sind(theta) - rcv.nv(index,2).*cosd(theta));
index = rcv.BClabel==0; % for left boundary
BCvec(index) = 0;

% construct flux kernel K = nx * K,x + nz * K,z such that K * φ = q
nxmat = repmat(rcv.nv(:,1),1,rcv.N);
nzmat = repmat(rcv.nv(:,2),1,rcv.N);
Kflux = Kdudx.*nxmat + Kdudz.*nzmat;

% construct kernel - flux and pressure BC
kernel = zeros(size(Ku));
index = rcv.BClabel~=1; % for bottom and left boundary
kernel(~index,:) = Ku(~index,:); % pressure BC at top
kernel(index,:) = Kflux(index,:); %dp/dn BC at bottom

% solve BEM problem
phi = kernel\BCvec;

%% compute solution inside the domain
nxo = 1000;
nzo = nxo/Lx;
xovec = linspace(1e-6,Lx-1e-6,nxo);
zovec = linspace(-0.1,1,nzo);
[xo,zo] = meshgrid(xovec,zovec);
% compute potential field kernel
[Ku_o,KDu_o] = compute_disp_stress_kernels_force(rcv,[xo(:),zo(:)]);
% compute pressure in the medium
p = Ku_o * phi;
dpdx = KDu_o(:,:,1)*phi;
dpdz = KDu_o(:,:,2)*phi;

% index points inside landslide
plotindex = inpolygon(xo,zo,rcv.x(:,1),rcv.x(:,2));
% plot pressure field
figure(2),clf
toplot = reshape(p,nzo,nxo);
toplot(~plotindex) = nan;
pcolor(xovec,zovec,toplot), shading interp, hold on
contour(xovec,zovec,toplot,-1:0.1:0,'k-','LineWidth',0.1)
plotpatch2d(rcv)
axis tight equal, box on
cb=colorbar;cb.Location='eastoutside';
cb.Label.Interpreter='latex';cb.Label.String='$\frac{p}{\rho g}$';
cb.LineWidth=1;cb.Label.Rotation=0;
clim([-1,0])
colormap(hot(10))
xlabel('x'), ylabel('z')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

figure(3),clf
subplot(2,1,1)
toplot = reshape(dpdx,nzo,nxo);
toplot(~plotindex) = nan;
pcolor(xovec,zovec,toplot), shading interp, hold on
contour(xovec,zovec,toplot,0:0.1:1,'k.-','LineWidth',0.1)
plotpatch2d(rcv)
axis tight equal
cb=colorbar;cb.Location='eastoutside';
cb.Label.Interpreter='latex';cb.Label.String='$\frac{1}{\rho g}\frac{\partial p}{\partial x}$';
cb.LineWidth=1;cb.Label.Rotation=0;
clim([0,1])
xlabel('x'), ylabel('z')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')

subplot(2,1,2)
toplot = reshape(dpdz,nzo,nxo);
toplot(~plotindex) = nan;
pcolor(xovec,zovec,toplot), shading interp, hold on
contour(xovec,zovec,toplot,0:0.1:1,'k.-','LineWidth',0.1)
plotpatch2d(rcv)
axis tight equal
cb=colorbar;cb.Location='eastoutside';
cb.Label.Interpreter='latex';cb.Label.String='$\frac{1}{\rho g}\frac{\partial p}{\partial z}$';
cb.LineWidth=1;cb.Label.Rotation=0;
clim([0,1]*1)
colormap(turbo(10))
xlabel('x'), ylabel('z')
set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','both')


