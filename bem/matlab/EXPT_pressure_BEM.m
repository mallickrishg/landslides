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
% at the bottom boundary p = -ρgh(x)
% at the left boundary p = -ρg(h-z)
% BC labels are stored in rcv.Vpl 
% 0 - traction BC
% 1 - velocity BC
% 2 - mixed BC
BCvec = zeros(rcv.N,1);
BCvec(rcv.Vpl==2) = -rhog*flipud(rcv.xc(rcv.Vpl==0,2));
BCvec(rcv.Vpl==1) = -rhog*(1-rcv.xc(rcv.Vpl==1,2));

nxmat = repmat(rcv.nv(:,1),1,rcv.N);
nzmat = repmat(rcv.nv(:,2),1,rcv.N);
% construct flux kernel K = nx * K,x + nz * K,z such that K * φ = q
Kflux = Kdudx.*nxmat + Kdudz.*nzmat;
% for pressure problem, BC is simply p = ρgh
kernel = zeros(size(Ku));
index = rcv.Vpl~=1;
kernel(index,:) = Ku(index,:);
kernel(~index,:) = Ku(~index,:);

% if left boundary has flux BC
% BCvec(2*npts+1:end) = 0; % if this dp/dx = 0 at left boundary
% kernel(~index,:) = Kflux(~index,:);

% solve BEM problem
phi = kernel\BCvec;

%% compute solution inside the domain
nxo = 100;
nzo = 20;
xovec = linspace(1e-6,Lx-1e-6,nxo);
zovec = linspace(0,1,nzo);
[xo,zo] = meshgrid(xovec,zovec);
% compute potential field kernel
[Ku_o,KDu_o] = compute_disp_stress_kernels_force(rcv,[xo(:),zo(:)]);
% compute pressure in the medium
p = Ku_o * phi;
dpdx = KDu_o(:,:,1)*phi;
dpdz = KDu_o(:,:,2)*phi;

% plot pressure field
figure(2),clf
pcolor(xovec,zovec,reshape(p,nzo,nxo)), shading interp, hold on
contour(xovec,zovec,reshape(p,nzo,nxo),-1:0.1:0,'k--','LineWidth',0.1)
plotpatch2d(rcv)
axis tight equal, box on
cb=colorbar;cb.Location='eastoutside';cb.Label.String='p';
clim([-1,0])
colormap(hot(10))

figure(3),clf
subplot(2,1,1)
pcolor(xovec,zovec,reshape(dpdx,nzo,nxo)), shading interp, hold on
contour(xovec,zovec,reshape(dpdx,nzo,nxo),-1:0.1:1,'k.-','LineWidth',0.1)
plotpatch2d(rcv)
axis tight equal
cb=colorbar;cb.Location='eastoutside';cb.Label.String='dp/dx';
clim([-1,1])
subplot(2,1,2)
pcolor(xovec,zovec,reshape(dpdz,nzo,nxo)), shading interp, hold on
contour(xovec,zovec,reshape(dpdz,nzo,nxo),0:0.05:1,'k.-','LineWidth',0.1)
plotpatch2d(rcv)
axis tight equal
cb=colorbar;cb.Location='eastoutside';cb.Label.String='dp/dz';
clim([0,1]*1)
colormap(turbo(20))
