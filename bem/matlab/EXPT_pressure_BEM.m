% script to compute the pressure field inside a landslide geometry subject
% to mixed boundary conditions:
% (1) at the top we impose p = 0
% (2) on the left side of the domain dp/dn =  0
% (3) at the bottom dp/dn = -ρg(nx.sinθ - nz.cosθ)
% _____________________________________________________________________________________________
%                                   EXPORT 
% _____________________________________________________________________________________________
% pressure at the boundary: p(x,z) 
% body force vector 
% fx(x,z) = dp/dx + ρg.sinθ, 
% fz(x,z) = dp/dz - ρg.cosθ
% 
% AUTHORS:
% Rishav Mallick, JPL, Caltech
% March 2025

clear

addpath functions/
import('geometry.*')
% assume landlside properties
rhog = 1; % landslide density - ρg
theta = 30;% slop in degrees
% load landslide geometry from file
eM = geometry.LDhs(1,0.25);% dummy parameters
load('mesh/rcv.mat','rcv');
mesh = geometry.shearZoneReceiver('mesh/geometry',eM);

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
% BCvec(index) = -rhog*(rcv.nv(index,1).*sind(theta) - rcv.nv(index,2).*cosd(theta));

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
nxo = 500;
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

%% compute solution on triangular mesh
% compute potential field kernel
[Ku_o,KDu_o] = compute_disp_stress_kernels_force(rcv,mesh.xc);
% compute pressure in the medium
p_mesh = Ku_o * phi;
% compute body force for equilibrium equations
dpdx_mesh = KDu_o(:,:,1)*phi + rhog*sind(theta);
dpdz_mesh = KDu_o(:,:,2)*phi - rhog*cosd(theta);

figure(4),clf
subplot(3,1,1)
toplot = reshape(p,nzo,nxo);
toplot(~plotindex) = nan;
plotshz2d(mesh,p_mesh), hold on, shading flat
contour(xovec,zovec,toplot,-1:0.1:0,'k-','LineWidth',0.1)
plotpatch2d(rcv)
axis tight equal, box on
cb=colorbar;cb.Location='eastoutside';
cb.Label.Interpreter='latex';cb.Label.String='$p$';
cb.LineWidth=1;%cb.Label.Rotation=0;
clim([-1,0])

subplot(3,1,2)
toplot = reshape(dpdx,nzo,nxo);
toplot(~plotindex) = nan;
plotshz2d(mesh,dpdx_mesh), hold on, shading flat
% contour(xovec,zovec,toplot,-1:0.1:1,'k.-','LineWidth',0.1)
plotpatch2d(rcv)
axis tight equal, box on
cb=colorbar;cb.Location='eastoutside';
cb.Label.Interpreter='latex';cb.Label.String='$\frac{\partial p}{\partial x} + \rho g\sin\theta$';
cb.LineWidth=1;%cb.Label.Rotation=0;
clim([-1,1])

subplot(3,1,3)
toplot = reshape(dpdz,nzo,nxo);
toplot(~plotindex) = nan;
plotshz2d(mesh,dpdz_mesh), hold on, shading flat
% contour(xovec,zovec,toplot,0:0.1:1,'k.-','LineWidth',0.1)
plotpatch2d(rcv)
axis tight equal, box on
cb=colorbar;cb.Location='eastoutside';
cb.Label.Interpreter='latex';cb.Label.String='$\frac{\partial p}{\partial z}- \rho g\cos\theta$';
cb.LineWidth=1;%cb.Label.Rotation=0;
clim([-1,1]*1)
% mycolors = [1 0 0; 1 1 0; 0 0 1];
colormap(turbo(20))

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20,'Linewidth',1.5)

%% export solutions
% (1) pressure 'p' evaluated on the boundary, 
% (2) body force dp/dx + ρg.sinθ, dp/dz - ρg.cosθ evaluated at mesh.xc
[Ku,~] = compute_disp_stress_kernels_force(rcv,rcv.xc);
Trcv = table(rcv.xc(:,1),rcv.xc(:,2),Ku*phi);
Tmesh = table(mesh.xc(:,1),mesh.xc(:,2),dpdx_mesh,dpdz_mesh);

writetable(Trcv,'mesh/sol_pressure.dat','WriteVariableNames',false)
writetable(Tmesh,'mesh/sol_bodyforce.dat','WriteVariableNames',false)


