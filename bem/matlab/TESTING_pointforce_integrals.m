% script to use numerical quadrature to compute response to body force
% distribution inside the landslide body evaluated inside the landslide AND
% at the boundary.
% 
% AUTHORS:
% Rishav Mallick, JPL, Caltech
% March 2025

clear

% load landslide geometry from file
eM = geometry.LDhs(1,0.25);% dummy parameters
mesh = geometry.shearZoneReceiver('mesh/geometry',eM);
load('mesh/rcv.mat','rcv');
forcedata = importdata('mesh/sol_bodyforce.dat');

% compute kernels evaluated at interior points
[Gu,Strain] = compute_disp_stress_kernels_forcevolume(mesh,mesh.xc);
Gux = Strain(:,:,1);
Guz = Strain(:,:,2);

% compute kernels evaluated at the boundary
[Gu_r,Strain] = compute_disp_stress_kernels_forcevolume(mesh,rcv.xc);
Gux_r = Strain(:,:,1);
Guz_r = Strain(:,:,2);

% plot random realization of the kernel
% L = 0.3;
% r = sqrt((mesh.xc(:,1)-0.5).^2 + (mesh.xc(:,2)-0.3).^2);
% phi = exp(-(r.^2)./(2*L^2));

% body force vector from file
phi_x = forcedata(:,3);
phi_z = forcedata(:,4);

figure(11),clf
% plot potential field: sqrt(Vx^2 + Vz^2)
subplot(4,1,1)
plotshz2d(mesh,sqrt((Gu*phi_x).^2 + (Gu*phi_z).^2)), hold on
plotpatch2d(rcv,sqrt((Gu_r*phi_x).^2 + (Gu_r*phi_z).^2))
quiver(mesh.xc(:,1),mesh.xc(:,2),Gu*phi_x,Gu*phi_z,0,'w')
quiver(rcv.xc(:,1),rcv.xc(:,2),Gu_r*phi_x,Gu_r*phi_z,0,'k','LineWidth',2)
axis tight equal, box on
cb=colorbar;cb.Label.Interpreter='latex';cb.Label.String='$v$';
clim([0 1]*0.3)

% plot dv_x/dx = G,x.φx
subplot(4,1,2)
plotshz2d(mesh,Gux*phi_x), hold on
plotpatch2d(rcv,Gux_r*phi_x)
axis tight equal, box on
clim([-1 1]*0.2)
cb=colorbar;cb.Label.Interpreter='latex';cb.Label.String='$\frac{\partial v_x}{\partial x}$';

% plot dv_z/dz = G,z.φz
subplot(4,1,3)
plotshz2d(mesh,Guz*phi_z), hold on
plotpatch2d(rcv,Guz_r*phi_z)
axis tight equal, box on
clim([-1 1]*0.2)
cb=colorbar;cb.Label.Interpreter='latex';cb.Label.String='$\frac{\partial v_z}{\partial z}$';

% plot dv_x/dz + dv_z/dx = G,z.φx + G,x.φz
subplot(4,1,4)
plotshz2d(mesh,Guz*phi_x + Gux*phi_z), hold on
plotpatch2d(rcv,Guz_r*phi_x + Gux_r*phi_z)
axis tight equal, box on
cb=colorbar;cb.Label.Interpreter='latex';cb.Label.String='$\frac{\partial v_x}{\partial z}+\frac{\partial v_z}{\partial x}$';
clim([-1 1]*0.2)
colormap(turbo(100))
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20,'Linewidth',1.5)
