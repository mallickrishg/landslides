clear

% load landslide geometry from file
eM = geometry.LDhs(1,0.25);% dummy parameters
mesh = geometry.shearZoneReceiver('mesh/geometry',eM);
load('mesh/rcv.mat','rcv');
forcedata = importdata('mesh/sol_bodyforce.dat');

% compute kernels
[Gu,Strain] = compute_disp_stress_kernels_forcevolume(mesh,mesh.xc);
Gux = Strain(:,:,1);
Guz = Strain(:,:,2);

% compute kernels
[Gu_r,Strain] = compute_disp_stress_kernels_forcevolume(mesh,rcv.xc);
Gux_r = Strain(:,:,1);
Guz_r = Strain(:,:,2);

% plot random realization of the kernel
L = 0.3;
r = sqrt((mesh.xc(:,1)-0.5).^2 + (mesh.xc(:,2)-0.3).^2);
phi = exp(-(r.^2)./(2*L^2));
phi = forcedata(:,3);

figure(11),clf
subplot(3,1,1)
plotshz2d(mesh,Gu*phi), hold on
plotpatch2d(rcv,Gu_r*phi)
axis tight equal, box on
colorbar

subplot(3,1,2)
plotshz2d(mesh,Gux*phi), hold on
plotpatch2d(rcv,Gux_r*phi)
axis tight equal, box on
colorbar

subplot(3,1,3)
plotshz2d(mesh,Guz*phi), hold on
plotpatch2d(rcv,Guz_r*phi)
axis tight equal, box on
colorbar

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20,'Linewidth',1.5)
