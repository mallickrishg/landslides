clear

addpath functions/
eta = 1;
a = 1;
dip = 90;

nx = 100;
ny = 100;
xvec = linspace(-2,2,nx);
yvec = linspace(-2,2,ny);
[xg,yg] = meshgrid(xvec,yvec);

[ux,uy,sxy,sxx,syy] = calc_disp_stress_forcestokes(xg(:),yg(:),a,dip,eta);

figure(1),clf
subplot(2,2,1)
toplot = reshape(ux(:,1),ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
axis tight equal
colorbar
clim([-1 1]*0.4)
title('u_x (f_x)')
subplot(2,2,2)
toplot = reshape(ux(:,2),ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
axis tight equal
colorbar
clim([-1 1]*0.2)
title('u_x (f_y)')
subplot(2,2,3)
toplot=reshape(uy(:,1),ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
axis tight equal
colorbar
clim([-1 1]*0.2)
title('u_y (f_x)')
subplot(2,2,4)
toplot=reshape(uy(:,2),ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
axis tight equal
colorbar
clim([-1 1]*0.5)
title('u_y (f_y)')
set(findall(gcf, '-property', 'Fontsize'), 'FontSize', 15,'Linewidth',1.5)
colormap(bluewhitered)

figure(2),clf
subplot(3,2,1)
toplot = reshape(sxy(:,1),ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
axis tight equal
colorbar
subplot(3,2,2)
toplot = reshape(sxy(:,2),ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
axis tight equal
colorbar
subplot(3,2,3)
toplot=reshape(sxx(:,1),ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
axis tight equal
colorbar
subplot(3,2,4)
toplot=reshape(sxx(:,2),ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
axis tight equal
colorbar
subplot(3,2,5)
toplot = reshape(syy(:,1),ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
axis tight equal
colorbar
subplot(3,2,6)
toplot = reshape(syy(:,2),ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
axis tight equal
colorbar
set(findall(gcf, '-property', 'CLim'), 'FontSize', 15,'Linewidth',1.5,'CLim',[-1,1]*0.5)
colormap(turbo(100))

% plot line plots
xpos = 1;
[ux,uy,~,~,~] = calc_disp_stress_forcestokes(xvec(:).*0+xpos,yvec(:),a,dip,eta);
figure(10),clf
subplot(2,1,1)
plot(yvec,ux,'LineWidth',2)
axis tight, grid on
subplot(2,1,2)
plot(yvec,uy,'LineWidth',2)
axis tight, grid on
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 15,'Linewidth',1.5)