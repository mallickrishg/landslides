clear

addpath functions/
eta = 1;
a = 0.5;
dip = 0; % atan2d(Δy,Δx) where Δy=y2-y1,Δx=x2-x1

nx = 100;
ny = 100;
xvec = linspace(-2,2,nx);
yvec = linspace(-2,2,ny);
[xg,yg] = meshgrid(xvec,yvec);

[ux,uy,sxy,sxx,syy] = calc_disp_stress_forcestokes(xg(:),yg(:),a,dip,eta);
[u,ex,ey] = calc_disp_stress_forcefault(xg(:),yg(:),a,dip);

nskip = 13;
figure(1),clf
subplot(2,2,1)
toplot = reshape(ux(:,1),ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
quiver(xg(1:nskip:end)',yg(1:nskip:end)',ux(1:nskip:end,1),uy(1:nskip:end,1),'k')
contour(xvec,yvec,toplot,10,'k-')
plot(-a.*[-cosd(dip),cosd(dip)],a.*[sind(dip),-sind(dip)],'ko-','LineWidth',2)
axis tight equal
colorbar
clim([-1 1]*0.4)
title('u_x (f_n)')
subplot(2,2,2)
toplot = reshape(ux(:,2),ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
quiver(xg(1:nskip:end)',yg(1:nskip:end)',ux(1:nskip:end,2),uy(1:nskip:end,2),'k')
contour(xvec,yvec,toplot,10,'k-')
plot(-a.*[-cosd(dip),cosd(dip)],a.*[sind(dip),-sind(dip)],'ko-','LineWidth',2)
axis tight equal
colorbar
clim([-1 1]*0.4)
title('u_x (f_s)')
subplot(2,2,3)
toplot=reshape(uy(:,1),ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
quiver(xg(1:nskip:end)',yg(1:nskip:end)',ux(1:nskip:end,1),uy(1:nskip:end,1),'k')
plot(-a.*[-cosd(dip),cosd(dip)],a.*[sind(dip),-sind(dip)],'ko-','LineWidth',2)
axis tight equal
colorbar
clim([-1 1]*0.4)
title('u_y (f_n)')
subplot(2,2,4)
toplot=reshape(uy(:,2),ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
plot(-a.*[-cosd(dip),cosd(dip)],a.*[sind(dip),-sind(dip)],'ko-','LineWidth',2)
axis tight equal
colorbar
clim([-1 1]*0.4)
title('u_y (f_s)')
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

figure(3),clf
subplot(3,1,1)
toplot=reshape(u,ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
plot(-a.*[-cosd(dip),cosd(dip)],a.*[sind(dip),-sind(dip)],'ko-','LineWidth',2)
axis tight equal
colorbar
clim([-1 1]*0.4)
title('u (Laplace)')

subplot(3,1,2)
toplot=reshape(ex,ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
plot(-a.*[-cosd(dip),cosd(dip)],a.*[sind(dip),-sind(dip)],'ko-','LineWidth',2)
axis tight equal
colorbar
clim([-1 1]*0.4)
title('\sigma_x (Laplace)')

subplot(3,1,3)
toplot=reshape(ey,ny,nx);
pcolor(xvec,yvec,toplot), shading interp, hold on
contour(xvec,yvec,toplot,10,'k-')
plot(-a.*[-cosd(dip),cosd(dip)],a.*[sind(dip),-sind(dip)],'ko-','LineWidth',2)
axis tight equal
colorbar
clim([-1 1]*0.4)
title('\sigma_y (Laplace)')

% plot line plots
xpos = 0;
[ux,uy,~,~,~] = calc_disp_stress_forcestokes(xvec(:).*0+xpos,yvec(:),a,dip,eta);
figure(10),clf
subplot(2,1,1)
plot(yvec,ux,'LineWidth',2)
axis tight, grid on
subplot(2,1,2)
plot(yvec,uy,'LineWidth',2)
axis tight, grid on
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 15,'Linewidth',1.5)