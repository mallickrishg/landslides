% Construct meshes for landslide problem
% and store as files that can be read by geometry.receiver()
% - topography is given by steady-state solution through local stress
% balance assuming constant yield stress and perturbed by a sinusoid
% - basal slip surface can have sinusoidal topography as well
% - labels added for each boundary are as follows:
% 0 - velocity BC OR dp/dn = 0
% 1 - traction BC OR p = 0
% 2 - mixed BC OR dp/dn is prescribed as non-zero
% 
% NOTE: the normal vector to the boundary points outward so make sure that
% all potential field and derivatives are modified accordingly.
% 
% AUTHORS:
% Rishav Mallick, JPL, Caltech 
% March, 2025

clear

addpath functions/
import('geometry.*')

% assume landlside properties
rhog = 1; % landslide density - ρg
theta = 30; % slope of landslide in degrees

% construct landslide geometry outline
npts = 50;
Lx = 5;
xpts = linspace(0,Lx,npts+1)';
% number of points needed to resolve left boundary
npts_left = npts/Lx;

% lambert-w function topography
topovals = cscd(theta) * (1 + lambertw(-exp(-1 + (xpts/Lx-1) * rhog * sind(theta) * tand(theta)))) / rhog + ...
            0.00*sin(8*pi*xpts/Lx);
ztopo = topovals/max(topovals);
zleft = linspace(0,1,npts_left+1)';
% sinusoidal topography for bottom slip surface
zbottom = 0.0*sin(1*pi*xpts/Lx);

% create outline mesh
xg = [xpts(1:end-1);...
     flipud(xpts(2:end));...
     zeros(npts_left,1)];
zg = [ztopo(1:end-1);...
     flipud(zbottom(2:end));...
     zleft(1:end-1)];
xgp = [xg;0];
zgp = [zg;1];
W = sqrt((xgp(1:end-1)-xgp(2:end)).^2 + (zgp(1:end-1)-zgp(2:end)).^2);
dip = atan2d(-zgp(1:end-1)+zgp(2:end),-xgp(1:end-1)+xgp(2:end));
% add labels for BC type 
% 0 - velocity BC
% 1 - traction BC
% 2 - mixed BC
labels = [ones(npts,1);2.*ones(npts,1);zeros(npts_left,1)];

fparams = [labels,xg,zg,W,dip,W];
eM = geometry.LDhs(1,0.25);% dummy parameters
rcv = geometry.receiver(fparams,eM);

% use outline to construct interior mesh
model = createpde;

R = [2,length(xg),xg',zg']';
g = decsg(R);
geometryFromEdges(model,g);
meshmodel = generateMesh(model,"Hmax",0.2,"Hgrad",1.1,"GeometricOrder","linear");

figure(1),clf
pdeplot(meshmodel), hold on
quiver(rcv.xc(:,1),rcv.xc(:,2),rcv.nv(:,1),rcv.nv(:,2),0.1,'k')
quiver(rcv.xc(:,1),rcv.xc(:,2),rcv.dv(:,1),rcv.dv(:,2),0.1,'r')
plotpatch2d(rcv)
axis tight equal

%% save mesh files in appropriate folder
p = meshmodel.Nodes';
t = meshmodel.Elements';

writetable(table(p),'mesh/geometry_vertices.dat','WriteVariableNames',false)
writetable(table(t),'mesh/geometry_triangulation.dat','WriteVariableNames',false)

save('mesh/rcv.mat','rcv')