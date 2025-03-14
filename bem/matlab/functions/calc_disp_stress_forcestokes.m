function [ux,uy,sxy,sxx,syy] = calc_disp_stress_forcestokes(X,Y,a,dip,eta)
% the elementary Stokeslet solutions are for a line centered at (0,0) 
% extending from [-a <= y <= a]
% X,Y - are actually [X-Xs,Y-Ys]
% forces applied are in x-normal,(-y)-shear direction - relative to a
% vertical line
% convert from normal,shear coordinates to local x-y coordinates
% 
% AUHTORS:
% Rishav Mallick, JPL, Caltech, 2025


beta = dip-90;

% first rotate from [X,Y] to [x,y] by dip angle
R = [cosd(beta),-sind(beta);...
     sind(beta), cosd(beta)];
rot_coords = [X,Y]*R;
xo = rot_coords(:,1);
yo = rot_coords(:,2);
% set source x position to 0
xs = 0;

uxfun = @(fx,fy) (1/4).*pi.^(-1).*eta.^(-1).*(4.*a.*fx+((-1).*a.*fx+(-1).*fy.*xo+fy.* ...
    xs+fx.*yo).*log((xo+(-1).*xs).^2+(a+(-1).*yo).^2)+(-1).*(fy.*((-1) ...
    .*xo+xs)+fx.*(a+yo)).*log((xo+(-1).*xs).^2+(a+yo).^2));

% remedy for branchcuts
atanfun = -(atan((a-yo)./(xo-xs)) + atan((a+yo)./(xo-xs)));

% ux,uy functions
uyfun = @(fx,fy) (1/4)/pi/eta.*(8.*a.*fy +...
    4.*fy.*(xo-xs).*(atanfun)+...
    ((-1).*a.*fy+(-1).*fx.*xo+fx.*xs+fy.*yo).* ...
    log((xo+(-1).*xs).^2+(a+(-1).*yo).^2)+(-1).*(fx.*((-1).*xo+xs)+ ...
    fy.*(a+yo)).*log((xo+(-1).*xs).^2+(a+yo).^2));

% stress component functions
sxyfun = @(fx,fy) (1/2).*pi.^(-1).*(xo+(-1).*xs).*(2.*a.*((xo+(-1).*xs).^2+(a+(-1).* ...
    yo).^2).^(-1).*(fy.*(a.^2+(xo+(-1).*xs).^2)+2.*fx.*((-1).*xo+xs).* ...
    yo+(-1).*fy.*yo.^2).*((xo+(-1).*xs).^2+(a+yo).^2).^(-1)+(-1).*fy.* ...
    (xo+(-1).*xs).^(-1).*(atan((xo+(-1).*xs).^(-1).*(a+(-1).*yo))+ ...
    atan((xo+(-1).*xs).^(-1).*(a+yo))));

sxxfun = @(fx,fy) (1/4).*pi.^(-1).*((-4).*a.*(xo+(-1).*xs).*((xo+(-1).*xs).^2+(a+( ...
    -1).*yo).^2).^(-1).*(fx.*(a.^2+(xo+(-1).*xs).^2)+2.*fy.*(xo+(-1).* ...
    xs).*yo+(-1).*fx.*yo.^2).*((xo+(-1).*xs).^2+(a+yo).^2).^(-1)+(-1) ...
    .*fy.*log((xo+(-1).*xs).^2+(a+(-1).*yo).^2)+fy.*log((xo+(-1).*xs) ...
    .^2+(a+yo).^2));

syyfun = @(fx,fy) (1/4).*pi.^(-1).*(4.*a.*(xo+(-1).*xs).*((xo+(-1).*xs).^2+(a+(-1).* ...
    yo).^2).^(-1).*(fx.*(a.^2+(xo+(-1).*xs).^2)+2.*fy.*(xo+(-1).*xs).* ...
    yo+(-1).*fx.*yo.^2).*((xo+(-1).*xs).^2+(a+yo).^2).^(-1)+fy.*log(( ...
    xo+(-1).*xs).^2+(a+(-1).*yo).^2)+(-1).*fy.*log((xo+(-1).*xs).^2+( ...
    a+yo).^2));

% compute displacements in rotated coordinates
ux_r = [uxfun(1,0),uxfun(0,-1)];
uy_r = [uyfun(1,0),uyfun(0,-1)];
% compute stress in rotated coordinates
sxy_r = [sxyfun(1,0),sxyfun(0,-1)];
sxx_r = [sxxfun(1,0),sxxfun(0,-1)];
syy_r = [syyfun(1,0),syyfun(0,-1)];

% rotate back to original coordinates
U1 = [ux_r(:,1),uy_r(:,1)]*R';
U2 = [ux_r(:,2),uy_r(:,2)]*R';
ux = [U1(:,1),U2(:,1)];
uy = [U1(:,2),U2(:,2)];
 
S1_rot = R(1,1)^2 * sxx_r(:,1) + 2 * R(1,1) * R(1,2) * sxx_r(:,1) + R(1,2)^2 * syy_r(:,1);
S2_rot = R(1,1) * R(2,1) * sxx_r(:,1) + (R(1,1) * R(2,2) + R(1,2) * R(2,1)) * sxy_r(:,1) + R(1,2) * R(2,2) * syy_r(:,1);
S3_rot = R(2,1)^2 * sxx_r(:,1) + 2 * R(2,1) * R(2,2) * sxy_r(:,1) + R(2,2)^2 * syy_r(:,1);
sxy = [S2_rot,sxy_r(:,2)];
sxx = [S1_rot,sxx_r(:,2)];
syy = [S3_rot,syy_r(:,2)];

end