function [Disp,Strain] = compute_disp_stress_kernels_force(rcv,obs)
% compute displacement & displacement gradient kernels for a discretized
% line mesh using uniform force greens function solutions
% AUTHOR
% Rishav Mallick, JPL, 2023

x = obs(:,1);
y = obs(:,2);

Disp = zeros(length(x),rcv.N);
Strain = zeros(length(x),rcv.N,2);
x2c = rcv.xc(:,1);
x3c = rcv.xc(:,2);

% shift by infinitesimal amount in the normal direction
dr = -1e-9;

for i = 1:rcv.N
    a = rcv.W(i)/2;
    [u,~,~] = calc_disp_strain_forcefault(x-x2c(i),y-x3c(i),a,rcv.dip(i));
    
    % shift coincident strain kernel by epsilon
    r = sqrt((x-x2c(i)).^2 + (y-x3c(i)).^2);
    index = r == 0;
    xmod = x-x2c(i);
    ymod = y-x3c(i);
    xmod(index) = xmod(index) + dr*rcv.nv(i,1);
    ymod(index) = ymod(index) + dr*rcv.nv(i,2);
    [~,dudx,dudz] = calc_disp_strain_forcefault(xmod,ymod,a,rcv.dip(i));

    Disp(:,i) = u;
    Strain(:,i,1) = dudx;
    Strain(:,i,2) = dudz;
end

end