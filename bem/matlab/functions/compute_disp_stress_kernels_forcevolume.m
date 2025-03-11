function [Disp,Strain] = compute_disp_stress_kernels_forcevolume(mesh,obs)

N = length(obs(:,1));
x = obs(:,1);
z = obs(:,2);

% kernels
Disp = zeros(N,mesh.N);
Gux = zeros(N,mesh.N);
Guz = zeros(N,mesh.N);
parfor iter=1:mesh.N
    % compute quadrature points and weights
    Q = quadtriangle(25,'domain',mesh.vert(mesh.tri(iter,:),:),'Type','nonproduct');
    % integrate numerically
    Disp(:,iter) = pointforce_u(x,z,Q.Points(:,1)',Q.Points(:,2)')*Q.Weights;
    Gux(:,iter) = pointforce_dudx(x,z,Q.Points(:,1)',Q.Points(:,2)')*Q.Weights;
    Guz(:,iter) = pointforce_dudz(x,z,Q.Points(:,1)',Q.Points(:,2)')*Q.Weights;
end

% store displacement gradient kernels
Strain = zeros(N,mesh.N,2);
Strain(:,:,1) = Gux;
Strain(:,:,2) = Guz;

end

%% define point source functions to be integrated numerically
function u = pointforce_u(x,y,xs,ys)
r = sqrt((x-xs).^2 + (y-ys).^2);
u = 1/2/pi.*log(r);
end
function dudx = pointforce_dudx(x,y,xs,ys)
r = sqrt((x-xs).^2 + (y-ys).^2);
dudx = 1/2/pi.*(x-xs)./(r.^2);
end
function dudz = pointforce_dudz(x,y,xs,ys)
r = sqrt((x-xs).^2 + (y-ys).^2);
dudz = 1/2/pi.*(y-ys)./(r.^2);
end