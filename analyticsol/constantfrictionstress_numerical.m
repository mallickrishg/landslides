clear
% Define parameters
L = 1;    % Domain length
a = 1;    % Example parameter 'a'
b = 0.5;    % Example parameter 'b'

% Initial mesh
x = linspace(0, L-0.01, 1e3);

% functions needed to solve ODE
odefun = @(x, y) odefun_param(x, y, a, b);
bcfun = @(ya, yb) bcfun_param(ya, yb);
guess = @(x) guess_param(x,L);

% Initial guess for solution
solinit = bvpinit(x, guess);

% Solve using bvp4c
options = bvpset('RelTol', 1e-7, 'AbsTol', 1e-9,'Nmax',1e5); % Higher accuracy
sol = bvp5c(odefun, bcfun, solinit,options);

%% Plot solution
figure(1), clf
subplot(2,1,1)
plot(sol.x,sol.y(1, :), '-', 'LineWidth', 3), hold on
dummy = guess(sol.x);
plot(sol.x,dummy(1,:),'k-','LineWidth',2)
xlabel('x');
ylabel('y(x)');
grid on, axis tight
% ylim([0 1.5])
set(gca,'FontSize',15,'LineWidth',1.5)

subplot(2,1,2)
plot(sol.x, sol.y(2, :), '-', 'LineWidth', 3);
xlabel('x');
ylabel('dy/dx');
grid on;
ylim([-1 1]*9)
set(gca,'FontSize',15,'LineWidth',1.5)


%% write out all the ODE functions
% Nested function: ODE system
function dydx = odefun_param(x,y,a,b)
% y1 -> h(x)
% y2 -> dh/dx
delta = 1e-4; % adding regularization
dydx = [y(2);...
       (1 - y(2)^2 - a*(b - y(2))) / (y(1) + delta)];
end

% Nested function: Boundary conditions
function res = bcfun_param(ya, yb)
% res = [ya(1); yb(1) - 1/2*(1-sqrt(5))]; % Enforcing y(-L) = 0 and y'(L) = 0
res = [ya(1)-1;...
       yb(1)-0];% Enforcing y(-L) = 0 and y(L) = 0
end

% Nested function: Initial guess
function y_guess = guess_param(x,L)
y_guess = [1 - (x./L).^2; -2*(x./L)]; % A simple quadratic guess
end

