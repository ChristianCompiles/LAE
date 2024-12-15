
% driver script for approximating the linear advection equation
% solution via an SBP-SAT scheme on a domain with interface boundary.
%
%    u_t + a u_x = 0, a is a positive real constant on [x_l,x_r]
%    u(x,0) = u_0(x) % initial condition
%    u(a,t) = g(t)   % boundary data

close all % to close old plot windows

% set the domain, make it symmetric, assuming boundary interface is at 0
x_l = -2.0;
x_r = 2.0;

% resolution and get x points
N = 150; % size of both left and right subdomain

% % create uniform gridpoints

% Interface at x=0
dx_u = (0.0 - x_l) / (N - 1);
dx_v = (x_r - 0.0) / (N - 1);

x_u = transpose(x_l:dx_u:0.0);
x_v = transpose(0.0:dx_v:x_r);

% create the SBP operator pair
% [P, D] = sbp42(N, dx_u);
%[P, D] = sbp63(N, dx_u); % for the time being assume that N is the same on both sides so that this is "probably" kosher
[P, D] = sbp84(N, dx_u);

% store P inverse for convenience
Pinv = zeros(N,N);
for j = 1:N
   Pinv(j,j) = 1 / P(j,j);
end

% boundary matrix to build the SAT on the left with E=diag(1,0,...,0,0)
% OBS! if you switch a to be negative then the boundary data is set at the
%      right and this matrix would be E=diag(0,0,...,0,1)

E_r = zeros(N,N); % e_N in Cognata
E_r(N,N) = 1;

E_l = zeros(N,N); % e_0 in Cognata
E_l(1,1) = 1;

a = 1.5; % Left domain wavespeed.
b = 1; % Right domain wavespeed.

% Set the SAT penalty parameter. For LAE must be <= -1/2.
sigma_l = b/2;
sigma_r = sigma_l - b;

% % setup the manufatured solution and boundary term
% u_ex = @(x,t) 2 + sin(2 * pi * (x - a * t));
% g = @(t) 2 + sin(2 * pi * (x_l - a * t));

% setup the manufactured solution
u_ex = @(x,t) sin(4*x-a*pi*(-1+3*t));
%u_ex = @(x,t) exp(-20*((x - a * t) + 1.5).^2);
%g = @(t) 0.0; % same as f(t) = 0; 

% Setup boundary term
g_l = @(t) u_ex(x_l, t); % manufactured solution only on left subdomain.

% % Try something fancier that emphasizes what the boundary condition does
% u_ex = @(x,t) zeros(size(x));
% g = @(t) boundary_condition_sine_sector(t);

% Set the initial condition.
% U = u_ex(x_u, 0.0);
% V = u_ex(x_v, 0.0);
U = zeros(size(x_u));
V = zeros(size(x_v));

% set the time step size
CFL = 0.8; % for SBP 42 or 63 or 84
dt = min(CFL * dx_u / abs(a), CFL * dx_v / abs(b));
t = 0.0;
t_final = 4.0; % 1.0 for error analysis
k = 0; % # of time steps (only used for plotting data)


% Do the time loop
while t < t_final
   % Avoid stepping over the final time because we use a while loop
   if t + dt > t_final
     dt = t_final - t;
   end

   [U,V] = step_by_rk3(t, dt, U, V, Pinv, D, E_l, E_r, a, b, sigma_l, sigma_r, g_l);
   
   k = k+1;
   t = t + dt;
   % plot on the fly to show a "movie"
   % OBS! somewhat ugly way to do this
   if mod(k, 5) == 0
      % plot(x,u_ex(x,t),'-k','LineWidth',1.5)
      plot(x_u, U, '-m', 'LineWidth', 1.5)
      hold on
      plot(x_v, V, '--k', 'LineWidth', 1.5)
   
      plot(x_l, g_l(t), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r')
      hold off
      % xlim([x(1) x(end)])
       ylim([-1.5 1.5])
      %ylim([-0.5 1.5])
      xlabel('$x$', 'interpreter', 'latex')
      % legend('$u^{\textrm{exact}}$', '$u$', 'interpreter', 'latex')
      legend('$U$', '$V$', '$g(t)$', 'interpreter', 'latex', 'location', 'southwest')
      title(['$t = $',num2str(t)], 'interpreter', 'latex')
      set(gca, 'fontsize', 24)

      pause(0.005)
   end
end

% %% Commands to check the operator spectrum
% 
% % draw the stability region for RK3 (or forward Euler)
% [x, y] = meshgrid(-3:0.1:1, -3:0.1:3);
% z = x + 1i * y;
% 
% % % forward Euler stability region
% % FE = 1 + z; 
% % FEModulus = abs(FE);
% % contourf(x, y, -FEModulus, [-1 -1])
% 
% % Runge-Kutta 3 stability region
% RK3 = 1 + z + z.^2/2 + z.^3/6;
% RK3Modulus = abs(RK3);
% contourf(x, y, -RK3Modulus, [-1 -1])
% 
% % setup nice axes and labels
% hold on
% 
% plot([x(1) x(end)], [0 0], '-k')
% plot([0 0], [y(1) y(end)], '-k')
% xlabel('Re(\lambda)')
% ylabel('Im(\lambda)')
% set(gca, 'FontSize',16)
% % axis equal
% 
% % Compute operator spectrum that collects everything that "hits" U
%   operator = -a * D + a * sigma * Pinv * E;

%operator_u = -a * D - a * Pinv * E;

%   lamb = eig(operator);
% % raw spectra
%   % plot(real(lamb), imag(lamb), 'ro', 'MarkerFaceColor', 'r')
% % spectra scaled by dt
%   CFL = 1.0;
%   dt = CFL * dx / a;
%   plot(real(lamb)*dt, imag(lamb)*dt, 'ro', 'MarkerFaceColor', 'r')
% 
% hold off
