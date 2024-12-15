function [Un, Vn] = step_by_rk3(t, dt, Un, Vn, Pinv, D, E_l, E_r, a, b, sigmaL, sigmaR, bndy_data_u)
%
% Un = stepByRK3(dt, Un, Pinv, D, E, a, sigma, bndy_data)
% 
%  Third order in time, three stage low-storage Runge-Kutta method from
%    Williamson (1980) "Low-storage Runge-Kutta schemes" 
%    Journal of Computational Physics
%    https://doi.org/10.1016/0021-9991(80)90033-9
% 
%    INPUT: dt - explicit time step size
%           Un - solution vector at the current time
%           Pinv - inverse of the integration matrix P
%           D - SBP differentiation matrix
%           E - convenience matrix for the SAT, e.g., E = diag(1, 0, ... ,0)
%           a - wave speed of the linear advection equation
%           sigma - penalty term of the SAT which must be <= -0.5
%           bndy_data(t) - (possibly) time dependent boundary data function
% 
%    OUTPUT: Un - solution vector updated one explicit time step dt

   % coefficients for the low-storage Runge-Kutta 3 method
   % due to Williamson
   coeffs_a = [0;   -5/9; -153/128];
   coeffs_b = [0;    1/3;  3/4];
   coeffs_g = [1/3; 15/16; 8/15];
   
   % Declare memory
   dUdt = zeros(size(Un));
   dVdt = zeros(size(Un));
   G = zeros(size(Un));
   H = zeros(size(Un));
   for k = 1:3
      % Compute the local time for the current Runge-Kutta stage
      t_loc = t + coeffs_b(k) * dt;
      % Assemble the right-hand-operator
      dUdt = rhs_lae_u(Un, Vn, Pinv, D, E_l, E_r, a, b, sigmaL, bndy_data_u, t_loc);
      dVdt = rhs_lae_v(Un, Vn, Pinv, D, E_l, a, b, sigmaR);
      % G,H are a running combination of previous and current Runge-Kutta stages
      G = coeffs_a(k) * G + dUdt;
      H = coeffs_a(k) * H + dVdt;
      % Update the solution
      Un = Un + coeffs_g(k) * dt * G;
      Vn = Vn + coeffs_g(k) * dt * H;
   end

end