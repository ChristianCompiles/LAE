function dUdt = rhs_lae_u(U, V, Pinv, D, E_l, E_r, a, b, sigma, g, t)
% 
%  Implement right-hand-side computation for LAE
%  on left subdomain using an SBP-SAT scheme.
% 
%    INPUT: U - solution vector
%           V - solution vector
%           Pinv - inverse of the integration matrix P
%           D - SBP differentiation matrix
%           E_l - zeros matrix with mat(1,1) = 1.
%           E_r - zeros matrix with mat(N,N) = 1.
%           a - wave speed in left subdomain.
%           b - wave speed in right subdomain (needed to calculate c)
%           sigma - penalty term of the SAT
%           g(t) - (possibly) time dependent boundary data function
%           t - current time
% 
%    OUTPUT: dUdt - RHS of the SBP-SAT linear advection scheme on left
%    subdomain

   % declare memory for output variable
   N = length(U);
   dUdt = zeros(N, 1);
   c = a/b;

   % SAT to penalize U solution against V solution at interface boundary.
   SAT_r = sigma * Pinv * E_r * (U*c - [zeros(N-1,1); V(1)]);
   
   sigmaX = -1; % <= -0.5
   
   % SAT to penalize U solution against initial condition
   % at start of left subdomain.
   SAT_l = sigmaX * a * Pinv * E_l * (U - [g(t); zeros(N-1,1)]);

   dUdt = -a * D * U + SAT_r + SAT_l;

end