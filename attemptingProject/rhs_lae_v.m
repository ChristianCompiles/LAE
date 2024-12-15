function dVdt = rhs_lae_v(U, V, Pinv, D, E_l, a, b, sigma)
% 
%  Implement right-hand-side computation for LAE 
%  on right subdomain using an SBP-SAT scheme.
% 
%    INPUT: U - solution vector
%           V - solution vector
%           Pinv - inverse of the integration matrix P
%           D - SBP differentiation matrix
%           E_l - zeros matrix with mat(1,1) = 1.
%           a - wave speed in left subdomain (needed to calculate c)
%           b - wave speed in right subdomain
%           sigma - penalty term of the SAT
%           t - current time
% 
%    OUTPUT: dVdt - RHS of the SBP-SAT linear advection scheme on right
%    subdomain.

   % declare memory for output variable
   N = length(V);
   dVdt = zeros(N, 1);
   c = a/b;

   % Penalize V solution against U solution at interface.
   SAT = sigma * Pinv * E_l * (V - [c*U(end); zeros(N-1,1)]);

   dVdt = -b * D * V + SAT;

end