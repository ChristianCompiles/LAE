

% setup the SBP 2-1 operator on an interval [x_l, x_r]

% setup the domain and resolution

x_r = 1.0;
x_l = 0.0;

% resolution and get x points
N = 10;

% create uniform gridpoints
dx = (x_r - x_l) / (N - 1);
x = transpose(x_l:dx:x_r);

% create the P and D matrices (also construct Q for convenience)

% integration matrix
P = diag(dx * [0.5, ones(1, N-2), 0.5], 0);

% derivative matrix
D = diag(0.5 * ones(1, N-1), 1) + diag(-0.5 * ones(1, N-1), -1);

% modify the first and last rows using the biased stencil
D(1,:) = [-1, 1, zeros(1,N-2)];
D(N,:) = [zeros(1,N-2), -1, 1];

% do not forget to divide by delta x
D = D / dx;

% compute Q directly
Q = P*D;

%%
% check the SBP condition as Q + tranpose(Q) = B = diag(-1,0,...,0,1)
B = Q + transpose(Q)

% check the accuracy conditions by differentiating polynomials.
% The 2-1 SBP operator should differentiate up to linear polynomials exactly

D * ones(N, 1) % will give a vector of zeros

D * x % will give a vector of ones

D * x.^2 % will be exact in the interior but inexact at the boundary




