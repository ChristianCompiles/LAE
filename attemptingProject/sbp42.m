function [P, D] = sbp42(N, dx)
%
%  Create a 4-2 SBP operator that is fourth order
%  in the interior with a second order boundary closure
%  on the interval [a, b].
%
%  This creates full matrices for ease of use but
%  could use sparse matrices if speed is desired
%  (see commented code at the end of the file).
%
%  Coefficients come from the paper where, note, the integration
%  matrix P is referred to as H:
%
%  Mattsson, Nordström (2004)
%   "Summation by parts operators for finite difference 
%    approximations of second derivatives."
%   Journal of Computational Physics 199, pp. 503-540.
%  
%  INPUT:  N  - number of uniform discretization points
%          dx - grid spacing defined by (b-a) / (N-1) 
%
%  OUTPUT: P - integration matrix
%          D - derivative matrix
%

    % Check that it is possible to create the operator
    if N < 8
       error("Must use N >= 8 to construct the 4-2 SBP operator.")
    end

    % For convenience create P matrix and its inverse
    P = dx * diag([17/48, 59/48, 43/48, 49/48, ...
                   ones(1, N-8), ...
                   49/48, 43/48, 59/48, 17/48]);
    % Pinv = diag(1 ./ diag(P)); % unused
    
    % Standard fourth order central stencil for the interior
    D = -1/12 * diag(ones(N-2,1), 2) ...
        +2/3  * diag(ones(N-1,1), 1) ...
        -2/3  * diag(ones(N-1,1),-1) ...
        +1/12 * diag(ones(N-2,1),-2);

    % Create the somewhat ugly boundary closures on the left
    D(1:4, 1:6) = [-24/17, 59/34, -4/17, -3/34, 0, 0;
                   -1/2, 0, 1/2, 0, 0, 0;
		              4/43, -59/86, 0, 59/86, -4/43, 0; 
		              3/98, 0, -59/98, 0, 32/49, -4/49];

    % Reverse the order and flip the sign of the boundary closure on the right
    D(N-3:N,N-5:N) = flipud(fliplr(-D(1:4,1:6)));

    % Scale by the grid spacing
    D = D / dx;

    % Q = P*D; % unused to create the unweighted differencing matrix
end

%%
% original sparse implementation saved for reference.
    % For convenience create P matrix and its inverse
    % P = dx * diag([17/48, 59/48, 43/48, 49/48, ...
    %                ones(1, N-8), ...
    %                49/48, 43/48, 59/48, 17/48]);
    % % Pinv = diag(1 ./ diag(P)); % unused
    % P = sparse(P);
    % Standard fourth order central stencil for the interior
    % D = -1/12 * diag(ones(N-2,1), 2) ...
    %     +2/3  * diag(ones(N-1,1), 1) ...
    %     -2/3  * diag(ones(N-1,1),-1) ...
    %     +1/12 * diag(ones(N-2,1),-2);
    % 
    % % Create the somewhat ugly boundary closures on the left
    % D(1:4, 1:6) = [-24/17, 59/34, -4/17, -3/34, 0, 0;
    %                -1/2, 0, 1/2, 0, 0, 0;
		%               4/43, -59/86, 0, 59/86, -4/43, 0; 
		%               3/98, 0, -59/98, 0, 32/49, -4/49];
    % 
    % % Reverse the order and flip the sign of the boundary closure on the right
    % D(N-3:N,N-5:N) = flipud(fliplr(-D(1:4,1:6)));
    % 
    % % Scale by the grid spacing
    % D = sparse(D / dx);
