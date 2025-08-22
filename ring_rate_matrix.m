function Q = ring_rate_matrix(L,f)

% Inputs. The function takes as input an integer colony size N ≥ 2, as well as
% switching and recruitment parameters eps and r, which are positive real
% numbers.
%
% Outputs. It returns an (N+1)-by-(N+1) matrix Q, which is the generator,
% or transition rate matrix, of the Follmer--Kirman model.
%
% Author: Jacob Calvert (calvert@gatech.edu)
% Date: June 4, 2024

% The size of rate matrix Q needs to be N+1 because you can have 0 ants on
% the first path.
Q = zeros(N+1,N+1);

% Calculate the forward rates (there are only nearest-neighbor jumps)
for x=1:N
    Q(x,x+1) = eps*(N-(x-1)) + r*(x-1)*(N-(x-1))/(N-1);
end

% Calculate the reverse rates (again, nearest-neighbor jumps only)
for y=2:N+1
    Q(y,y-1) = eps*(y-1) + r*(y-1)*(N-(y-1))/(N-1);
end

% The diagonal elements of Q are set to make the row sums equal to 0
Q = Q - diag(diag(Q));
Q = Q - diag(sum(Q,2));

end