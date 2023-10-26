% Two-Asset Rainbow Option Pricing using Backward Difference Method

% Initialising required Parameters
S1 = 95;    % Initial price of asset 1
S2 = 105;   % Initial price of asset 2
K = 100;            % Strike price
r = 0.05;           % Risk-free interest rate
sigma1 = 0.2;       % Volatility of asset 1
sigma2 = 0.2;       % Volatility of asset 2
rho = 0.1;          % Correlation between asset prices
T = 1;              % Time to expiration
N1 = 100;           % Number of asset 1 steps
N2 = 100;           % Number of asset 2 steps
Nt = 100;           % Number of time steps
dt = T/Nt;          % Time step size

% Initialization
dS1 = S1 / N1;
dS2 = S2 / N2;
S1_vec = linspace(0, 2*S1, N1+1)';
S2_vec = linspace(0, 2*S2, N2+1)';

V = zeros(N1+1, N2+1);  % option value matrix

% New boundary conditions
V(:, end) = max(min(S1_vec, S2_vec)-K, 0);  % S2 = 2*S2
V(end, :) = max(min(S1_vec, S2_vec)-K, 0);  % S1 = 2*S1

% Compute coefficient matrices
a1 = 0.5 * dt * (sigma1 * S1_vec / dS1).^2;
a2 = 0.5 * dt * (sigma2 * S2_vec / dS2).^2;
b1 = 0.5 * dt * rho * sigma1 .* sigma2 .* S1_vec .* S2_vec / (dS1 * dS2);
b2 = b1;
c = 1 - (a1 + a2 + b1 + b2 + r * dt);

% Create the tridiagonal matrix A1 using sparse matrix to save space
diag1 = [a1(3:N1); 0] ;
diag2 = c(2:N1);
diag3 = b1(2:N1);
A1 = spdiags([diag3 diag2 diag1], [-1 0 1], N1-1, N1-1);  

% Create the tridiagonal matrix A2 using sparse matrix to save space
diag1 = [a2(3:N2); 0] ;
diag2 = c(2:N2);
diag3 = b2(2:N2);
A2 = spdiags([diag3 diag2 diag1], [-1 0 1], N2-1, N2-1);


for t = Nt:-1:1
    % Compute the option value at each time step
    V_intermediate = A1 * V(2:N1, 2:N2) * A2';
    % Apply the max-min style payoff
    V_intermediate = max(min(S1_vec(2:N1), S2_vec(2:N2)) - K, V_intermediate);
    % Update the option value matrix
    V(2:N1, 2:N2) = V_intermediate;
end

% optionValue = V(N1+1, N2+1)
backward_dif_option_value = interp2(S2_vec, S1_vec, V, S2, S1)


