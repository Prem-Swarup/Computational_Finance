% This code uses the binomial tree method to price rainbow-style two-asset options by 
% modeling stock price movements, calculating option values at each time step, and 
% applying backward induction to determine the final option price.


% Set initial parameters
T = 1;              % Time to expiration
N = 5;             % Number of time steps
dt = T/N;           % Time increment
r = 0.05;           % Risk-free interest rate
K = 100;            % Strike price
S1_0 = 95;           % Initial price of stock 1
S2_0 = 105;          % Initial price of stock 2
rho = 0.1;          % Correlation between stock prices
sigma1 = 0.2;       % Volatility of stock 1
sigma2 = 0.2;       % Volatility of stock 2

% Calculate up and down factors
u1 = exp(sigma1*sqrt(dt));
u2 = exp(sigma2*sqrt(dt));
d1 = 1/u1;
d2 = 1/u2;

% Calculate stock prices at each time step
i = (1:N+1).';
S1_t = S1_0 * u1.^(i) .* d1.^(N-i);
S2_t = S2_0 * u2.^(i) .* d2.^(N-i);

% Initialize option value matrix
V = zeros(N+1, N+1, N+1);

% Set final option values at expiration
V(:,:,N+1) = max(min(S1_t, S2_t.') - K, 0);

% Calculate option values at earlier time steps using backward induction
p_uu = 0.25*(1 + rho); %probability that both assets go up
p_ud = 0.25*(1 - rho); %probability that 1st asset go up and other down
p_du = 0.25*(1 - rho); %probability that 2nd asset go up and other down
p_dd = 0.25*(1 + rho); %probability that both assets go up

% Applying backward induction to determine option price
for n = N-1:-1:0
    for i = 1:n+1
        for j = 1:n+1
            V(i,j,n+1) = exp(-r*dt)*(p_uu*V(i+1,j+1,n+2) + p_ud*V(i+1,j,n+2) + p_du*V(i,j+1,n+2) + p_dd*V(i,j,n+2));
        end
    end
end

% The final option value is the value at time 0
rainbow_call_price = V(1,1,1)

