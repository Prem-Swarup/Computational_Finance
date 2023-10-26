% This code uses the binomial tree method to price digital-style two-asset options by 
% modeling stock price movements, calculating option values at each time step, and 
% applying backward induction to determine the final option price.

% initial parameters as given
T = 1;                      % Time to expiration in years
N = 70;                    % Number of time steps
dt = T/N;                   % Length of each time step in years
r = 0.05;                   % Risk-free interest rate
K = 100;                    % Strike price
S1_0 = 95;                  % Initial price of stock 1
S2_0 = 105;                 % Initial price of stock 2
rho = 0.1;                  % Correlation between stock prices
sigma1 = 0.2;               % Volatility of stock 1
sigma2 = 0.2;               % Volatility of stock 2

u1 = exp(sigma1*sqrt(dt));  % Up factor for stock 1
u2 = exp(sigma2*sqrt(dt));  % Up factor for stock 2
d1 = 1/u1;                  % Down factor for stock 1
d2 = 1/u2;                  % Down factor for stock 2

% Calculating stock prices at each time step
i = (0:N).';
S1_v = S10 * u1.^(i-1) .* d1.^(N-i+1);  % stock price of asset 1 at various levels of up and down
S2_v = S20 * u2.^(i-1) .* d2.^(N-i+1);  % stock price of asset 2 at various levels of up and down

% Initialize option value matrix
V = zeros(N+1, N+1, N+1);

% Set final option values at expiration
for i = 1:N+1
    for j = 1:N+1
        if S1_v(i) <= K && S2_v(j) <= K
            V(i,j,N+1) = 1;  % Option value is 1 if both stock prices are below the strike price
        else
            V(i,j,N+1) = 0;  % Option value is 0 otherwise
        end
    end
end

% Calculate option values at earlier time steps using backward induction
puu = 0.25*(1 + rho);  % Probability of up-up movement
pud = 0.25*(1 - rho);  % Probability of up-down movement
pdu = 0.25*(1 - rho);  % Probability of down-up movement
pdd = 0.25*(1 + rho);  % Probability of down-down movement
for n = N-1:-1:0
    for i = 1:n+1
        for j = 1:n+1
            V(i,j,n+1) = exp(-r*dt)*(puu*V(i+1,j+1,n+2) + pud*V(i+1,j,n+2) + pdu*V(i,j+1,n+2) + pdd*V(i,j,n+2));
            % Calculate option value at time step n using discounted expected value of future options
        end
    end
end

% The final option value is the value at time 0
digital_option_price = V(1,1,1)


