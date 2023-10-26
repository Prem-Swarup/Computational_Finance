% Monte Carlo Simulation for Two-Asset Option Pricing
% It generates multiple asset price by taking account in volatality and
% correlation, whose mean should converge to actual value of option

% Initialising required parameters
S10 = 95; % Initial stock price of asset 1
S20 = 105; % Initial stock price of asset 2
K = 100; % Strike price of the option
r = 0.05; % Risk-free interest rate
T = 1; % Time to maturity
sigma1 = 0.2; % Volatility of asset 1
sigma2 = 0.2; % Volatility of asset 2
rho = 0.1; % Correlation between the two assets
d = 2; % Number of assets
N = 100000; % Number of simulations

% Cholesky decomposition
sigma_sq = [sigma1^2 rho*sigma1*sigma2; rho*sigma1*sigma2 sigma2^2];
L = chol(sigma_sq, 'lower');

% Simulating stock prices by taking account in volatality by generating random numbers normally
Z = randn(N, d);
S = [S10*exp((r-0.5*sigma1^2)*T + sqrt(T)*Z(:,1)), S20*exp((r-0.5*sigma2^2)*T + sqrt(T)*(L*Z(:,1:2)')')];

% Computing option payoffs
payoff = max(min(S(:,1), S(:,2)) - K, 0);

% Computing option price
option_price = exp(-r*T)*mean(payoff)
