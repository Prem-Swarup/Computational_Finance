% Monte Carlo Method for Pricing Two-Asset Digital Option

% Parameters
S1_0 = 95; % Initial stock price of asset 1
S2_0 = 105; % Initial stock price of asset 2
K = 100; % Strike price of the option
r = 0.05; % Risk-free interest rate
T = 1; % Time to maturity
sigma1 = 0.2; % Volatility of asset 1
sigma2 = 0.2; % Volatility of asset 2
rho = 0.1; % Correlation between the two assets
num_simulations = 1000; % Number of simulations



% Cholesky decomposition
sigma_sq = [sigma1^2 rho*sigma1*sigma2; rho*sigma1*sigma2 sigma2^2];
L = chol(sigma_sq, 'lower');


% Simulating stock prices by taking account in volatality by generating random numbers normally
Z = randn(num_simulations, 2);
S = [S1_0*exp((r-0.5*sigma1^2)*T + sqrt(T)*Z(:,1)), S2_0*exp((r-0.5*sigma2^2)*T + sqrt(T)*(L*Z(:,1:2)')')];


% Computing option payoffs 
payoff = (S(:,1) <= K & S(:,2) <= K);

% Computing option price
monte_carlo_option_price = exp(-r * T) * mean(payoff)

