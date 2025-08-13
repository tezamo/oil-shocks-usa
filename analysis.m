%% 0) Data
% Load the data
data = readmatrix('data.xlsx', 'Range', 'B2:E577');  

% Define variables
[T, k] = size(data);  % T is the number of observations, k is the number of variables

%% 1) Obtaining optimal Lag Legth using AIC information criteria
maxLags = 20;  % Assume a maximum of 20 lags for this example
aic = zeros(maxLags, 1);
bic = zeros(maxLags, 1);

% Loop through different lag lengths
for p = 1:maxLags
    X = createLaggedMatrix(data, p);
    X = X(p+1:end, :);  % Remove rows with NaN values
    Y = data(p+1:end, :);  % Adjust the range accordingly
    
    beta = (X\Y);  % Estimate the coefficients by OLS
    
    residuals = Y - X*beta;
    Sigma = cov(residuals);
    
    aic(p) = log(det(Sigma)) + (2*k*p)/T;
    bic(p) = log(det(Sigma)) + (k*p*log(T))/T;
end

% Find the optimal lag length
[~, optimalLagAIC] = min(aic);
[~, optimalLagBIC] = min(bic);

% Display the optimal lag length
disp(['Optimal Lag Length (AIC): ', num2str(optimalLagAIC)]);
disp(['Optimal Lag Length (BIC): ', num2str(optimalLagBIC)]);

%% 2) VAR

% Use the optimal lag length based on AIC
p = 20;  

% Create the lagged matrix
X = createLaggedMatrix(data, p);
X = X(p+1:end, :);  % Remove rows with NaN values
Y = data(p+1:end, :);  % Adjust the range accordingly

% Add a row of ones before X to create the intercept term.
X = [ones(size(X, 1), 1), X];

% Estimate the coefficients by OLS
beta = (X\Y);

% Reshape beta to get the coefficient matrices A_i
k = size(data, 2);  % Number of variables
A = reshape(beta(2:end, :), [k, k, p]);

% Estimate the residuals
residuals = Y - X*beta;

% Estimate the covariance matrix of the residuals
Sigma = cov(residuals);

% Display the estimated coefficient matrices and the covariance matrix of the residuals
disp('Estimated Coefficient Matrices:');
disp(A);
disp('Covariance Matrix of Residuals:');
disp(Sigma);

%% 3) Triangular var

% Cholesky decomposition of the covariance matrix of residuals
L = chol(Sigma, 'lower');

% Triangular VAR model estimation

% The lower triangular matrix L contains the temporal relations between the variables
% We can use the previously estimated A (coefficient matrices) and L to construct the triangular VAR model

% Now, for each lag, we can compute the triangular VAR coefficient matrices :
B = cell(p, 1);  % B will store the triangular VAR coefficient matrices
for i = 1:p
    B{i} = L \ A(:,:,i) / L;  % This solves the equation L * B{i} * L' = A(:,:,i) for B{i}
end

% Display the estimated triangular VAR coefficient matrices
disp('Estimated Triangular VAR Coefficient Matrices:');
for i = 1:p
    disp(['B{', num2str(i), '} = ']);
    disp(B{i});
end

%% 4) Assess the stability of your results : 4-1) Eigenvalues of the Coefficient Matrices

% Eigenvalues of the Companion Matrix
% Create the companion matrix
k = size(A, 1);  % Number of variables
C = [reshape(A, [k, k*p]); eye(k*(p-1)), zeros(k*(p-1), k)];

% Get the eigenvalues of the companion matrix
eigenvalues_C = eig(C);

% Display the eigenvalues and their magnitudes
disp('Eigenvalues of the companion matrix:');
disp(eigenvalues_C);
disp('Magnitudes of the eigenvalues of the companion matrix:');
disp(abs(eigenvalues_C));

% Eigenvalues of the Individual Coefficient Matrices A_i and B_i
% Loop through each lag and compute the eigenvalues of A_i and B_i
for i = 1:p
    % Eigenvalues of A_i
    eigenvalues_Ai = eig(A(:,:,i));
    disp(['Eigenvalues of A{', num2str(i), '}']);
    disp(eigenvalues_Ai);
    disp('Magnitudes of the eigenvalues:');
    disp(abs(eigenvalues_Ai));
    
    % Eigenvalues of B_i
    eigenvalues_Bi = eig(B{i});
    disp(['Eigenvalues of B{', num2str(i), '}']);
    disp(eigenvalues_Bi);
    disp('Magnitudes of the eigenvalues:');
    disp(abs(eigenvalues_Bi));
end

%% 4.2) Residual Analysis
% Check for Zero Mean
% Assuming residuals is your residuals matrix from the VAR model
mean_res = mean(residuals, 1);  % Compute the mean of residuals for each variable
disp('Mean of residuals for each variable:');
disp(mean_res);

%% 4.3) Constant Variance Checking (Homoscedasticity)
figure;
plot(residuals);
title('Residuals Over Time');
legend('OilSupplySurpriseSeries', 'OilSupplyNewsShock', 'CPI', 'IPI'); 

%% 4.4) Impulse Response Functions (IRFs)
% Variables
% 1: Oil Supply Surprise Series
% 2: Oil Supply News Shock
% 3: CPI
% 4: IPI

% Parameters
p = 20;  % Number of lags
k = size(data, 2);  % Number of variables
steps = 10;  % Number of steps ahead to compute the IRFs

% Cholesky decomposition for identifying the shocks
L = chol(Sigma, 'lower');

% Initialize the IRFs matrix
% Dimensions are: number of steps, number of variables, number of shocks
irfs = zeros(steps, k, k);

% Compute the IRFs
for shock = 1:k  % Loop over each shock
    impulse = L(:, shock);  % One-unit shock
    response = impulse;  % Initial response is the shock itself
    irfs(1, :, shock) = response';  % Store the initial response
    for step = 2:steps  % Loop over each step ahead
        % Compute the response at each step
        response = reshape(A, [k, k*p]) * [response; zeros((p-1)*k, 1)];
        irfs(step, :, shock) = response';  % Store the response
    end
end

% Variables Names
var_names = {
    'Oil Supply Surprise Series', 
    'Oil Supply News Shock', 
    'CPI', 
    'IPI'
};

% Specific graphs to plot
graph_specs = {
    [3, 1]; % CPI to Oil Supply Surprise Series
    [3, 2]; % CPI to Oil Supply News Shock
    [4, 1]; % IPI to Oil Supply Surprise Series
    [4, 2]  % IPI to Oil Supply News Shock
};

% Plot the specified IRFs
for i = 1:size(graph_specs, 1)
    var = graph_specs{i}(1);
    shock = graph_specs{i}(2);
    figure;
    plot(1:steps, squeeze(irfs(:, var, shock)));
    title(['Response of ', var_names{var}, ' to Shock in ', var_names{shock}]);
end

%% 5) Create plots of all your results, 5.1) Plotting Coefficient Matrices:
% A and B are the estimated coefficient matrices from VAR and Triangular VAR respectively

% For VAR
figure;
for i = 1:p
    subplot(4,5,i);
    imagesc(A(:,:,i));
    colorbar;
    title(['A_{', num2str(i), '}']);
end
sgtitle('Coefficient Matrices from VAR');

% For Triangular VAR
figure;
for i = 1:p
    subplot(4,5,i);
    imagesc(B{i});
    colorbar;
    title(['B_{', num2str(i), '}']);
end
sgtitle('Coefficient Matrices from Triangular VAR');

%% 5.2) Plotting Covariance Matrix of Residuals:
% Sigma is the covariance matrix of residuals from VAR

figure;
imagesc(Sigma);
colorbar;
title('Covariance Matrix of Residuals');

%% Define a function to create a lagged matrix 
function lagged = createLaggedMatrix(data, lags)
    [T, k] = size(data);
    lagged = NaN(T, k*lags);
    for i = 1:lags
        lagged(i+1:end, ((i-1)*k+1):(i*k)) = data(1:end-i, :);
    end
end