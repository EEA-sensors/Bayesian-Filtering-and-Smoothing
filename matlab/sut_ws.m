function [WM, XI, WC] = sut_ws(n,alpha,beta,kappa)
%SUT_WS  Scaled unscented transform weights and unit sigma points
%
% Syntax:
%   [WM, XI, WC] = sut_ws(n,alpha,beta,kappa)
%
% In:
%   n     - Dimensionality of random variable
%   alpha - Transformation parameter  (optional, default 0.5)
%   beta  - Transformation parameter  (optional, default 2)
%   kappa - Transformation parameter  (optional, default 3-n)
%
% Out:
%   WM - Weights for mean calculation
%   XI - Unit sigma points
%   WC - Weights for covariance calculation

    %
    % Check which arguments are there
    %
    if nargin < 1
        error('At least dimensionality n required.');
    end
    if nargin < 2
        alpha = [];
    end
    if nargin < 3
        beta = [];
    end
    if nargin < 4
        kappa = [];
    end
    
    %
    % Apply default values
    %
    if isempty(alpha)
        alpha = 1;
    end
    if isempty(beta)
        beta = 2;
    end
    if isempty(kappa)
        kappa = 3 - n;
    end

    %
    % Compute the weights
    %
    lambda = alpha^2 * (n + kappa) - n;        
    WM = zeros(1,2*n+1);
    WC = zeros(1,2*n+1);
    for j=1:2*n+1
        if j==1
            wm = lambda / (n + lambda);
            wc = lambda / (n + lambda) + (1 - alpha^2 + beta);
        else
            wm = 1 / (2 * (n + lambda));
            wc = wm;
        end
        WM(j) = wm;
        WC(j) = wc;
    end

    %
    % Form the unit sigma points
    %
    XI = sqrt(n + lambda) * [zeros(n,1) eye(n) -eye(n)];
end
