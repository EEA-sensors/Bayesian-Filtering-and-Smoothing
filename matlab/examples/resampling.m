%RESAMPLING  Perform resampling of weighted samples to uniform weighting
%
% Syntax:
%   ind = resampling(W,scheme,M,normalize)
%
% In:
%           W - 1xN vector of positive weights
%      scheme - Resampling scheme 'multinomial', 'stratified', or
%              'systematic' (default: 'multinomial')
%           M - Number of indices to draw (default: length(W))
%   normalize - Whether weight should be normalized, if not, then they
%               are assumed to already sum to 1 (default: true)
%
% Out:
%   ind - 1xM vector of indices.
%
% Description:
%   Perform resampling using a selected method. The function returns a
%   vector of indices such that the resampled version of original (W,X)
%   is given by (1/M,X(ind)). The resampling method options are
%   'multinomial', 'stratified', or 'systematic'. These particular
%   implementations are based on "An introduction to Sequential Monte Carlo"
%   by Nicolas Chopin and Omiros Papaspiliopoulos, Springer, 2020 and they
%   have computational complexities of O(N), where N = length(W).

% Copyright:
%   Copyright (c) 2021 Simo Särkkä
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

function ind = resampling(W,scheme,M,normalize)

    if nargin < 2
        scheme = [];
    end
    if nargin < 3
        M = [];
    end
    if nargin < 4
        normalize = [];
    end
    
    if isempty(scheme)
        scheme = 'multinomial';
    end
    if isempty(M)
        M = length(W);
    end
    if isempty(normalize)
        normalize = true;
    end

    if normalize
        W = W ./ sum(W);
    end    
    
    if strcmp(scheme,'multinomial')
        tmp = cumsum(-log(rand(1,M + 1)));
        us = tmp(1:end-1) / tmp(end);
    elseif strcmp(scheme,'stratified')
        us = (rand(1,M) + (0:(M-1))) / M;
    elseif strcmp(scheme,'systematic')
        us = (rand + (0:(M-1))) / M;
    else
        error('Unknown resampling scheme %s.\n',scheme);
    end
    
    m = 1;
    s = W(1);
    M = length(us);
    ind = zeros(1,M);
    for n=1:M
        while s < us(n)
            m = m + 1;
            s = s + W(m);
        end
        ind(n) = m;
    end
end


