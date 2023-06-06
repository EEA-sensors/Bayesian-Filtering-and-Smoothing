% Syntax:
%   [W,XI] = gh_ws(n,p)
%
% Description:
%   Generate Gauss-Hermite cubature order p based
%   weights W and sigma points SX such that
%
%     int g(x) N(x | 0,I) dx =~ sum W(i) g(SX(i))

function [ W, XI, xi1, W1 ] = gh_ws(n,p)

    %
    % Form Probabilists' Hermite polynomials of
    % order p-1 and p
    %
    Hpm = 1;
    Hp  = [1 0];
    for i=1:p-1
        tmp = Hp;
        Hp = [Hp 0] - [0 0 i*Hpm];
        Hpm = tmp;
    end

    %
    % Single dimensional weights and points
    %
    xi1 = roots(Hp)';
    W1 = factorial(p)./(p^2*polyval(Hpm,xi1).^2);
     
    %Generate all p^n collections of indexes by
    %transforming numbers 0...p^n-1) into p-base system
    %and by adding 1 to each digit
    num = 0:(p^n-1);
    ind = zeros(n,p^n);
    for i=1:n
        ind(i,:) = rem(num,p)+1;
        num = floor(num / p);
    end

    %Form the unit sigma points and weights
    XI = xi1(ind);
    W  = prod(W1(ind),1); % ND weights
