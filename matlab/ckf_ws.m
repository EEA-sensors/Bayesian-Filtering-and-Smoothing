function [W, XI] = ckf_ws(n)
    XI = sqrt(n) * [eye(n) -eye(n)];
    W  = ones(1,2*n)/(2*n);
end
