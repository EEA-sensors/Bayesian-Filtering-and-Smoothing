function U = sym_set( n, gen )
% U = sym_set( n, gen )

    if nargin < 3
        nonzero = 0;
    end

    if isempty(gen)
        if nonzero
            U = [];
        else
            U = zeros(n,1);
        end
    else
        U = [];
        
        for i=1:n
            u = zeros(n,1);
            u(i) = gen(1);

            if length(gen) > 1
                if abs(gen(1) - gen(2)) < eps
                    V = sym_set(n-i, gen(2:end));
                    for j=1:size(V,2)
                        u(i+1:end) = V(:,j);
                        U = [U u -u];
                    end
                else
                    V = sym_set(n-1, gen(2:end));
                    for j=1:size(V,2)
                        u([1:i-1 i+1:end]) = V(:,j);
                        U = [U u -u];
                    end
                end
            else
                U = [U u -u];
            end
        end
    end
end

