function [Fx2,Fxx2] = check_ders(x,f,Fx,Fxx,h)
    if nargin < 4
       Fxx = {};
    end
    if nargin < 5
       h = [];
    end

    if isempty(h)
        h = 1e-6;
    end

    thrs = [1e-6 1e-4 1e-2 Inf];
    strs = {'looks great','quite fine','barely passes','failed'};
    
    n = size(x,1);

    %
    % Evaluate the function
    %
    f0 = f(x);
    d = length(f0);

    fprintf('** Checking derivative of %d -> %d function:\n',n,d);
    disp(f);
    
    %
    % Evaluate and check Fx
    %    
    if ~isnumeric(Fx)
      Fx1 = Fx(x);
    else
      Fx1 = Fx;
    end
    
    Fx2 = zeros(size(Fx1));
    
    for i=1:n
        dx = zeros(size(x));
        dx(i) = h;
        Fx2(:,i) = (f(x + dx) - f0) / h;
    end
    
    err = norm(Fx1 - Fx2);
    ind = find(err < thrs, 1);
        
    fprintf('Checking Fx, error %g, %s.\n',err,strs{ind});
    
    if ind == length(thrs)
        analytical = Fx1
        numerical  = Fx2
        error('Derivative check for Fx failed.');
    end

    if ~isempty(Fxx)
        %
        % Check Fxx
        %
        if isa(Fxx,'function_handle')
            Fxx1 = Fxx(x);
            Fxx2 = Fxx1;
        else        
            Fxx1 = cell(1,d);
            Fxx2 = cell(1,d);
            for i=1:d
                if ~isnumeric(Fxx{i})
                    Fxx1{i} = Fxx{i}(x);
                else
                    Fxx1{i} = Fxx{i};
                end
                Fxx2{i} = zeros(n);
            end
        end
        
        if ~isnumeric(Fx)
            Fx0 = Fx(x);
        else
            Fx0 = Fx;
        end
        
        for i=1:n 
            dx = zeros(size(x));
            dx(i) = h;
            if ~isnumeric(Fx)
                dFx_dxi = (Fx(x + dx) - Fx0) / h;
            else
                dFx_dxi = (Fx2 - Fx0) / h; % todo replace with: df2 = @(x) ( f(x+2*h) - 2*f(x+h) + f(x))/h^2;
            end
            for j=1:d
                Fxx2{j}(i,:) = dFx_dxi(j,:);
            end
        end
        
        for i=1:d
            err = norm(Fxx1{i} - Fxx2{i});
            ind = find(err < thrs, 1);
        
            fprintf('Checking Fxx{%d}, error %g, %s.\n',i,err,strs{ind});
    
            if ind == length(thrs)
                analytical = Fxx1{i}
                numerical  = Fxx2{i}
                error('Derivative check for Fxx failed.');
            end
        end
    else
        Fxx2 = {};
    end
end
