%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Estimate pendulum state with a particle filter as in Example 11.10 of
% the book
%
% Simo Sarkka and Lennart Svensson (2023), Bayesian Filtering and Smoothing,
% 2nd ed., Cambridge University Press.
% 
% See LICENSE provided with the software.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%
    % Simulate data
    %
    pendulum_sim;

    %%
    % Filter
    %

    % You can comment this out
    rng(1, 'twister');

    m = m0;
    P = P0;
    N = 10000;

    QL = chol(Q,'lower');
    
    %
    % Initial sample set
    %
    SX = repmat(m,1,N) + chol(P,'lower') * randn(size(m,1), N);

    %
    % Do the filtering and store the histories
    %  
    MM  = zeros(size(m,1),length(Y));
    SSX = zeros(size(m,1),N,length(Y)); % filter history
    for k=1:length(Y)

        % Propagate through the dynamic model
        SX = [SX(1,:)+SX(2,:)*DT; SX(2,:)-g*sin(SX(1,:))*DT];
    
        % Add the process noise
        SX = SX + QL * randn(size(SX));
    
        % Compute the weights
        my = sin(SX(1,:));
        W  = exp(-1/(2*R)*(Y(k) - my).^2); % Constant discarded
        W  = W ./ sum(W);
    
        % Do resampling
        ind = resampling(W,'stratified');
        SX   = SX(:,ind);      
    
        SSX(:,:,k) = SX;
        % Mean estimate
        m = mean(SX,2);
    
        MM(:,k) = m;
        if rem(k,100)==0
            fprintf('%d/%d\n',k,length(Y));
        end
    end

    h = plot(T,Y,'k.',T,X(1,:),'r-',T,MM(1,:),'b--');
    set(h,'Linewidth',5);
    fprintf('BF estimate.\n');
    legend('Measurements','True','Estimate');
  
    rmse_bf = sqrt(mean((X(1,:)-MM(1,:)).^2))

  
    %%
    % Plot the final filtering figure
    %
    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,MM(1,:),'r');
  
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
    legend('True angle','Measurements','PF estimate');
