%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Estimate cluttered pendulum state with a particle filter as in Examples
% 11.11 and 15.5 of the book
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
    pendulum_sim2;

    %%
    % GHKF
    %
    m = m0;
    P = P0;
    
    n = size(m,1);

    p = 5;
    [W, XI] = gh_ws(n,p);
    
    %
    % Do the filtering
    %  
    MM = zeros(size(m,1),length(Y));
    PP = zeros(size(P,1),size(P,2),length(Y));
    for k=1:length(Y)

        % Form the sigma points for dynamic model
        SX = repmat(m,1,size(XI,2)) + chol(P,'lower') * XI;

        % Propagate through the dynamic model
        HX = [SX(1,:)+SX(2,:)*DT; SX(2,:)-g*sin(SX(1,:))*DT];
  
        % Compute the predicted mean and covariance
        m = zeros(size(m));
        P = zeros(size(P));
        for i=1:size(HX,2)
            m = m + W(i) * HX(:,i);
        end
        for i=1:size(HX,2)
            P = P + W(i) * (HX(:,i) - m) * (HX(:,i) - m)';
        end
        P = P + Q;

        % Form sigma points for measurement step and
        % propagate throught the measurement model
        SX = repmat(m,1,size(XI,2)) + chol(P,'lower') * XI;
        HY = sin(SX(1,:)); 
    
        % Compute the updated mean and covariance
        mu = zeros(size(HY,1),1);
        S  = zeros(size(HY,1),size(HY,1));
        C  = zeros(size(SX,1),size(HY,1));
        for i=1:size(SX,2)
            mu = mu + W(i) * HY(:,i);
        end
        for i=1:size(SX,2)
            S = S + W(i) * (HY(:,i) - mu) * (HY(:,i) - mu)';
            C = C + W(i) * (SX(:,i) - m) * (HY(:,i) - mu)';
        end
        S = S + R;
    
        % Compute the gain and updated mean and covariance  
        K = C/S;
        m = m + K*(Y(k) - mu);
        P = P - K*S*K';
    
        MM(:,k) = m;
        PP(:,:,k) = P;
    end  
    ghkf_MM = MM;
    ghkf_PP = PP;

    clf;
    h = plot(T,Y,'k.',T,X(1,:),'r-',T,MM(1,:),'b--');
    set(h,'Linewidth',5);
    title('GHKF estimate');
    legend('Measurements','True','Estimate');

    rmse_ghkf = sqrt(mean((X(1,:)-MM(1,:)).^2))

    %%
    % GHRTS
    %
    ms = m;
    Ps = P;
    MS = zeros(size(m,1),length(Y));
    PS = zeros(size(P,1),size(P,2),length(Y));
    MMS(:,k) = m;
    PPS(:,:,k) = P;
    for k=size(MM,2)-1:-1:1
        m = MM(:,k);
        P = PP(:,:,k);
    
        % Form the sigma points for dynamic model
        SX = repmat(m,1,size(XI,2)) + chol(P,'lower') * XI;

        % Propagate through the dynamic model
        HX = [SX(1,:)+SX(2,:)*DT; SX(2,:)-g*sin(SX(1,:))*DT];
 
        % Compute the predicted mean and covariance
        % and the cross-covariance D.
        mp = zeros(size(m));
        Pp = zeros(size(P));
        D  = zeros(size(P));
        for i=1:size(HX,2)
            mp = mp + W(i) * HX(:,i);
        end
        for i=1:size(HX,2)
            Pp = Pp + W(i) * (HX(:,i) - mp) * (HX(:,i) - mp)';
            D = D + W(i) * (SX(:,i) - m) * (HX(:,i) - mp)';
        end
        Pp = Pp + Q;

        Gk = D/Pp;
        ms = m + Gk*(ms - mp);
        Ps = P + Gk*(Ps - Pp)*Gk';
        MMS(:,k) = ms;
        PPS(:,:,k) = Ps;
    end
  
    ghrtss_MMS = MMS;
    ghrtss_PPS = PPS;

    h = plot(T,Y,'k.',T,X(1,:),'r-',...
             T,MM(1,:),'b--',T,MMS(1,:),'g--');
    set(h,'Linewidth',5);
    title('GHKF and GHRTS estimates');
    legend('Measurements','True','GHKF','GHRTS');
    
    rmse_ghrts = sqrt(mean((X(1,:)-MMS(1,:)).^2))


    %%
    % Particle filter
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
    
        % Draw indicator
        c = rand(1,size(SX,2)) < cp;
        
        % Compute the weights
        ind0 = find(c == 0);
        ind1 = find(c == 1);
        
        my = sin(SX(1,ind0));
        W(ind0) = exp(-1/(2*R)*(Y(k) - my).^2); % Constant discarded
        W(ind1) = 1/4;
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

    pf_MM = MM;

    h = plot(T,Y,'k.',T,X(1,:),'r-',T,MM(1,:),'b--');
    set(h,'Linewidth',5);
    fprintf('BF estimate.\n');
    legend('Measurements','True','Estimate');
  
    rmse_bf = sqrt(mean((X(1,:)-MM(1,:)).^2))

  
%%
% Backward simulation smoother
%

    % You can comment this out
    rng(1, 'twister');

    NS = 100;
    SM_SSX = zeros(size(m,1),NS,length(Y));
  
    for i=1:NS
        ind = floor(rand * N + 1);
        xn = SSX(:,ind,end);
        SM_SSX(:,i,end) = xn;
        for k=length(Y)-1:-1:1
            SX = SSX(:,:,k);
            mu = [SX(1,:)+SX(2,:)*DT; SX(2,:)-g*sin(SX(1,:))*DT];

            % Evaluate the transition density
            DX = xn-mu;
            E = 0.5 * sum(DX .* (Q \ DX),1) + 0.5 * size(mu,1) * log(2*pi) + 0.5 * log(det(Q));
            W = exp(-E);
            W  = W ./ sum(W);

            % Draw a new sample
            ind = resampling(W,'multinomial',1);
            xn = SSX(:,ind,k);
            SM_SSX(:,i,k) = xn;
        end
        if rem(i,10)==0
            plot(T,X(1,:),'r-',T,squeeze(SM_SSX(1,i,:)),'b--');
            title(sprintf('Doing backward simulation %d/%d\n',i,NS));
            drawnow;
        end
    end

    MMS = squeeze(mean(SM_SSX,2));
    ps_MMS = MMS;
  
    h = plot(T,Y,'k.',T,X(1,:),'r-',T,MMS(1,:),'b--');
    set(h,'Linewidth',5);
    fprintf('BS estimate.\n');
    legend('Measurements','True','Estimate');
  
    rmse_bs = sqrt(mean((X(1,:)-MMS(1,:)).^2))
  
%%
% Plot the filtering result
%

    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,pf_MM(1,:),'r',T,ghkf_MM(1,:),'--');
  
    legend('True angle','Measurements','PF estimate','GHKF estimate');
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 


%%
% Plot the smoothing result
%

    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,ps_MMS(1,:),'r',T,ghrtss_MMS(1,:),'--');
  
    legend('True angle','Measurements','PS estimate','GHRTSS estimate');
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
