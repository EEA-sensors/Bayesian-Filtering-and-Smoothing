%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Estimate pendulum state with PF and BS-PS as in Examples 7.1
% and 11.1 of the book
%
% Simo Sarkka (2013), Bayesian Filtering and Smoothing,
% Cambridge University Press. 
%
% Last updated: $Date: 2013/08/26 12:58:41 $.
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Simulate data
%
    pendulum_sim;

%%
% Filter
%

    m = m0;
    P = P0;
    N = 1000;

    QL = chol(Q,'lower');
    
    %
    % Initial sample set
    %
    chol_P = chol(P, 'lower');
    SX = m + chol_P*randn2(2, N);
  
    %
    % Do the filtering and store the histories
    %  
    MM  = zeros(size(m,1),length(Y));
    SSX = zeros(size(m,1),N,length(Y)); % filter history
    for k=1:length(Y)

        % Propagate through the dynamic model
        SX = [SX(1,:)+SX(2,:)*DT; SX(2,:)-g*sin(SX(1,:))*DT];
        % Add the process noise
        SX = SX + QL * randn2(size(SX));

        % Compute the weights
        my = sin(SX(1,:));
        W  = exp(-1/(2*R)*(Y(k) - my).^2); % Constant discarded
        W  = W ./ sum(W);
    
        % Do resampling
        ind = resampstr(W);
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
% Backward simulation smoother
%

%%
% Backward simulation smoother
%

    NS = 100;
    SM_SSX = zeros(size(m,1),NS,length(Y));
    W = 1 / N + zeros(1, N);
    ind = categ_rnd(W, NS);
    SM_s = SSX(:,ind,k);
    SM_SSX(:,:,k) = SM_s;
    
    
    for k=length(Y)-1:-1:1
        SX = SSX(:,:,k);
        mu = [SX(1,:)+SX(2,:)*DT; SX(2,:)-g*sin(SX(1,:))*DT];
        for i=1:NS
            SM = SM_s(:, i);
            W = mvnpdf(SM',mu',Q)';
            W  = W ./ sum(W);
            ind = randsample(N, 1, true, W);
            SM = SX(:,ind);
         
            SM_s(:, i) = SM;
            SM_SSX(:,i,k) = SM;
        end
    end
    
    MMS = squeeze(mean(SM_SSX,2));
  
    h = plot(T,Y,'k.',T,X(1,:),'r-',T,MMS(1,:),'b--');
    set(h,'Linewidth',5);
    fprintf('BS estimate.\n');
    legend('Measurements','True','Estimate');
  
    rmse_bs = sqrt(mean((X(1,:)-MMS(1,:)).^2))
  
%%
% Plot the filtering result
%

    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,MM(1,:),'r');
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
    legend('True Angle','Measurements','PF Estimate');
    
%%
% Plot the final smoothing result
%

    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,MMS(1,:),'r');
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
    legend('True Angle','Measurements','PS Estimate');
    
