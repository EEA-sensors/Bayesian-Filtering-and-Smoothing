%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pendulum parameter posterior estimation with GHKF and PMCMC
% as in Example 16.17 of the book
%
% Simo Sarkka and Lennart Svensson (2023), Bayesian Filtering and Smoothing,
% 2nd ed., Cambridge University Press.
% 
% See LICENSE provided with the software.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Simulate data
%
    pendulum_sim;

%%
% Estimate pendulum the posterior with GHKF
%

    p = 5; % Order of the method
    n = 2;

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
    W = prod(W1(ind),1); % ND weights

    
    f = @(x) [x(1)+x(2)*DT; x(2)-g*sin(x(1))*DT];
    h = @(x) sin(x(1));

    dR = 0.001;
    %dR = 0.01; % A bit less dense grid
    RR = 0.05:dR:0.15;
    pp = zeros(size(RR));
    for ri=1:length(RR)
        m = m0;
        P = P0;
        
        energy = 0;
        
        MM = zeros(size(m,1),size(Y,2));
        PP = zeros(size(P,1),size(P,2),size(Y,2));
        for k=1:size(Y,2)
            SX = repmat(m,1,size(XI,2)) + chol(P,'lower') * XI;
            m = zeros(size(m));
            for i=1:size(SX,2)
                SX(:,i) = f(SX(:,i));
                m = m + W(i) * SX(:,i);
            end
            P = Q;
            for i=1:size(SX,2)
                P = P + W(i) * (SX(:,i) - m) * (SX(:,i) - m)';
            end
            
            SX = repmat(m,1,size(XI,2)) + chol(P,'lower') * XI;
            SY = zeros(size(Y,1),size(SX,2));
            mu = zeros(size(Y,1),1);
            for i=1:size(SX,2)
                SY(:,i) = h(SX(:,i));
                mu = mu + W(i) * SY(:,i);
            end
            S = RR(ri);
            C = zeros(size(SX,1),size(SY,1));
            for i=1:size(SY,2)
                S = S + W(i) * (SY(:,i) - mu) * (SY(:,i) - mu)';
                C = C + W(i) * (SX(:,i) - m) * (SY(:,i) - mu)';
            end
            
            v = Y(:,k) - mu;
            K = C / S;
            m = m + K * v;
            P = P - K * S * K';
            
            energy = energy + 0.5 * log(det(2*pi*S)) + 0.5 * v'*(S\v);
            
            MM(:,k) = m;
            PP(:,:,k) = P;
        end

        pp(ri) = energy;
    end
    
    pp = pp - min(pp);
    pp = exp(-pp);
    pp = pp ./ sum(pp) / dR;
    plot(RR,pp);
    
%%
% PMCMC
%

    % You can comment this out
    rng(1, 'twister');

    N = 1000; % PF samples
    %N = 100; % PF samples
    nmc = 10000; % MCMC samples
    %nmc = 1000; % A bit less

    QL = chol(Q,'lower');

    
    doram = 1;
    samp = [];
    as = 0.234;
    S  = 0.25;
    
    theta = -2;

    en = 0;
    accepted = 0;
    
    for j=1:nmc

        %
        % Draw candidate
        %
        r = randn;
        new_theta = theta + S * r;

        %
        % Evalaute energy with PF
        %        
        m = m0;
        P = P0;
        
        %
        % Initial sample set
        %
        SX = repmat(m,1,N) + chol(P,'lower') * randn(size(m,1), N);
  
        %
        % Do the filtering and store the histories
        %  
        MM  = zeros(size(m,1),length(Y));
        SSX = zeros(size(m,1),N,length(Y)); % filter history
        
        new_en = 0;
        for k=1:length(Y)

            % Propagate through the dynamic model
            SX = [SX(1,:)+SX(2,:)*DT; SX(2,:)-g*sin(SX(1,:))*DT];
    
            % Add the process noise
            SX = SX + QL * randn(size(SX));
    
            % Compute the weights
            my = sin(SX(1,:));
            W  = 1/sqrt(2*pi*exp(new_theta)) * exp(-1/(2*exp(new_theta))*(Y(k) - my).^2);
            lh = mean(W);
            new_en = new_en - log(lh);
            W  = W ./ sum(W);
    
            % Do resampling
            ind = resampstr(W);
            SX  = SX(:,ind);      
    
            %if rem(k,100)==0
            %    fprintf('%d/%d\n',k,length(Y));
            %end
        end
        
        %
        % Accept or reject
        %
        if j == 1
            en = new_en;
        end
        a = min(1,exp(en - new_en));
        u = rand;

        if u <= a
            theta = new_theta;
            en    = new_en;
            accepted = accepted + 1;
%            fprintf('Accepted: %f\n',theta);
        else
%            fprintf('Rejected\n');
        end
        samp = [samp theta];

        %
        % Adapt
        %
        if doram && j > 10
            nu = j^(-0.9);
            S = sqrt(S * (1 + nu * (a - as) * r*r/r^2) * S);
        end
        
        if rem(j,10) == 0
            fprintf('%d/%d Acceptance rate %f\n',j,nmc,accepted/j);        
            subplot(2,1,1);
            plot(samp);
            subplot(2,1,2);
            hist(exp(samp),10)
            drawnow;
        end
    end
    
%%
% Plotting
%

    [HN,HX] = hist(exp(samp(100:end)),100);
    HN = HN / sum(HN) / (HX(2)-HX(1));
    HN = [0 HN 0];
    HX = [HX(1) HX HX(end)];

    clf;
    h = fill(HX,HN,'g');
    set(h,'FaceColor',0.5 * [1 1 1]);
    set(h,'EdgeColor',0.5 * [1 1 1]);
    hold on;
    h = plot(RR,pp,'-');
    set(h,'Color',0.0 * [1 1 1]);
    set(h,'LineWidth',1)
    ax = axis;
    h = plot([R R],[0 ax(end)],'-');
    set(h,'Color',0.7 * [1 1 1]);
    set(h,'LineWidth',2)

    xlabel('{\it R}');
    ylabel('{\it p}({\it{R}} |{\it y_{{\rm{1}}:T}} )')

    legend('PMCMC histogram','Gaussian filter estimate',...
           'True parameter value');
