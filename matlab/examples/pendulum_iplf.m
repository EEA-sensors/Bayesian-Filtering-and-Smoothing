%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Estimate pendulum state with IEKF and IERTS as in Examples 10.13 and 14.23
% of the book
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
% Select the sigma points
%
    gf_type = 1;

    if gf_type == 1
        name = 'ghkf'
        p = 3;
        [W,XI] = gh_ws(size(m0,1),p);
    elseif gf_type == 2
        name = 'ckf'
        [W,XI] = ckf_ws(size(m0,1));
    elseif gf_type == 3
        name = 'ukf'
        [W,XI] = sut_ws(size(m0,1));
    elseif gf_type == 4
        name = 'ukf5'
        [W,XI] = ut5_ws(size(m0,1));
    end


%%
% Filter
%

    iter = 10;

    m = m0;
    P = P0;
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

        mp = m;
        Pp = P;

        for it=1:iter
            % Compute the required statistics
            SX = repmat(m,1,size(XI,2)) + chol(P,'lower') * XI;
            HY = sin(SX(1,:)); 

            mup = zeros(size(HY,1),1);
            Py  = zeros(size(HY,1),size(HY,1));
            Pxy = zeros(size(SX,1),size(HY,1));
            for i=1:size(SX,2)
                mup = mup + W(i) * HY(:,i);
            end
            for i=1:size(SX,2)
                Py = Py + W(i) * (HY(:,i) - mup) * (HY(:,i) - mup)';
                Pxy = Pxy + W(i) * (SX(:,i) - m) * (HY(:,i) - mup)';
            end
            Py = Py + R;

            H = Pxy' / P;
            b = mup -  H * m;
            O = Py - H * P * H';

            S = H * Pp * H' + O;
            K = Pp * H' / S;
            m = mp + K * (Y(:,k) - H * mp - b);
            P = Pp - K * S * K';
            P = 0.5 * (P + P');
        end
        MM(:,k) = m;
        PP(:,:,k) = P;
    end
  
    subplot(2,1,1);
    h = plot(T,Y,'k.',T,X(1,:),'r-',T,MM(1,:),'b--');
    set(h,'Linewidth',5);
    title('IPLF estimate');
    legend('Measurements','True','Estimate');

    subplot(2,1,2);
    h = plot(T,squeeze(PP(1,1,:)),'b--');
    
    rmse_iplf = sqrt(mean((X(1,:)-MM(1,:)).^2))
  
%%
% Smoother
%

    iter = 10;

    Nm0  = m0;
    NP0  = P0;
    NMMS = MM; % Initialize with IPLF
    NPPS = PP;

    nx = size(m0,1);
    ny = size(Y,1);

    As = zeros(nx,nx,size(Y,2));
    as = zeros(nx,size(Y,2));
    Ls = zeros(nx,nx,size(Y,2));
    
    Hs = zeros(ny,nx,size(Y,2));
    bs = zeros(ny,size(Y,2));
    Os = zeros(ny,ny,size(Y,2));

    MMi = zeros(size(MM));
    PPi = zeros(size(PP));
    MMSi = zeros(size(MM));
    PPSi = zeros(size(PP));

    for it=1:iter

        % Linearize the whole thing
        for k=1:length(Y)
            if k > 1
                m = NMMS(:,k-1);
                P = NPPS(:,:,k-1);
            else
                m = Nm0;
                P = NP0;
            end

            SX = repmat(m,1,size(XI,2)) + chol(P,'lower') * XI;
            HX = [SX(1,:)+SX(2,:)*DT; SX(2,:)-g*sin(SX(1,:))*DT];
      
            mum = zeros(size(m));
            for i=1:size(HX,2)
                mum = mum + W(i) * HX(:,i);
            end
            Px  = zeros(size(P));
            Pxx = zeros(size(P));
            for i=1:size(HX,2)
                Px = Px + W(i) * (HX(:,i) - mum) * (HX(:,i) - mum)';
                Pxx = Pxx + W(i) * (SX(:,i) - m) * (HX(:,i) - mum)';
            end
            Px = Px + Q;
            
            A = Pxx' / P;
            a = mum - A * m;
            L = Px - A * P * A';

            % Note that we store A_k to As(:,:,k+1) etc.
            As(:,:,k) = A;
            as(:,k) = a;
            Ls(:,:,k) = L;
        end

        for k=1:length(Y)
            m = NMMS(:,k);
            P = NPPS(:,:,k);

            SX = repmat(m,1,size(XI,2)) + chol(P,'lower') * XI;
            HY = sin(SX(1,:)); 

            mup = zeros(size(HY,1),1);
            Py  = zeros(size(HY,1),size(HY,1));
            Pxy = zeros(size(SX,1),size(HY,1));
            for i=1:size(SX,2)
                mup = mup + W(i) * HY(:,i);
            end
            for i=1:size(SX,2)
                Py = Py + W(i) * (HY(:,i) - mup) * (HY(:,i) - mup)';
                Pxy = Pxy + W(i) * (SX(:,i) - m) * (HY(:,i) - mup)';
            end
            Py = Py + R;

            H = Pxy' / P;
            b = mup - H * m;
            O = Py - H * P * H';

            Hs(:,:,k) = H;
            bs(:,k) = b;
            Os(:,:,k) = O;
        end

        % Run affine filter and smoother
        m = m0;
        P = P0;
        for k=1:length(Y)
            A = As(:,:,k);
            a = as(:,k);
            L = Ls(:,:,k);
            
            m = A * m + a;
            P = A * P * A' + L;
        
            H = Hs(:,:,k);
            b = bs(:,k);
            O = Os(:,:,k);

            mu = H * m + b;
            S = H * P * H' + O;
            K = P * H' / S;
            m = m + K * (Y(:,k) - mu);
            P = P - K * S * K';
        
            MMi(:,k) = m;
            PPi(:,:,k) = P;
        end
        
        ms = m;
        Ps = P;
        MMSi(:,end) = m;
        PPSi(:,:,end) = P;
        for k=size(MMi,2)-1:-1:0
            if k == 0
                m = m0;
                P = P0;
            else
                m = MMi(:,k);
                P = PPi(:,:,k);
            end

            A = As(:,:,k+1);
            a = as(:,k+1);
            L = Ls(:,:,k+1);

            mp = A * m + a;
            Pp = A * P * A' + L;
            G  = P * A' / Pp;

            ms = m + G * (ms - mp);
            Ps = P + G * (Ps - Pp) * G';
            Ps = 0.5 * (Ps + Ps');

            if k == 0
                Nm0 = ms;
                NP0 = Ps;
            else
                MMSi(:,k) = ms;
                PPSi(:,:,k) = Ps;
            end
        end
        
        NMMS = MMSi;
        NPPS = PPSi;

        MMS = MMSi;
        PPPS = PPSi;

        rmse_ipls = sqrt(mean((X(1,:)-MMSi(1,:)).^2));
        fprintf('IPLS RMSE %d/%d = %f\n',it,iter,rmse_ipls);
    end

    h = plot(T,Y,'k.',T,X(1,:),'r-',...
             T,MMi(1,:),'b--',T,MMSi(1,:),'g--');
    set(h,'Linewidth',5);
    title('IPLF and IPLS estimates');
    legend('Measurements','True','IPLF','IPLS');
    
    rmse_ipls = sqrt(mean((X(1,:)-MMSi(1,:)).^2))

%%
% Plot the filtering result
%

    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,MM(1,:),'r');
  
    legend('True angle','Measurements','IPLF estimate');
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
    

%%
% Plot the smoothing result
%

    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,MMS(1,:),'r');
  
    legend('True angle','Measurements','IPLS estimate');
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
