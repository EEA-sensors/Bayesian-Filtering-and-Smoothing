%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Estimate pendulum state with IEKF and IERTS as in Examples 7.10 and 13.6 of
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
    iter = 10;

    m = m0;
    P = P0;
    MM = zeros(size(m,1),length(Y));
    PP = zeros(size(P,1),size(P,2),length(Y));
    for k=1:length(Y)
        f = [m(1)+m(2)*DT; m(2)-g*sin(m(1))*DT];
        F = [1 DT; -g*cos(m(1))*DT 1];
        m = f;
        P = F*P*F' + Q;
    
        mp = m;

        for i=1:iter
            h = sin(m(1));
            H = [cos(m(1)) 0];
            v = Y(k) - h - H * (mp - m);
            S = H*P*H' + R;
            K = P*H'/S;
            m = mp + K*v;
        end
        P = P - K*S*K';
    
        MM(:,k) = m;
        PP(:,:,k) = P;
    end
  
    subplot(2,1,1);
    h = plot(T,Y,'k.',T,X(1,:),'r-',T,MM(1,:),'b--');
    set(h,'Linewidth',5);
    title('IEKF estimate');
    legend('Measurements','True','Estimate');

    subplot(2,1,2);
    h = plot(T,squeeze(PP(1,1,:)),'b--');
    
    rmse_iekf = sqrt(mean((X(1,:)-MM(1,:)).^2))
  
%%
% Smoother
%

    iter = 10;

    Nm0  = m0;
    NP0  = P0;
    NMMS = MM; % Initialize with IEKF (could do with basic EKF)
    NPPS = PP;

    for i=1:iter
        m = m0;
        P = P0;
        for k=1:length(Y)
            if k > 1
                nm = NMMS(:,k-1);
            else
                nm = Nm0;
            end

            F = [1 DT; -g*cos(nm(1))*DT 1];
            f = [nm(1)+nm(2)*DT; nm(2)-g*sin(nm(1))*DT] + F*(m - nm);
            m = f;
            P = F*P*F' + Q;
        
            nm = NMMS(:,k);

            H = [cos(nm(1)) 0];
            h = sin(nm(1)) + H*(m - nm);

            S = H*P*H' + R;
            K = P*H'/S;
            m = m + K*(Y(k) - h);
            P = P - K*S*K';
        
            MM(:,k) = m;
            PP(:,:,k) = P;
        end
        
        ms = m;
        Ps = P;
        MMS = MM;
        PPS = PP;
        for k=size(MM,2)-1:-1:0
            if k == 0
                m = m0;
                P = P0;
                nm = Nm0;
                nP = NP0;
            else
                m = MM(:,k);
                P = PP(:,:,k);
                nm = NMMS(:,k);
                nP = NPPS(:,:,k);
            end
                    
            F = [1 DT; -g*cos(nm(1))*DT 1];
            f = [nm(1)+nm(2)*DT; nm(2)-g*sin(nm(1))*DT] + F*(m - nm);

            mp = f;
            Pp = F*P*F'+Q;
            Ck = P*F'/Pp;

            ms = m + Ck*(ms - mp);
            Ps = P + Ck*(Ps - Pp)*Ck';
            Ps = 0.5 * (Ps + Ps');

            if k == 0
                Nm0 = ms;
                NP0 = Ps;
            else
                MMS(:,k) = ms;
                PPS(:,:,k) = Ps;
            end
        end
        
        NMMS = MMS;
        NPPS = PPS;

        rmse_ieks = sqrt(mean((X(1,:)-MMS(1,:)).^2));
        fprintf('IEKS RMSE %d/%d = %f\n',i,iter,rmse_ieks);
    end

    h = plot(T,Y,'k.',T,X(1,:),'r-',...
             T,MM(1,:),'b--',T,MMS(1,:),'g--');
    set(h,'Linewidth',5);
    title('IEKF and IEKS estimates');
    legend('Measurements','True','IEKF','IEKS');
    
    rmse_ieks = sqrt(mean((X(1,:)-MMS(1,:)).^2))

%%
% Plot the filtering result
%

    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,MM(1,:),'r');
  
    legend('True angle','Measurements','IEKF estimate');
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
    
%%
% Plot the smoothing result
%

    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,MMS(1,:),'r');
  
    legend('True angle','Measurements','IERTSS estimate');
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
   
