%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Estimate pendulum state with EKF and ERTS as in Examples 7.6 and 13.2 of
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
% Check the derivatives that we have computed
%
    f_fun = @(m) [m(1)+m(2)*DT; m(2)-g*sin(m(1))*DT];
    F_fun = @(m) [1 DT; -g*cos(m(1))*DT 1];

    x = randn(2,1);
    check_ders(x, f_fun, F_fun);

    h_fun = @(m) sin(m(1));
    H_fun = @(m) [cos(m(1)) 0];
    check_ders(x, h_fun, H_fun);


%%
% Filter
%

    m = m0;
    P = P0;
    MM = zeros(size(m,1),length(Y));
    PP = zeros(size(P,1),size(P,2),length(Y));
    for k=1:length(Y)
        f = [m(1)+m(2)*DT; m(2)-g*sin(m(1))*DT];
        F = [1 DT; -g*cos(m(1))*DT 1];
        m = f;
        P = F*P*F' + Q;
    
        h = sin(m(1));
        H = [cos(m(1)) 0];
        S = H*P*H' + R;
        K = P*H'/S;
        m = m + K*(Y(k) - h);
        P = P - K*S*K';
    
        MM(:,k) = m;
        PP(:,:,k) = P;
    end
  
    subplot(2,1,1);
    h = plot(T,Y,'k.',T,X(1,:),'r-',T,MM(1,:),'b--');
    set(h,'Linewidth',5);
    title('EKF estimate');
    legend('Measurements','True','Estimate');

    subplot(2,1,2);
    h = plot(T,squeeze(PP(1,1,:)),'b--');
    
    
    rmse_ekf = sqrt(mean((X(1,:)-MM(1,:)).^2))
  
%%
% Smoother
%

    ms = m;
    Ps = P;
    MMS = zeros(size(m,1),length(Y));
    PPS = zeros(size(P,1),size(P,2),length(Y));
    MMS(:,end) = m;
    PPS(:,:,end) = P;
    for k=size(MM,2)-1:-1:1
        m = MM(:,k);
        P = PP(:,:,k);
        f = [m(1)+m(2)*DT; m(2)-g*sin(m(1))*DT];
        F = [1 DT; -g*cos(m(1))*DT 1];
    
        mp = f;
        Pp = F*P*F'+Q;
        Gk = P*F'/Pp;
        ms = m + Gk*(ms - mp);
        Ps = P + Gk*(Ps - Pp)*Gk';
        MMS(:,k) = ms;
        PPS(:,:,k) = Ps;
    end

    h = plot(T,Y,'k.',T,X(1,:),'r-',...
             T,MM(1,:),'b--',T,MMS(1,:),'g--');
    set(h,'Linewidth',5);
    title('EKF and ERTS estimates');
    legend('Measurements','True','EKF','ERTS');
    
    rmse_erts = sqrt(mean((X(1,:)-MMS(1,:)).^2))

%%
% Plot the filtering result
%

    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,MM(1,:),'r');
  
    legend('True angle','Measurements','EKF estimate');
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
    

%%
% Plot the smoothing result
%

    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,MMS(1,:),'r');
  
    legend('True angle','Measurements','ERTSS estimate');
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
    