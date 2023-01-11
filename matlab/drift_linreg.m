%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Drifted linear regression demonstration from Chapter 3 of the book
%
% Simo Sarkka and Lennard Svensson (2023), Bayesian Filtering and Smoothing,
% Cambridge University Press. 
%
% See LICENSE provided with the software.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
% Simulate data
%

%    randn('state',12);  % Use this to get book's data
    rng(1, 'twister');
    
    dt = 0.01;
    sd = 0.2;
    t = (0:dt:2);
    x = sin(pi*t);
    y = x + sd*randn(size(x));
    
    h = plot(t,y,'.',t,x,'-');
    
    axis([0 2 -1.5 1.5]);
    
    set(h,'Markersize',7);
    set(h,'LineWidth',3);
    set(h(1),'Color',[0.0 0.0 0.0]);
    set(h(2),'Color',[0.7 0.7 0.7]);
    
    
    h = legend('Measurements','True Signal','Location','NE');
    xlabel('{\it t}');
    
    rmse = sqrt(mean((x - y).^2))
  
%%
% Kalman filter 1
%
    m0 = [0;0];
    P0 = eye(2);
    m = m0;
    P = P0;
    kf1_MM = zeros(size(m0,1),length(y));
    kf1_PP = zeros(size(P0,1),size(P0,1),length(y));
    Q = 0.1*eye(2)*dt;
    for k=1:length(y)
        P = P + Q;
        H = [1 t(k)];
        S = H*P*H'+sd^2;
        K = P*H'/S;
        m = m + K*(y(k)-H*m);
        P = P - K*S*K';
        
        kf1_MM(:,k) = m;
        kf1_PP(:,:,k) = P;
    end
    
    est = kf1_MM(1,:) + t.*kf1_MM(2,:);
    h = plot(t,y,'.',t,x,'-',t,est,'b-');
    
    axis([0 2 -1.5 1.5]);
    
    set(h,'Markersize',7);
    set(h(2),'LineWidth',4);
    set(h(3),'LineWidth',1.5);
    set(h(1),'Color',[0.0 0.0 0.0]);
    set(h(2),'Color',[0.7 0.7 0.7]);
    set(h(3),'Color',[0.0 0.0 0.0]);
 
    h = legend('Measurements','True Signal','Estimate','Location','NE');
    
    xlabel('{\it t}');
    
    rmse = sqrt(mean((x - est).^2))

%%
% Kalman filter 2
%
    m0 = [0;0];
    P0 = eye(2);
    m = m0;
    P = P0;
    kf2_MM = zeros(size(m0,1),length(y));
    kf2_PP = zeros(size(P0,1),size(P0,1),length(y));
    q = 10;
    A = [1 dt; 0 1];
    Q = q * [dt^3/3 dt^2/2; dt^3/3 dt];
    
    for k=1:length(y)
        m = A*m;
        P = A*P*A' + Q;      
        H = [1 0];
        S = H*P*H'+sd^2;
        K = P*H'/S;
        m = m + K*(y(k)-H*m);
        P = P - K*S*K';
        
        kf2_MM(:,k) = m;
        kf2_PP(:,:,k) = P;
    end
    
    est = kf2_MM(1,:);
    h = plot(t,y,'.',t,x,'-',t,est,'b-');
    
    axis([0 2 -1.5 1.5]);
    
    set(h,'Markersize',7);
    set(h(2),'LineWidth',4);
    set(h(3),'LineWidth',1.5);
    set(h(1),'Color',[0.0 0.0 0.0]);
    set(h(2),'Color',[0.7 0.7 0.7]);
    set(h(3),'Color',[0.0 0.0 0.0]);
    
    
    h = legend('Measurements','True Signal','Estimate','Location','NE');
    
    xlabel('{\it t}');
    
    rmse = sqrt(mean((x - est).^2))

%%
% Kalman smoother
%
    ms = kf2_MM(:,end);
    Ps = kf2_PP(:,:,end);
    rts_m = zeros(size(m,1),length(y));
    rts_P = zeros(size(P,1),size(P,2),length(y));
    rts_m(:,end) = ms;
    rts_P(:,:,end) = Ps;
    for k=size(kf2_MM,2)-1:-1:1
        mp = A*kf2_MM(:,k);
        Pp = A*kf2_PP(:,:,k)*A'+Q;
        Ck = kf2_PP(:,:,k)*A'/Pp;
        ms = kf2_MM(:,k) + Ck*(ms - mp);
        Ps = kf2_PP(:,:,k) + Ck*(Ps - Pp)*Ck';
        rts_m(:,k) = ms;
        rts_P(:,:,k) = Ps;
    end
    
    est = rts_m(1,:);
    h = plot(t,y,'.',t,x,'-',t,est,'b-');
    
    axis([0 2 -1.5 1.5]);
    
    set(h,'Markersize',7);
    set(h(2),'LineWidth',4);
    set(h(3),'LineWidth',1.5);
    set(h(1),'Color',[0.0 0.0 0.0]);
    set(h(2),'Color',[0.7 0.7 0.7]);
    set(h(3),'Color',[0.0 0.0 0.0]);
    
    
    h = legend('Measurements','True Signal','Estimate','Location','NE');
    
    xlabel('{\it t}');
    
    rmse = sqrt(mean((x - est).^2))
    
%%
% Prediction
%
    psteps = 50;
    
    PMM = zeros(size(rts_m,1),size(rts_m,2)+psteps);
    PPP = zeros(size(rts_P,1),size(rts_P,2),size(rts_P,3)+psteps);
    ind = size(rts_m,2);
    PMM(:,1:ind) = rts_m;
    PPP(:,:,1:ind) = rts_P;
    
    PT = [t t(end)+(1:psteps)*dt];
    m = rts_m(:,end);
    P = rts_P(:,:,end);
    for k=ind+1:ind+psteps
        m = A*m;
        P = A*P*A' + Q;
        PMM(:,k) = m;
        PPP(:,:,k) = P;    
    end
    
    est = PMM(1,:);
    h = plot(t,y,'.',t,x,'-',PT,est,'b-');
    
    axis([0 PT(end) -1.5 1.5]);
    
    set(h,'Markersize',7);
    set(h(2),'LineWidth',4);
    set(h(3),'LineWidth',1.5);
    set(h(1),'Color',[0.0 0.0 0.0]);
    set(h(2),'Color',[0.7 0.7 0.7]);
    set(h(3),'Color',[0.0 0.0 0.0]);
    
    h = legend('Measurements','True Signal','Estimate','Location','SW');
    
    xlabel('{\it t}');
