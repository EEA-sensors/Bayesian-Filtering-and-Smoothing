%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Track car state with Kalman filter and Rauch-Tung-Striebel
% smoother as in Examples 6.8 and 12.4 of the book
%
% Simo Sarkka and Lennart Svensson (2023), Bayesian Filtering and Smoothing,
% Cambridge University Press. 
%
% See LICENSE provided with the software.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Set the parameters
%
    q = 1;
    dt = 0.1;
    s = 0.5;
    A = [1 0 dt 0;
        0 1 0 dt;
        0 0 1 0;
        0 0 0 1];
    Q = q*[dt^3/3 0 dt^2/2 0;
           0 dt^3/3 0 dt^2/2;
           dt^2/2 0 dt 0;
           0 dt^2/2 0 dt];
    
    H = [1 0 0 0;
         0 1 0 0];
    R = s^2*eye(2);
    m0 = [0;0;1;-1];
    P0 = eye(4);

%%
% Simulate data
%

    randn('state',33);  % Use this to get the book's data
%    rng(14, 'twister');

    steps = 100;
    X = zeros(size(A,1),steps);
    Y = zeros(size(H,1),steps);
    x = m0;
    for k=1:steps
        q = chol(Q)'*randn(size(A,1),1);
        x = A*x + q;
        y = H*x + s*randn(2,1);
        X(:,k) = x;
        Y(:,k) = y;
    end
    
    plot(X(1,:),X(2,:),'-',Y(1,:),Y(2,:),'.',X(1,1),X(2,1),'*');
    legend('Trajectory','Measurements');
    xlabel('{\it x}_1');
    ylabel('{\it x}_2');

    % Store the data in JSON
    %{
    jstruct = struct('T',T,'X',X,'Y',Y);
    json = jsonencode(jstruct,'PrettyPrint',true);
    filename = 'car_track.json';
    fid = fopen(filename, 'w');
    fwrite(fid, json);
    fclose(fid);
    %}

%%
% Kalman filter
%
    m = m0;
    P = P0;
    kf_m = zeros(size(m,1),size(Y,2));
    kf_P = zeros(size(P,1),size(P,2),size(Y,2));
    for k=1:size(Y,2)
        m = A*m;
        P = A*P*A' + Q;
        
        S = H*P*H' + R;
        K = P*H'/S;
        m = m + K*(Y(:,k) - H*m);
        P = P - K*S*K';
        
        kf_m(:,k) = m;
        kf_P(:,:,k) = P;
    end
    
    rmse_raw = sqrt(mean(sum((Y - X(1:2,:)).^2,1)))
    rmse_kf = sqrt(mean(sum((kf_m(1:2,:) - X(1:2,:)).^2,1)))
    
    clf;
    h=plot(X(1,:),X(2,:),'-',Y(1,:),Y(2,:),'o',...
        kf_m(1,:),kf_m(2,:),'-');
    legend('True Trajectory','Measurements','Filter Estimate');
    xlabel('{\it x}_1');
    ylabel('{\it x}_2');

%%
% RTS smoother
%
    ms = m;
    Ps = P;
    rts_m = zeros(size(m,1),size(Y,2));
    rts_P = zeros(size(P,1),size(P,2),size(Y,2));
    rts_m(:,end) = ms;
    rts_P(:,:,end) = Ps;
    for k=size(kf_m,2)-1:-1:1
        mp = A*kf_m(:,k);
        Pp = A*kf_P(:,:,k)*A'+Q;
        Gk = kf_P(:,:,k)*A'/Pp; 
        ms = kf_m(:,k) + Gk*(ms - mp);
        Ps = kf_P(:,:,k) + Gk*(Ps - Pp)*Gk';
        rts_m(:,k) = ms;
        rts_P(:,:,k) = Ps;
    end
    
    rmse_rts = sqrt(mean(sum((rts_m(1:2,:) - X(1:2,:)).^2,1)))
    
    clf;
    h=plot(X(1,:),X(2,:),'-',Y(1,:),Y(2,:),'o',...
        rts_m(1,:),rts_m(2,:),'-');
    legend('True Trajectory','Measurements','Smoother Estimate');
    xlabel('{\it x}_1');
    ylabel('{\it x}_2');

    