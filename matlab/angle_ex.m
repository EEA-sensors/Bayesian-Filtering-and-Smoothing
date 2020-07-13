%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Supplemental Matlab code for Exercise 5.5 in the book
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

%% Generate data

  % Lock seed
  randn('state',123);

  % Implement RMSE (true and estimate)
  rmse = @(X,EST) sqrt(mean(sum((X-EST).^2)));
  
  % Create a bit curved trajectory and angle
  % measurements from two sensors
  S1 = [-1.5;0.5]; % Position of sensor 1
  S2 = [1;1];      % Position of sensor 2
  sd = 0.05;       % Standard deviation of measurements
  dt = 0.01;       % Sampling period
  x0 = [0;0;1;0];  % Initial state

  a = zeros(1,500);
  a(1,50:100)  = pi/2/51/dt + 0.01*randn(1,51);
  a(1,200:250) = pi/2/51/dt + 0.01*randn(1,51);
  a(1,350:400) = pi/2/51/dt + 0.01*randn(1,51);
  x = x0;
  t = 0;
  X = [];
  Theta = [];
  T = [];
  for i=1:500
    F = [0 0  1    0;...
         0 0  0    1;...
         0 0  0   a(i);...
         0 0 -a(i) 0];
    x = expm(F*dt)*x;
    y1 = atan2(x(2)-S1(2), x(1)-S1(1)) + sd * randn;
    y2 = atan2(x(2)-S2(2), x(1)-S2(1)) + sd * randn;
    t  = t + dt;
    X = [X x];
    T = [T t];
    Theta = [Theta [y1;y2]];
  end
  steps = size(Theta,2);


%% Dynamic model

  % Parameters of the dynamic model
  qc = 0.1;
  
  % This is the transition matrix
  A  = [1 0 dt 0;
        0 1 0 dt;
        0 0 1 0;
        0 0 0 1];

  % This is the process noise covariance
  Q = [qc*dt^3/3 0 qc*dt^2/2 0;
       0 qc*dt^3/3 0 qc*dt^2/2;
       qc*dt^2/2 0 qc*dt 0;
       0 qc*dt^2/2 0 qc*dt];
  
   
%% Baseline solution
   
  % Baseline solution. The estimates
  % of x_k are stored as columns of
  % the matrix EST1.

  %if 0 %% <--- Uncomment to disable
  
  fprintf('Running base line solution.\n');
    
  % Initialize to true value
  m1 = x0;     
  EST1 = zeros(4,steps);
  
  % Set up figure
  figure(1); clf
  
  % Loop through steps
  for k=1:steps
      
    % Compute crossing of the measurements
    dx1 = cos(Theta(1,k));
    dy1 = sin(Theta(1,k));
    dx2 = cos(Theta(2,k));
    dy2 = sin(Theta(2,k));
    d = [dx1 dx2; dy1 dy2]\[S2(1)-S1(1);S2(2)-S1(2)];
    
    % Crossing
    cross_xy = S1 + [dx1;dy1]*d(1);
    
    % Compute estimate
    m1(3:4) = [0;0];
    m1(1:2) = cross_xy;
    EST1(:,k) = m1;
    
    % Animate
    if rem(k,10) == 1
      len = 3;
      dx1 = len*cos(Theta(1,k));
      dy1 = len*sin(Theta(1,k));
      dx2 = len*cos(Theta(2,k));
      dy2 = len*sin(Theta(2,k));
      clf;
      plot(X(1,:),X(2,:),'r-',...
           m1(1),m1(2),'bo',...
           EST1(1,1:k),EST1(2,1:k),'b--',...
           [S1(1);S1(1)+dx1],[S1(2);S1(2)+dy1],'k--',...
           [S2(1);S2(1)+dx2],[S2(2);S2(2)+dy2],'k--');
      axis([-2 2 -2.5 1.5]);
      
      % Pause and draw
      drawnow;
      pause(.1)
    end
  end
    
  % Compute error
  err1 = rmse(X,EST1)

  % Plot baseline
  figure(1); clf

    plot(X(1,:),X(2,:),'--',...
         EST1(1,:),EST1(2,:),'-',...
         S1(1),S1(2),'kx',S2(1),S2(2),'ko')
    legend('True trajectory','Baseline estimate','Sensor 1','Sensor 2');
    xlabel('x'); ylabel('y'); title('\bf Baseline Solution')
    axis([-2 2 -2.5 1.5]);
  
    fprintf('This is the BL solution. Press enter.\n');
    pause;
  
  %end %% <--- Uncomment to disable
 
%% EKF solution
  
  % EKF solution. The estimates
  % of x_k are stored as columns of
  % the matrix EST2.
  
  %if 0 %% <--- Uncomment to disable

  fprintf('Running EKF solution.\n');

  m2 = x0;            % Initialize to true value
  P2 = eye(4);        % Some uncertainty
  R  = sd^2*eye(2);   % The joint covariance
  EST2 = zeros(4,steps);
  
  % Set up figure
  figure(1); clf
  
  % Loop through steps
  for k=1:steps
    
    % Compute estimate here
    warning('You should implement EKF here.')
    m2 = [0;0;0;0];
    EST2(:,k) = m2;
    
    % Animate
    if rem(k,10) == 1
      len = 3;
      dx1 = len*cos(Theta(1,k));
      dy1 = len*sin(Theta(1,k));
      dx2 = len*cos(Theta(2,k));
      dy2 = len*sin(Theta(2,k));
      clf;
      plot(X(1,:),X(2,:),'r-',...
           m2(1),m2(2),'bo',...
           EST2(1,1:k),EST2(2,1:k),'b--',...
           [S1(1);S1(1)+dx1],[S1(2);S1(2)+dy1],'k--',...
           [S2(1);S2(1)+dx2],[S2(2);S2(2)+dy2],'k--');
      axis([-2 2 -2.5 1.5]);
      
      % Pause and draw
      drawnow;
      pause(.1)
    end
  end

  % Compute error
  err2 = rmse(X,EST2)
  
  % Plot EKF
  figure(1); clf

    plot(X(1,:),X(2,:),'--',...
         EST2(1,:),EST2(2,:),'-',...
         S1(1),S1(2),'kx',S2(1),S2(2),'ko')
    legend('True trajectory','EKF estimate','Sensor 1','Sensor 2');
    xlabel('x'); ylabel('y'); title('\bf EKF Solution')
    axis([-2 2 -2.5 1.5]);
  
  
  fprintf('This will be the EKF solution. Press enter.\n');
  pause;
  
  %end %% <--- Uncomment to disable

