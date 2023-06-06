%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulate cluttered pendulum data for the examples in the book
%
% Simo Sarkka and Lennart Svensson (2023), Bayesian Filtering and Smoothing,
% 2nd ed., Cambridge University Press.
% 
% See LICENSE provided with the software.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Simulate simple pendulum. Note that the system easily diverges, but
% it should not matter. Also draw some clutter measurements.
    
    % Remember to comment these out when benchmarking!
    fprintf('!! Using fixed random stream !!\n');
    rng(1,'twister');
    % <---
    
    DT = 0.01;
    g  = 9.81;
    Q  = 0.01*[DT^3/3 DT^2/2; DT^2/2 DT];
    R  = 0.1;
    m0 = [0;0];
    P0 = eye(2);
        
    steps = 500;

    QL = chol(Q,'lower');
   
    T = [];
    X = [];
    Y = [];
    t = 0;
    x = [1.5;0];
    for k=1:steps
        x = [x(1)+x(2)*DT;
             x(2)-g*sin(x(1))*DT];
        w = QL * randn(2,1);
        x = x + w;
        y = sin(x(1)) + sqrt(R)*randn;
        t = t + DT;
        X = [X x];
        Y = [Y y];
        T = [T t];
  end
  
  cp = 0.5;
  C = rand(size(Y)) < cp;
  ind = find(C);
  Y(ind) = 4*rand(size(ind))-2;
  
  % Plot the data    
  plot(T,Y,'g.',T,X(1,:),'r-');

