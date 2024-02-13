%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulate pendulum data for the examples in the book
%
% Simo Sarkka and Lennart Svensson (2023), Bayesian Filtering and Smoothing,
% 2nd ed., Cambridge University Press.
% 
% See LICENSE provided with the software.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%
    % Simulate simple pendulum. Note that the system easily
    % diverges, but it should not matter.
    %
    
    % Remember to comment these out when benchmarking!
    fprintf('!! Using fixed random stream !!\n');
    rng(1,'twister');
    % <--- 
    
    DT = 0.01;
    g  = 9.81;
    Q  = 0.01*[DT^3/3 DT^2/2; DT^2/2 DT];
    R  = 0.1;
%    m0 = [1.6;0]; % Slightly off
%    P0 = 0.1*eye(2);
    m0 = [0;0];
    P0 = eye(2);
    
    steps = 500;
    
    QL = chol(Q,'lower');
    
    T = zeros(1,steps);
    X = zeros(2,steps);
    Y = zeros(1,steps);
    t = 0;
    x = [1.5;0];
    for k=1:steps
        x = [x(1)+x(2)*DT;
            x(2)-g*sin(x(1))*DT];
        w = QL * randn(2,1);
        x = x + w;
        y = sin(x(1)) + sqrt(R)*randn;
        t = t + DT;
        T(k) = t;
        X(:,k) = x;
        Y(:,k) = y;
    end
    
    % Plot the data
    clf;
    plot(T,Y,'g.',T,X(1,:),'r-');

    % Store the data in JSON
    jstruct = struct('T',T,'X',X,'Y',Y);
    json = jsonencode(jstruct,'PrettyPrint',true);
    filename = 'pendulum.json';
    fid = fopen(filename, 'w');
    fwrite(fid, json);
    fclose(fid);


