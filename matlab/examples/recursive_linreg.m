%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Recursive linear regression demonstration from Chapter 3 of the book
%
% Simo Sarkka and Lennart Svensson (2023), Bayesian Filtering and Smoothing,
% 2nd ed. Cambridge University Press. 
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
    sd = 0.1;
    t = (0:dt:1);
    x = 1 + 0.5*t;
    y = x + sd*randn(size(x));
    
    h = plot(t,y,'.',t,x);
    
    axis([0 1 0.5 2]);
    
    set(h,'Markersize',7);
    set(h,'LineWidth',4);
    set(h(1),'Color',[0.0 0.0 0.0]);
    set(h(2),'Color',[0.7 0.7 0.7]);
        
    h = legend('Measurement','True Signal');
    xlabel('{\it t}');
    ylabel('{\it y}');
      
%%
% Batch linear regression
%
  m0 = [0;0];
  P0 = 1*eye(2);
  n  = length(y);
  H  = [ones(length(t),1) t'];
  mb = inv(inv(P0) + 1/sd^2*H'*H)*(1/sd^2*H'*y'+inv(P0)*m0)
  Pb = inv(inv(P0) + 1/sd^2*H'*H);
  
  h = plot(t,y,'.',t,x,t,mb(1)+mb(2)*t,'-');

  axis([0 1 0.5 2]);
  
  set(h,'Markersize',7);
  set(h(2),'LineWidth',4);
  set(h(3),'LineWidth',1.5);
  set(h(1),'Color',[0.0 0.0 0.0]);
  set(h(2),'Color',[0.7 0.7 0.7]);
  set(h(3),'Color',[0.0 0.0 0.0]);
  
  h = legend('Measurement','True Signal','Estimate');
  xlabel('{\it t}');
  ylabel('{\it y}');
  
%%
% Recursive solution
%
  m = m0;
  P = P0;
  MM = zeros(size(m0,1),length(y));
  PP = zeros(size(P0,1),size(P0,1),length(y));
  for k=1:length(y)
      H = [1 t(k)];
      S = H*P*H'+sd^2;
      K = P*H'/S;
      m = m + K*(y(k)-H*m);
      P = P - K*S*K';
      
      MM(:,k) = m;
      PP(:,:,k) = P;
  end
  % Note that this last estimate is exatly the same as mb:
  m

  h = plot(t,MM(1,:),'b-',[0 1],[mb(1) mb(1)],'b--',...
           t,MM(2,:),'r-',[0 1],[mb(2) mb(2)],'r--');
  
  set(h,'Markersize',10);
  set(h,'LineWidth',2);
  set(h(1),'Color',[0.0 0.0 0.0]);
  set(h(2),'Color',[0.0 0.0 0.0]);
  set(h(3),'Color',[0.5 0.5 0.5]);
  set(h(4),'Color',[0.5 0.5 0.5]);

  
  h = legend('Recursive E[ {\it\theta}_1 ]','Batch E[ {\it\theta}_1 ]',...
         'Recursive E[ {\it\theta}_2 ]','Batch E[ {\it\theta}_2 ]','Location','SE');
  
  xlabel('{\it t}');
  ylabel('{\it y}');
  grid on;

%%
  h = semilogy(t,squeeze(PP(1,1,:)),'b-',[0 1],[Pb(1,1) Pb(1,1)],'b--',...
               t,squeeze(PP(2,2,:)),'r-',[0 1],[Pb(2,2) Pb(2,2)],'r--');

  set(h,'Markersize',10);
  set(h,'LineWidth',2);
  set(h(1:2),'Color',[0.0 0.0 0.0]);
  set(h(3:4),'Color',[0.5 0.5 0.5]);
  
  h = legend('Recursive Var[ {\it\theta}_1 ]','Batch Var[ {\it\theta}_1 ]',...
         'Recursive Var[ {\it\theta}_2 ]','Batch Var[ {\it\theta}_2 ]');     
  
  xlabel('{\it t}');
  ylabel('{\it y}');
  grid on;
  
