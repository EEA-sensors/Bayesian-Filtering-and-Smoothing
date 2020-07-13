%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Estimate pendulum state with SLF and SLRTS as in Examples 5.2
% and 9.2 of the book
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

%%
% Simulate data
%
    pendulum_sim;

%%
% Filter
%

    m = m0;
    P = P0;
    MM = zeros(size(m,1),length(Y));
    PP = zeros(size(P,1),size(P,2),length(Y));
    for k=1:length(Y)
        mf = [m(1)+m(2)*DT; m(2)-g*sin(m(1))*DT*exp(-P(1,1)/2)];
        mfdx = [...
            P(1,1)+DT*P(1,2) P(1,2)+DT*P(2,2);
            P(1,2)-g*DT*cos(m(1))*P(1,1)*exp(-P(1,1)/2) P(2,2)-g*DT*cos(m(1))*P(1,2)*exp(-P(1,1)/2)];
        m = mf;
        P = (mfdx/P)*mfdx'+Q;        
    
        mh = sin(m(1))*exp(-P(1,1)/2);
        mhdx = [cos(m(1))*P(1,1)*exp(-P(1,1)/2) cos(m(1))*P(1,2)*exp(-P(1,1)/2)];
    
        S = (mhdx/P)*mhdx'+R;
        K = mhdx'/S;
        m = m + K*(Y(k) - mh);
        P = P - K*S*K';
    
        MM(:,k) = m;
        PP(:,:,k) = P;
    end

    h = plot(T,Y,'k.',T,X(1,:),'r-',T,MM(1,:),'b--');
    set(h,'Linewidth',5);
    title('SLF estimate');
    legend('Measurements','True','Estimate');

    rmse_slf = sqrt(mean((X(1,:)-MM(1,:)).^2))
  
%%
% Smoother
%

    ms = m;
    Ps = P;
    MS = zeros(size(m,1),length(Y));
    PS = zeros(size(P,1),size(P,2),length(Y));
    MMS(:,k) = m;
    PPS(:,:,k) = P;
    for k=size(MM,2)-1:-1:1
        m = MM(:,k);
        P = PP(:,:,k);
        mf = [m(1)+m(2)*DT; m(2)-g*sin(m(1))*DT*exp(-P(1,1)/2)];
        mfdx = [...
            P(1,1)+DT*P(1,2) P(1,2)+DT*P(2,2);
            P(1,2)-g*DT*cos(m(1))*P(1,1)*exp(-P(1,1)/2) P(2,2)-g*DT*cos(m(1))*P(1,2)*exp(-P(1,1)/2)];
    
        mp = mf;
        Pp = (mfdx/P)*mfdx'+Q;
        Ck = mfdx'/Pp;
        ms = m + Ck*(ms - mp);
        Ps = P + Ck*(Ps - Pp)*Ck';
        MMS(:,k) = ms;
        PPS(:,:,k) = Ps;
    end
  
    h = plot(T,Y,'k.',T,X(1,:),'r-',...
        T,MM(1,:),'b--',T,MMS(1,:),'g--');
    set(h,'Linewidth',5);
    title('SLF and SLRTS estimates');
    legend('Measurements','True','SLF','SLRTS');
    
    rmse_slrts = sqrt(mean((X(1,:)-MMS(1,:)).^2))
    
%%
% Plot the filtering result
%

    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,MM(1,:),'r');
    legend('True Angle','Measurements','SLF Estimate');
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
    

%%
% Plot the smoothing result
%
    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,MMS(1,:),'r');
    legend('True Angle','Measurements','SLRTSS Estimate');
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
    
    