
function simufunct(Yzero,delta,a,a0,sigma,N,T,gamma)
%The code incorporate ideas presented in em.m from Higham paper

%Input Variables
%Yzero : Initial Asset Price 
%sigma : volitility
%N     : number of time steps
%T     : Trading Horizon Length
%gamma : Index for utility function (exponential in this case)

%dW: Brownian increments
randn('state',100)
dt = T/N;
dW = sqrt(dt)*randn(3,N);


%Terms involved in the ornstein-urenbeck process (mean reversion)
omega = sigma*sigma';
kappa = -1 * delta' * a;
A     = diag(a);
theta = trace(A*omega)/(2*delta'*a);
                          


%R :Step coefficient for Euler-Maruyama Method
%L :Steps in Euler-Maruyama Method
%Dt:Upgraded time step
R  = 1;
Dt = R*dt;
L  = N/R;


%Yem        : Set of Asset Values in one Trading Horizon
%Ytemp      : Stepwise Asset Value
%alphaem    : Set of Co-integration Factors in one Trading Horizon
%alphatemp  : Stepwise Co-integration Factor
%piem       : Set of Positions in one Trading Horizon
%wealthem   : Set of Wealth of Trader in one Trading Horizon
%wealthtemp : Stepwise Wealth of Trader
Yem        = zeros(3,L);
Ytemp      = Yzero;
alphaem    = zeros(1,L);
alphatemp  = a0 + a'*log(Ytemp);
piem       = zeros(3,L);
wealthem   = zeros(1,L);
wealthtemp = [0 0 0]';

for j = 1:L
    Winc = sum(dW(:,R*(j-1)+1:R*j),2);
    alphatemp = alphatemp + kappa*(theta - alphatemp)*Dt + a'*sigma*Winc;
    Ytemp = Ytemp + (Dt*alphatemp)*(delta.*Ytemp) + sigma * Winc;
    pitemp = (1/gamma)*((inv(omega)*delta)*alphatemp + (delta'*omega*delta)*((L-j)*Dt*alphatemp/2 + 0.25*trace(A*omega)*((L-j)*Dt)^2)*a);
    wealthtemp = wealthtemp + pitemp.*((Dt*alphatemp)*(delta.*Ytemp) + sigma * Winc)./Ytemp;
    Yem(:,j) = Ytemp;
    alphaem(j) = alphatemp;
    piem(:,j) = pitemp;
    wealthem(j) = sum(wealthtemp);
end

%wealth = sum(abs(piem),1);

figure

subplot(2,2,1)
for l = 1:3
    plot([0:Dt:T],[Yzero(l),Yem(l,:)]),hold on
end
hold off
xlabel('Time','FontSize',12)
ylabel('Asset Prices','FontSize',12,'Rotation',90)
legend('Y_{t}^{1}','Y_{t}^{2}','Y_{t}^{3}')

subplot(2,2,2)
plot([0:Dt:T],[0,alphaem])
xlabel('Time','FontSize',12)
ylabel('Co-Intergration Factor','FontSize',12,'Rotation',90)

subplot(2,2,3)
for k = 1:3
    plot([0:Dt:T],[0,piem(k,:)./Yem(k,:)]),hold on
end
hold off
legend('\pi^{1}/Y_{t}^{1}','\pi^{2}/Y_{t}^{2}','\pi^{3}/Y_{t}^{3}')
xlabel('Time','FontSize',12)
ylabel('Position','FontSize',12,'Rotation',90)

subplot(2,2,4)
plot([0:Dt:T],[0,wealthem])
xlabel('Time','FontSize',12)
ylabel('Wealth','FontSize',12,'Rotation',90)

end



