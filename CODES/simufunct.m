
function [Yem,piem,alphaem,wealthem,Dt,L] = simufunct(Yzero,delta,a,a0,sigma,N,T,gamma)
%The code incorporate ideas presented in em.m from Higham paper

%Input Variables
%Yzero : Initial Asset Price 
%sigma : volitility
%N     : number of time steps
%T     : Trading Horizon Length
%gamma : Index for utility function (exponential in this case)

%dW: Brownian increments
%randn('state',100)
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
R  = 4;
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
Yem(:,1)   = Yzero;

alphaem    = zeros(1,L);
alphatemp  = a0 + a'*log(Ytemp);
alphaem(1) = alphatemp;

piem       = zeros(3,L);
piem(:,1)  = (1/gamma)*(((omega^(-1))*delta)*alphatemp - a*((delta'*omega*delta)*((L)*Dt*alphatemp - 0.25*trace(A*omega)*((L)*Dt)^2)));

wealthem   = zeros(1,L);
wealthtemp = [0 0 0]';
wealthem(1)= 0;

for j = 1:L
    Winc       = sum(dW(:,R*(j-1)+1:R*j),2);
    
    Ytemp        = Ytemp + (Dt*alphatemp)*(delta.*Ytemp) + sigma * Winc.*Ytemp;
    pitemp       = (1/gamma)*(((omega^(-1))*delta)*alphatemp - a*((delta'*omega*delta)*((L-j)*Dt*alphatemp - 0.25*trace(A*omega)*((L-j)*Dt)^2)));
    wealthtemp   = wealthtemp + pitemp'*delta*alphatemp*Dt + pitemp'*sigma * Winc;
    alphatemp    = a0 + a'*log(Ytemp);
    Yem(:,j+1)   = Ytemp;
    alphaem(j+1) = alphatemp;
    piem(:,j+1)  = pitemp;
    wealthem(j+1)= sum(wealthtemp);
end
end



