randn('state',100)
T = 1;
N = 2^8;
dt = T/N;
Yzero = transpose([11.10 12.00 11.00]);
dW = sqrt(dt)*randn(3,N);

delta = [1 1 0]';
a = [-1 0 1]';
a0 = 0;
sigma = [0.200 0 0;0.0375 0.1452 0;0.0250 0.0039 0.0967];
omega = sigma*sigma';
kappa = -1 * delta' * a;
A = diag(a);
theta = trace(A*omega)/(2*delta'*a);

R = 1;
Dt = R*dt;
L = N/R;

Yem = zeros(3,L);
alphaem = zeros(1,L);
Ytemp = Yzero;
alphatemp = a0 + a'*log(Ytemp);

for j = 1:L
    Winc = sum(dW(:,R*(j-1)+1:R*j),2);
    alphatemp = alphatemp + kappa*(theta - alphatemp)*Dt + a'*sigma*dW(:,j);
    Ytemp = Ytemp + (Dt*alphatemp)*(delta.*Ytemp) + sigma * Winc;
    Yem(:,j) = Ytemp;
    alphaem(j) = alphatemp;
end

figure
subplot(2,1,1)
for l = 1:3
    plot([0:Dt:T],[Yzero(l),Yem(l,:)]),hold on
end
hold off
xlabel('Time','FontSize',12)
ylabel('Asset Prices','FontSize',12,'Rotation',90)
legend('Yt1','Yt2','Yt3')
subplot(2,1,2)
plot([0:Dt:T],[0,alphaem])
xlabel('Time','FontSize',12)
ylabel('Co-Intergration Factor','FontSize',12,'Rotation',90)