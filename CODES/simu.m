%The code incorporate ideas presented in em.m from Higham paper
%randn('state',100)

clear
delta = [1 1 0]';
a = [-1 0 1]';
a0 = 0;
sigma = [0.200 0 0;0.0375 0.1452 0;0.0250 0.0039 0.0967];
T = 1;
N = 2^13;
gamma = 0.6;

dt = T/N;
Yzero = transpose([11.10 12.00 11.00]);
dW = sqrt(dt)*randn(3,N);

omega = sigma*sigma';
kappa = -1 * delta' * a;
A = diag(a);
theta = trace(A*omega)/(2*delta'*a);


R = 4;
Dt = R*dt;
L = N/R;

Yem = zeros(3,L);
alphaem = zeros(1,L);
aem = zeros(1,L);
Ytemp = Yzero;
Yem(:,1) = Ytemp;
alphatemp = a0 + a'*log(Ytemp);
abt = alphatemp;
alphaem(1) = alphatemp;
aem(1) = abt;
piem = zeros(3,L);
piem(:,1) = (1/gamma)*(((omega^(-1))*delta)*alphatemp - a*((delta'*omega*delta)*(L*Dt*alphatemp - 0.25*trace(A*omega)*(L*Dt)^2)));
wealthem = zeros(1,L);
wealthtemp = [0 0 0]';
wealthem(1) =0 ;

for j = 1:L
    Winc = sum(dW(:,R*(j-1)+1:R*j),2);
    %pip = pitemp;
    %for k = 1:3
        %Ytemp(k) = - Ytemp(k)*((Dt*alphatemp)*(delta(k)) + sigma(k,:) * Winc -1)^-1;
    %end
    Ytemp = Ytemp + (Dt*alphatemp)*(delta.*Ytemp) + sigma * Winc.*Ytemp;
    %Ytemp = - Ytemp*((Dt*alphatemp)*(delta) + sigma * Winc-1)^-1;
    pitemp = (1/gamma)*(((omega^(-1))*delta)*alphatemp - a*((delta'*omega*delta)*((L-j)*Dt*alphatemp - 0.25*trace(A*omega)*((L-j)*Dt)^2)));
    wealthtemp = wealthtemp + pitemp'*delta*alphatemp*Dt + pitemp'*sigma * Winc;
    %wealthtemp = wealthtemp + sum(pitemp-pip,1);
    abt = abt + kappa*(theta-abt)*Dt+a'*sigma*Winc;
    alphatemp = a0 + a'*log(Ytemp);
    Yem(:,j+1) = Ytemp;
    alphaem(j+1) = alphatemp;
    aem(j+1) = abt;
    piem(:,j+1) = pitemp;
    wealthem(j+1) = sum(wealthtemp);
end

%wealth = sum(abs(piem),1);

figure

subplot(2,2,1)
for l = 1:3
    plot([0:Dt:T],[Yem(l,:)]),hold on
end
hold off
xlabel('Time','FontSize',12)
ylabel('Asset Prices','FontSize',12,'Rotation',90)
legend('Y_{t}^{1}','Y_{t}^{2}','Y_{t}^{3}')

subplot(2,2,2)
plot([0:Dt:T],alphaem)
xlabel('Time','FontSize',12)
ylabel('Co-Intergration Factor','FontSize',12,'Rotation',90)

subplot(2,2,3)
for k = 1:3
    plot([0:Dt:T],piem(k,:)./Yem(k,:)),hold on
end
hold off
legend('\pi^{1}/Y_{t}^{1}','\pi^{2}/Y_{t}^{2}','\pi^{3}/Y_{t}^{3}')
xlabel('Time','FontSize',12)
ylabel('Position','FontSize',12,'Rotation',90)

subplot(2,2,4)
plot([0:Dt:T],wealthem)
xlabel('Time','FontSize',12)
ylabel('Wealth1','FontSize',12,'Rotation',90)

figure

plot([0:Dt:T],aem)
xlabel('Time','FontSize',12)
ylabel('aem','FontSize',12,'Rotation',90)

[omg,loomg,ek,uk,uuk,mm] = VARMOD(Yem,Dt);

if all(imag(diag(ek)) == 0)
    disp('real');
    for k = 1:length(uk(1,:))
        if all(abs(ek(k,k)) >= abs(diag(ek))) == 1
            disp(k);
            sluk = uuk(k,:);
        end
    end
else
    disp('imag');
    for l = 1:length(uk(1,:))
        if all((imag(uuk(l,:)))<1.0e-7) == 1
            disp(l);
            sluk = uuk(l,:);
        end
    end
end

varstd = 0;

for stdnb1 = 1:length(sluk)
    for stdnb2 = 1:length(sluk)
        varstd = varstd + sluk(stdnb1)*sluk(stdnb2)*loomg(stdnb1,stdnb2);
    end
end
std = sqrt(varstd);

EXPATRAD(sluk,Yem,mm,std,L,Dt,T)