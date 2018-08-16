%The code incorporate ideas presented in em.m from Higham paper
randn('state',100)

clear
close all

delta = [1 1 0]';
a = [-1 0 1]';
a0 = 0;
sigma = [0.200 0 0;0.0375 0.1452 0;0.0250 0.0039 0.0967];
T = 1;
N = 2^13;
gamma = 0.1;

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

pitemp = (1/gamma)*(((omega^(-1))*delta)*alphatemp - a*((delta'*(omega^(-1))*delta)*(L*Dt*alphatemp - 0.25*trace(A*omega)*(L*Dt)^2)));
piem(:,1) = pitemp;
wealthem = zeros(1,L);
wealthtemp = [0 0 0]';
wealthem(1) =0 ;

possem = zeros(3,L);
posstemp = pitemp./Ytemp;
possem(:,1) = posstemp;

for j = 1:L
    Winc = sum(dW(:,R*(j-1)+1:R*j),2);
    %pip = pitemp;
    %for k = 1:3
        %Ytemp(k) = - Ytemp(k)*((Dt*alphatemp)*(delta(k)) + sigma(k,:) * Winc -1)^-1;
    %end
    Ytemp = Ytemp + (Dt*alphatemp)*(delta.*Ytemp) + sigma * Winc.*Ytemp;
    alphatemp = a0 + a'*log(Ytemp);
    %Ytemp = - Ytemp*((Dt*alphatemp)*(delta) + sigma * Winc-1)^-1;
    pitemp = (1/gamma)*(((omega^(-1))*delta)*alphatemp - a*((delta'*(omega^(-1))*delta)*((L-j)*Dt*alphatemp - 0.25*trace(A*omega)*((L-j)*Dt)^2)));
    posstemp = pitemp./Ytemp;
    %wealthtemp = wealthtemp + pitemp'*delta*alphatemp*Dt + pitemp'*sigma * Winc;
    wealthtemp = wealthtemp + sum(pitemp'*((Ytemp-Yem(:,j))./Ytemp));
    abt = abt + kappa*(theta-abt)*Dt+a'*sigma*Winc;
    disp(Yem(:,j));
    disp(Ytemp);
    Yem(:,j+1) = Ytemp;
    alphaem(j+1) = alphatemp;
    aem(j+1) = abt;
    piem(:,j+1) = pitemp;
    %wealthem(j+1) = sum(wealthtemp);
    wealthem(j+1) = sum(wealthtemp);
    possem(:,j+1) = posstemp;
end

%wealth = sum(abs(piem),1);

figure
plot([0:Dt:T],possem);


lomo = varm(3,1);
eslomo = estimate(lomo,log(Yem'));
loBB = eslomo.AR{1,1};
loomg2 = eslomo.Covariance;
losisi = chol(loomg2);
lobinv = (loBB^(-1));
losi = lobinv*losisi';
loomg = losi*losi'/Dt;
lolol = chol(loomg);

lokpii = lobinv*(eye(3)-loBB)/Dt;
[louk,loek] = eig(lokpii);
louuk = louk^-1;
locc = lobinv * eslomo.Constant/ Dt;
lovar = diag(loomg);

%AAPL,INTC,QCOM


%following variables naming convention follow previous ones without the 'lo'
mo = varm(3,1);
esmo = estimate(mo,Yem');
BB = esmo.AR{1,1};
omg2 = esmo.Covariance;
sisi = chol(omg2);
binv = (BB^(-1));
si = binv*sisi';
omg = si*si'/Dt;
lol = chol(omg);

kpii = binv*(eye(3)-BB)/Dt;
[uk,ek] = eig(kpii);
uuk = uk^-1;

%mm: estimated mean reversion mean
mm = (eye(3)-BB)^(-1) * esmo.Constant;

var = diag(omg);


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

if all(imag(sluk) < 1.0e-8) == 1
    sluk = real(sluk);
end

varstd = 0;

for stdnb1 = 1:length(sluk)
    for stdnb2 = 1:length(sluk)
        varstd = varstd + sluk(stdnb1)*sluk(stdnb2)*loomg(stdnb1,stdnb2);
    end
end
std = sqrt(varstd);

beta = 0;
bbeta = zeros(1,L);
bbeta(1) = 0;
Wtt = zeros(1,L);
Wtt(1) = 0;
Wt = 0;
cash = 0;
for lll = 1:L
    ytemp = Yem(:,lll);
    ycoi = sluk*ytemp;
    if lll == 1
        ytempl = ytemp;
    end
    if ycoi >= sluk*mm+std && beta == 0
        beta = -1;
        ytempl = ytemp;     
    elseif ycoi <=sluk*mm-std && beta == 0
        beta = 1;
        ytempl = ytemp;
    elseif ycoi <= sluk*mm+ 0.1*std && beta == -1
        cash = cash+ beta * sluk*(ytemp - ytempl);
        ytempl = 0;
        beta =0 ;
    elseif ycoi>= sluk*mm-0.1*std && beta == 1
        cash = cash + beta * sluk*(ytemp - ytempl);
        ytempl = 0;
        beta =0;
    end
    Wt = cash + beta * sluk*(ytemp - ytempl);
    Wtt(lll+1) = Wt;
    bbeta(lll+1) = beta;
end

figure

subplot(2,3,1)
for l = 1:3
    plot([0:Dt:T],[Yem(l,:)]),hold on
end
hold off
xlabel('Time','FontSize',12);
ylabel('Asset Prices','FontSize',12,'Rotation',90);
legend('Y_{t}^{1}','Y_{t}^{2}','Y_{t}^{3}');

subplot(2,3,2)
plot([0:Dt:T],alphaem);
xlabel('Time','FontSize',12);
ylabel('Co-Intergration Factor \alpha','FontSize',12,'Rotation',90);

subplot(2,3,3)
for k = 1:3
    plot([0:Dt:T],piem(k,:)./Yem(k,:)),hold on
end
hold off
legend('\pi^{1}/Y_{t}^{1}','\pi^{2}/Y_{t}^{2}','\pi^{3}/Y_{t}^{3}');
xlabel('Time','FontSize',12);
ylabel('Position','FontSize',12,'Rotation',90);

subplot(2,3,4)
plot([0:Dt:T],wealthem);
xlabel('Time','FontSize',12);
ylabel('Wealth1','FontSize',12,'Rotation',90);

subplot(2,3,5)
plot([0:Dt:T],Wtt);
xlabel('Time','FontSize',12)
ylabel('Wealth of pairs trading','FontSize',12,'Rotation',90);

subplot(2,3,6)
yyaxis right
plot([0:Dt:T],bbeta),hold on
ylabel('Position','FontSize',12,'Rotation',90);
yyaxis left
plot([0:Dt:T],sluk*Yem),hold on
xlabel('Time','FontSize',12)
ylabel('Cointegration factor of pairs trading','FontSize',12,'Rotation',90);
sdp = ones(1,L+1)*(sluk*mm+std);
sdm = ones(1,L+1)*(sluk*mm-std);
sdpp = ones(1,L+1)*(sluk*mm+0.1*std);
sdmm = ones(1,L+1)*(sluk*mm-0.1*std);
plot([0:Dt:T],sdp,'-'),hold on
plot([0:Dt:T],sdm,'-'),hold on
plot([0:Dt:T],sdpp,'-'),hold on
plot([0:Dt:T],sdmm,'-'),hold off