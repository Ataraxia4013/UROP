function PLOT(Yem,Dt,T,alphaem,piem,wealthem,sluk,std,mm,Wtt,bbeta,L)
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
end