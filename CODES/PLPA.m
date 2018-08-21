function PLPA(sluk,mm,std,Wtt,bbeta,Yem)
L = length(Yem);
figure
subplot(2,1,1)
yyaxis right
plot(1:L,bbeta),hold on
ylabel('Position','FontSize',12,'Rotation',90);
yyaxis left
plot(1:L,sluk*Yem),hold on
xlabel('Time','FontSize',12)
ylabel('Cointegration factor of pairs trading','FontSize',12,'Rotation',90);
sdp = ones(1,L)*(sluk*mm+std);
sdm = ones(1,L)*(sluk*mm-std);
sdpp = ones(1,L)*(sluk*mm+0.1*std);
sdmm = ones(1,L)*(sluk*mm-0.1*std);
sdsll = ones(1,L)*(sluk*mm-2*std);
sdslu = ones(1,L)*(sluk*mm+2*std);
plot(1:L,sdp,'-'),hold on
plot(1:L,sdm,'-'),hold on
plot(1:L,sdpp,'-'),hold on
plot(1:L,sdsll,'r-'),hold on
plot(1:L,sdslu,'r-'),hold on
plot(1:L,sdmm,'-'),hold off


subplot(2,1,2)
plot(1:L,Wtt);
xlabel('Time','FontSize',12)
ylabel('Wealth of pairs trading','FontSize',12,'Rotation',90);
end