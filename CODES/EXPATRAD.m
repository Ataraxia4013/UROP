function EXPATRAD(sluk,Yem,mm,std,L,Dt,T)
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
    Wt = cash + beta * ycoi;
    if ycoi >= sluk*mm+std && beta == 0
        beta = 1;
        ytempl = ytemp;
        Wt = cash + beta * ycoi;       
    elseif ycoi <=sluk*mm-std && beta == 0
        beta = -1;
        ytemps = ytemp;
        Wt = cash + beta * ycoi;
    elseif ycoi <= sluk*mm+ 0.1*std && beta == 1
        cash = cash+ beta * sluk*(ytempl - ytemp);
        beta =0 ;
    elseif ycoi>= sluk*mm-0.1*std && beta == -1
        cash = cash + beta * sluk*(ytemps - ytemp);
        beta =0;
    else
    end
    
    Wtt(lll+1) = Wt;
    bbeta(lll+1) = beta;
end

figure
plot([0:Dt:T],Wtt);
xlabel('wtt','FontSize',12)
figure
yyaxis right
plot([0:Dt:T],bbeta),hold on
yyaxis left
plot([0:Dt:T],sluk*Yem),hold on
xlabel('slukyem','FontSize',12)
sdp = ones(1,L+1)*(sluk*mm+std);
sdm = ones(1,L+1)*(sluk*mm-std);
sdpp = ones(1,L+1)*(sluk*mm+0.1*std);
sdmm = ones(1,L+1)*(sluk*mm-0.1*std);
plot([0:Dt:T],sdp,'-'),hold on
plot([0:Dt:T],sdm,'-'),hold on
plot([0:Dt:T],sdpp,'-'),hold on
plot([0:Dt:T],sdmm,'-'),hold off
end