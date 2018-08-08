function [Wtt,bbeta] = EXPATRAD(sluk,Yem,mm,std,L)
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
end