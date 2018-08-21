function [Wtt,bbeta] = EXPATRAD(sluk,Yem,mm,std)
L = length(Yem);
beta = 0;
bbeta = [];
Wtt = [];
cash = 0;
flag = 0;
for lll = 1:L
    ytemp = Yem(:,lll);
    ycoi = sluk*ytemp;
    if lll == 1
        ytempl = ytemp;
    end
    if ycoi >= sluk*mm+std && beta == 0 && flag == 0
        beta = -1;
        ytempl = ytemp;     
    elseif ycoi <=sluk*mm-std && beta == 0 && flag == 0
        beta = 1;
        ytempl = ytemp;
    elseif ycoi <= sluk*mm+ 0.1*std && beta == -1 && flag == 0
        cash = cash+ beta * sluk*(ytemp - ytempl);
        ytempl = 0;
        beta = 0 ;
    elseif ycoi>= sluk*mm-0.1*std && beta == 1 && flag == 0
        cash = cash + beta * sluk*(ytemp - ytempl);
        ytempl = 0;
        beta = 0;
    elseif ycoi <= sluk*mm-2*std && beta == 1 && flag == 0
        cash = cash + beta * sluk*(ytemp - ytempl);
        ytempl = 0;
        beta = 0;
        flag = 1;
    elseif ycoi >= sluk*mm+2*std && beta == -1 && flag == 0
        cash = cash + beta * sluk*(ytemp - ytempl);
        ytempl = 0;
        beta = 0;
        flag = -1;
    elseif ycoi <= sluk*mm+1.1*std && beta == 0 && flag == -1
        flag = 0;
    elseif ycoi >= sluk*mm-1.1*std && beta == 0 && flag == 1
        flag = 0;
    end
    Wt = cash + beta * sluk*(ytemp - ytempl);
    Wtt = [Wtt Wt];
    bbeta = [bbeta beta];
end
end