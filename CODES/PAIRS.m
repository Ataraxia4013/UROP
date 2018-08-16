function [beta,ytempl] = PAIRS(sluk,ytemp,mm,std,cash,ytempl,beta)
ycoi = sluk*ytemp;
if ycoi >= sluk*mm+std & beta == 0
    beta = -1;
    ytempl = ytemp;     
elseif ycoi <=sluk*mm-std & beta == 0
    beta = 1;
    ytempl = ytemp;
elseif ycoi <= sluk*mm+ 0.1*std & beta == -1
    cash = cash+ beta * sluk*(ytemp - ytempl);
    ytempl = 0;
    beta =0 ;
elseif ycoi>= sluk*mm-0.1*std & beta == 1
    cash = cash + beta * sluk*(ytemp - ytempl);
    ytempl = 0;
    beta =0;
end
end