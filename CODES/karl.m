function [mu,esvarm] = karl(A,B,C,D,Yem)
mu = Yem(:,1);
muu = [];
esvarm = zeros(length(Yem(:,1)),length(Yem(:,1)));
for l = 1:length(Yem(1,:))-1
    mum = A + B*mu;
    esvarmm = B*esvarm*B' + C;
    kgain = esvarmm*((D+esvarmm)^(-1));
    mu = mum+kgain*(Yem(:,l+1)-mum);
    esvarm = esvarmm - kgain*esvarmm;
    muu = [muu mu];
end

end