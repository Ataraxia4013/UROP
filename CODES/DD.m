Dsum = 0;

for i = 1:length(Yem(1,:))
    if i == 1
        Dsum = Dsum + Yem(:,i).*Yem(:,i) - 2 * Yem(:,i).* mean(Yem(:,1:i),2)+zeros(length(Yem(:,1)),length(Yem(:,1)))+mean(Yem(:,1:i),2).*mean(Yem(:,1:i),2);
    else
        Dsum = Dsum + Yem(:,i).*Yem(:,i) - 2 * Yem(:,i).* mean(Yem(:,1:i),2)+cov((Yem(:,1:i))')+mean(Yem(:,1:i),2).*mean(Yem(:,1:i),2);
    end
    
end

DDD = Dsum/(length(Yem)+1);
[omg,loomg,ek,uk,uuk,mm,sluk,std,louuk,loek,con,BB,omg2] = VARMOD(Yem,Dt);

mu = Yem(:,1);
muu = [];
esvarm = zeros(3,3);
for l = 1:length(Yem(1,:))-1
    [mu,esvarm,kgain] = karl(con,BB,omg2,DDD,mu,esvarm,Yem(:,l+1));
    muu = [muu mu];
end
