function DDD = Des(Yem)
Dsum = 0;
for i = 1:length(Yem(1,:))
    if i == 1
        Dsum = Dsum + Yem(:,i).*Yem(:,i) - 2 * Yem(:,i).* mean(Yem(:,1:i),2)+zeros(length(Yem(:,1)),length(Yem(:,1)))+mean(Yem(:,1:i),2).*mean(Yem(:,1:i),2);
    else
        Dsum = Dsum + Yem(:,i).*Yem(:,i) - 2 * Yem(:,i).* mean(Yem(:,1:i),2)+cov((Yem(:,1:i))')+mean(Yem(:,1:i),2).*mean(Yem(:,1:i),2);
    end
    
end
DDD = Dsum/(length(Yem)+1);
end