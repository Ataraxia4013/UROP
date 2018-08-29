function [omg,loomg,ek,uk,uuk,mm,sluk,std,louuk,loek,con,BB,omg2 ] = VARMOD(Yem,Dt)
%lomo  : var model initiation 
%eslomo: log var model estimation 
%loBB  : B of Yt = A + B*Yt-1 + noise
%loomg2: estimated variance-covariance matrix from the model
%losisi: cholesky decomposition of estimated variance loomg2
%lobinv: Inverse of loBB
%losi  : adjusted standard deviation 
%loomg : adjusted final variance-covariance matrix of the data
%lolol : adjusted standard deviation of the data
%lokpii: estimated kappa value of mean reverting
%louk  : eigenvector of lokpii
%loek  : eigenvalues of lokpii
%louuk : inverse of louk
%lovar : adjusted variance of data


lomo = varm(length(Yem(:,1)),1);
eslomo = estimate(lomo,log(Yem'));
loBB = eslomo.AR{1,1};
loomg2 = eslomo.Covariance;
losisi = chol(loomg2);
losi = losisi';
loomg = losi*losi'/Dt;

lokpii = (eye(length(Yem(:,1)))-loBB)/Dt;
[louk,loek] = eig(lokpii);
louuk = louk^-1;
%locc = lobinv * eslomo.Constant/ Dt;

%AAPL,INTC,QCOM


%following variables naming convention follow previous ones without the 'lo'
mo = varm(length(Yem(:,1)),1);
esmo = estimate(mo,Yem');
BB = esmo.AR{1,1};
omg2 = esmo.Covariance;
sisi = chol(omg2);
si = sisi';
omg = si*si'/Dt;
con = esmo.Constant;

kpii = (eye(length(Yem(:,1)))-BB)/Dt;
[uk,ek] = eig(kpii);
uuk = uk^-1;

%mm: estimated mean reversion mean
mm = (eye(length(Yem(:,1)))-BB)^(-1) * con;


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

end
