function [omg,loomg,ek,uk,uuk] = VARMOD(Yem,Dt)
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


lomo = varm(3,1);
eslomo = estimate(lomo,log(Yem'));
loBB = eslomo.AR{1,1};
loomg2 = eslomo.Covariance;
losisi = chol(loomg2);
lobinv = (loBB^(-1));
losi = lobinv*losisi';
loomg = losi*losi'/Dt;
lolol = chol(loomg);

lokpii = lobinv*(eye(3)-loBB)/Dt;
[louk,loek] = eig(lokpii);
louuk = louk^-1;
%locc = lobinv * eslomo.Constant/ Dt;
lovar = diag(loomg);

%AAPL,INTC,QCOM


%following variables naming convention follow previous ones without the 'lo'
mo = varm(3,1);
esmo = estimate(mo,Yem');
BB = esmo.AR{1,1};
omg2 = esmo.Covariance;
sisi = chol(omg2);
binv = (BB^(-1));
si = binv*sisi';
omg = si*si'/Dt;
lol = chol(omg);

kpii = binv*(eye(3)-BB)/Dt;
[uk,ek] = eig(kpii);
uuk = uk^-1;

%mm: estimated mean reversion mean
mm = (eye(3)-BB)^(-1) * esmo.Constant;

var = diag(omg);
end
