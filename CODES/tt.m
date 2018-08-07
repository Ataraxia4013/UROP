
mo = varm(3,1);
esmo = estimate(mo,log(Yem'));
BB = esmo.AR{1,1};
omg2 = esmo.Covariance;
sisi = chol(omg2);
binv = (BB^(-1));
si = binv*sisi'/Dt;
omg = si*si';
lol = chol(omg);
cc = binv * esmo.Constant/ Dt;
BOI = binv*(eye(3)-BB)/Dt;
var = diag(omg);
deeet = -(0.5 * [0.04 0.0225 0.0015])';
