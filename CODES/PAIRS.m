close all;
M1 = csvread('s1.csv',1,1);
M2 = csvread('s2.csv',1,1);
M3 = csvread('s3.csv',1,1);

Yda = [M1(1:65000,1),M2(1:65000,1)];

TA = 52;
traini = 20;
testi = 1;
TD = TA+traini;
T = 1;
N = 437*traini;
Dt = T/N;

[omg,loomg,ek,uk,uuk,mm,sluk,std,louuk,loek,con,BB,omg2] = VARMOD(Yda(437*TA+1:437*TD,:)',Dt);

DDD = Des(Yda(437*TA+1:437*TD,:)');
[mu,esvarm] = karl(con,BB,omg2,DDD,Yda(437*TA+1:437*TD,:)');
[Wtt,bbeta] = EXPATRAD(sluk,Yda(437*TD+1:437*(TD+testi),:)',mu,std,200);
PLPA(sluk,mu,std,Wtt,bbeta,Yda(437*TD+1:437*(TD+testi),:)')

disp('ek = ')
disp(ek)
disp('sluk = ')
disp(sluk)
