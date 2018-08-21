M1 = csvread('1.csv',1,1);
M2 = csvread('2.csv',1,1);
M3 = csvread('3.csv',1,1);
M4 = csvread('4.csv',1,1);

T=1;
Yda = [M1(1:1500,1),M3(1:1500,1)];
Dt = T/500;
[omg,loomg,ek,uk,uuk,mm,sluk,std] = VARMOD(Yda(1:1000,:)',Dt);
[Wtt,bbeta] = EXPATRAD(sluk,Yda(1001:1500,:)',mm,std);
PLPA(sluk,mm,std,Wtt,bbeta,Yda(1001:1500,:)')