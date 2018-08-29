close all;
TA = 0;
traini = 5000;
testi = 1000;
TD = TA+traini;
T = 1;
N = traini;
Dt = T/N;

[omg,loomg,ek,uk,uuk,mm,sluk,std,louuk,loek,con,BB,omg2] = VARMOD(Yem(:,(TA+1):TD),Dt);
[Wtt,bbeta] = EXPATRAD(sluk,Yem(:,(TA+1):TD+testi),mm,std);
PLPA(sluk,mm,std,Wtt,bbeta,Yem(:,(TA+1):TD+testi))

disp('loomg = ')
disp(loomg)
 