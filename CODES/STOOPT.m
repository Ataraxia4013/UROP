function [pitemp] = STOOPT(Yem,Dt,n,N)
pitemp       = (1/gamma)*(((omega^(-1))*delta)*alphatemp - a*((delta'*omega^(-1)*delta)*((N-n)*Dt*alphatemp - 0.25*trace(A*omega)*((N-n)*Dt)^2)));
end