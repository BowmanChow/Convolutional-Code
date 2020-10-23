clear
x  = ones(20,1);
b = 0.5;
rho = 0.95;
NoiseSigma = 0.3;
[y,a,beta,n] = channel(x,b,rho,NoiseSigma);

