function [sigma1, sigma2, sigma3]=sigma(lame, mu, R, A, B, K)
%fonction calculant les termes diagonaux de la matrice sigma
%Les paramètres a renseigner sont:
%lame: le coef de lame, mu: le coef de poisson, R: la position,
%A, B, K: voire TP2 MMC

sigma1=lame*A+2*mu*((A-K)/2)-(B./R.^2);
sigma2=lame*A+2*((A-K)/2+(B./R.^2));
sigma3=lame*A+2*mu*K;