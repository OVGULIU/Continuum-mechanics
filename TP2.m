%%%%%%%%%%%%%%%%%%%
%%MMC TP2        %%
%%Yann LE GUILLY %%
%%%%%%%%%%%%%%%%%%%
clc; clear;

%%%%%%%%%%%%%%%%%%
%%%parametres:%%%%
%%%%%%%%%%%%%%%%%%

%dimensions
b=0.01; 
dR=0.01; %epaisseur

%materiau
E=210e9; %module d Young acier [pa]
mu=0.27; %coef poisson acier

%pressions
P0=1e5; %pression atmospherique [pa]
dP=5e5; %difference de pression [pa]
P1=P0+dP; %pression a l'interieur du cylindre
Pm=(P1+P0)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lame=(E*mu)/((1+mu)*(1-2*mu)); %coef de lame
a=b+dR;
L=100*b;

R=linspace(b,a,100); %on va se deplacer entre b et a
B=(dP/2*mu)*((a^2*b^2)/(b^2-a^2)); %B est le meme dans tout les cas


%%%%%%%%%%%%%%%%
% ENCASTREMENT %
%%%%%%%%%%%%%%%%

A=(1/(lame+mu))*((b^2+a^2)/(b^2-a^2)*(dP/2)-Pm);
K=0; %encastremement



[sigma1, sigma2, sigma3]=sigma(lame, mu, R, A, B, K);
[Mises, Tresca]=criteres(sigma1, sigma2, sigma3);

figure(1)
plot(R,Mises, 'r')
title('Critere de von Mises en fonction de r');

figure(2)
plot(R, Tresca, 'k')
title('Critere de Tresca en fonction de r');

figure(3) 
hold on;
plot(R, sigma1, 'r')
plot(R, sigma2, 'g')
plot(R, sigma3, 'k')
title('Encastrement: Sigma en fonction r');
legend('Sigma1','Sigma2', 'Sigma3')

ur=((A-K/2)*a)+(B/a^2);
display('Dans le cas de l encastrement, en a, ur vaut:')
display(ur)

uz=K*L;
display('Dans le cas de l encastrement, en a, uz vaut:')
display(uz)

epaisseur=abs(((A-K/2)*a)+(B/a^2)-((A-K/2)*b)+(B/b^2));
display('Dans le cas de l encastrement, l epaisseur finale vaut:')
display(epaisseur)
display('Appuyez sur une touche pour continuer...')
pause();

%%%%%%%%%%%%%%
%Bordslibres %
%%%%%%%%%%%%%%

%B toujours le meme
A=(((b^2+a^2)/(b^2-a^2)-(1/2))*(dP/2)-(3/2)*Pm)*(2/(3*lame+2*mu));
K=(-P0/2-(lame*A/2))/mu;

[sigma1, sigma2, sigma3]=sigma(lame, mu, R, A, B, K);
[Mises, Tresca]=criteres(sigma1, sigma2, sigma3);

figure(1)
plot(R,Mises, 'r')
title('Critere de von Mises en fonction de r');

figure(2)
plot(R, Tresca, 'k')
title('Critere de Tresca en fonction de r');


figure(3) 
hold on;
plot(R, sigma1, 'r')
plot(R, sigma2, 'g')
plot(R, sigma3, 'k')
title('Bords libres: Sigma en fonction r');
legend('Sigma1','Sigma2', 'Sigma3')

figure(4) 
hold on;

for i=1:length(R)
    ENER(i)=lame/2*((((A-K)/2)-(B/R(i)^2))+(((A-K)/2)+(B/R(i).^2))+K)^2+mu*((((A-K)/2)-(B/R(i)^2))+(((A-K)/2)+(B/R(i)^2))+K)^2;
    plot(ENER,'-or');
end
title('Bords libres: energie');


pause();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fermeture par un couvercle %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%B idem
%A=
K=(-b^2*P1+a^2*P0+((2*A*lame*b^2+a^2)/2))*(2/(b^2-a^2));

[sigma1, sigma2, sigma3]=sigma(lame, mu, R, A, B, K);
[Mises, Tresca]=criteres(sigma1, sigma2, sigma3);

figure(1)
plot(R,Mises, 'r')
title('Critere de von Mises en fonction de r');

figure(2)
plot(R, Tresca, 'k')
title('Critere de Tresca en fonction de r');

figure(3) 
hold on;
plot(R, sigma1, 'r')
plot(R, sigma2, 'g')
plot(R, sigma3, 'k')
title('Fermeture par un couvercle: Sigma en fonction r');
legend('Sigma1','Sigma2', 'Sigma3')
