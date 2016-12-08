function [Mises, Tresca]=criteres(sigma1, sigma2, sigma3)
%fonction calculant les criteres de von Mises et de Tresca
%les parametres sont les 3 termes diagonaux de la matrice
%des contraintes sigma

Mises=sqrt((1/2)*(sigma1-sigma2).^2+(sigma2+sigma3).^2+(sigma3-sigma1).^2);
Tresca=max([abs(sigma1-sigma2), abs(sigma2-sigma3), abs(sigma3-sigma1)]);