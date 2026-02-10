%
% Ajuste por polyfit y errores en los coeficientes.
%
% Sintaxis: [coefs,errs] = errormc(x,y,n)
% inputs:
% * x, y  --> datos a ajustar
% * n     --> grado del polinomio de ajuste
% outputs:
% * coefs --> coeficientes del polinomio
% * errs  --> sigmas de cada coeficiente
%

function [coefs,errs] = errormc(x,y,n)

[p,S] = polyfit(x,y,n); % n es el grado del polinomio por el que vamos a ajustar la curva y vs. x
Rinv = inv(S.R);
matriz_de_covarianza = (Rinv*Rinv')*(S.normr).^2./(S.df);
y = zeros(n+1,2);
y(:,1) = p;
for i = 1:n+1
    y(i,2) = sqrt(matriz_de_covarianza(i,i));
end
coefs = y(:,1);
errs = y(:,2);