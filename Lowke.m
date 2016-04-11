function [Vcyl Vsph] = Lowke(R)
% global Eo Q delta No R2
global Q R2 delta N No Eo

% Eo = 24.72*1e5;                                                            %_V/m
% No = 2.51e19*1e6;                                                          %_m^-3
% N  = No;                                                                   %_m^-3
B  = 2.08e12*1e4;                                                          %_V^2_m^-2
A    = .012;                                                               %_V^-1

Vcyl = Eo*delta*R.*log(R2./R).*(1+sqrt(log(Q)./(B*N*R*(Eo/No)^2)));        %_V
Vsph = Eo*delta*R.*(1+sqrt(log(Q)./(A*N*R*(Eo/No)))).^2;                   %_V

end