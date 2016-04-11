clear all
close all
clc
beep off

No              = 2.51e19*1e6;                                             %_m^-3s
N               = No;                                                      %_m^-3s
EE              = (1:500)*No*1e-21;                                        %_V/_m
vi              = air1(EE, 0, 10);                                         %_s^-1
va2             = air1(EE, 0, 2);                                          %_s^-1
yi              = air1(EE, 0, 11);                                         %_m^2/_V/_s

alphaA          = vi ./(yi.*EE);                                           %_m^-1
etaA            = va2./(yi.*EE);                                           %_m^-1

ysA             = (alphaA>=etaA).*(alphaA-etaA)./N*1e22;
xA              = 1./(EE./N*1e21);
[AA BA uA]      = ExpFit(xA(171:300),ysA(171:300));

yiA             = AA*exp(BA*xA);
fprintf('fitting function     : y = %1.2e * exp(%1.2e*x)\nnorm of the residuals: r = %1.2e\n',AA,BA,uA.normr);

% AA             = AA/1e22;
% BA             = abs(BA)/1e21;

vi              = morrowair(EE, 0, 1);                                     %_s^-1
va2             = morrowair(EE, 0, 2);                                     %_s^-1
yi              = morrowair(EE, 0, 4);                                     %_m^2/_V/_s

alphaM          = vi ./(yi.*EE);                                           %_m^-1
etaM            = va2./(yi.*EE);                                           %_m^-1

ysM             = (alphaM>=etaM).*(alphaM-etaM)./N*1e22;
xM              = 1./(EE./N*1e21);
[AM BM uM]      = ExpFit(xM(151:200),ysM(151:200));

yiM             = AM*exp(BM*xM);
fprintf('fitting function     : y = %1.2e * exp(%1.2e*x)\nnorm of the residuals: r = %1.2e\n',AM,BM,uM.normr);

% AM             = AM/1e22;
% BM             = abs(BM)/1e21;

plot(...
    EE/(No*1e-21), (alphaA>=etaA).*(alphaA-etaA)./N*1e22, 'r',...
    EE/(No*1e-21), yiA, 'r--',...
    EE/(No*1e-21), (alphaM>=etaM).*(alphaM-etaM)./N*1e22, 'b',...
    EE/(No*1e-21), yiM, 'b--'...
    )
box on
% set(gca,'Xscale','log','Yscale','log')
set(gca,'YMinorTick','on','XMinorTick','on','TickDir','out')
xlabel('E/N (Td)')
ylabel('\alpha/N  (10^{-18} cm^2)');
legend('air1','air1 fit','Morrowair','Morrowair fit');
legend('location','best')
legend('boxoff')