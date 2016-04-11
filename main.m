clear all
close all
clc 
beep  off
warning off

global p h gamma qe Q EE Eo N No alpha eta delta kB T R2

qe              = 1.60218e-19;                                             %_C 
kB              = 1.380658e-23;                                            %_J/_K
T               = 292.386136874;                                           %_K
No              = 2.51e19*1e6;                                             %_m^-3
po              = No*kB*T;                                                 %_Pa
N               = No;                                                      %_m^-3
p               = po;                                                      %_Pa
R2              = 1e0;                                                    %_m

delta           = N/No;                                                    %_
h               = 0;                                                       %_km

EE              = (1:500)*No*1e-21;                                        %_V/_m
% vi              = air1(EE, 0, 10);                                          %_s^-1
% va2             = air1(EE, 0, 2);                                           %_s^-1
% yi              = air1(EE, 0, 11);                                          %_m^2/_V/_s
vi              = morrowair(EE, 0, 1);                                     %_s^-1
va2             = morrowair(EE, 0, 2);                                     %_s^-1
yi              = morrowair(EE, 0, 4);                                     %_m^2/_V/_s


alpha           = vi ./(yi.*EE);                                           %_m^2
eta             = va2./(yi.*EE);                                           %_m^2


Eo              = 31*1e5;                                                  %_V/_m
% Eo              = 24.72*1e5;                                               %_V/m
Ebd             = 260*1e5;                                                 %_V/_m
Q               = 1e4;                                                     %_
gamma           = 1e-2;                                                    %_

d               = logspace(-7,-1, 1001);                                   %_m
R               = logspace(-7,-1, 1001);                                   %_m

Vc_RL                           = LowkeRiousset(R);                        %_V
[VcLowCyl VcLowSph]             = Lowke(R);                                %_V
[VcErf  VcApp A_RP B_RP x ys yi]= RioussetPasko(R);                        %_V
[VcPdLg VcApp A_RP B_RP x ys yi]= GibsonRioussetPasko(R);                  %_V
VcAnaCar                        = Raizer(d);                               %_V
VcNumCar                        = load('VcNumCar.dat');                    %_V
VcAnaCyl                        = VcApp;                                   %_V
VcNumCyl                        = load('VcNumCyl.dat');                    %_V
VcAnaSph                        = VcErf;                                   %_V
VcNumSph                        = load('VcNumSph.dat');                    %_V

% VcNumCar(1:2,:)                 = NumCarSolution(d,'morrowair.m');         %_V
% VcNumCar(3:4,:)                 = NumCarSolution(d,'air1.m');              %_V
% VcNumCar(5:6,:)                 = NumCarSolution(d,'exp');                 %_V
% VcNumCyl(1:2,:)                 = NumCylSolution(R,'morrowair.m');         %_V
% VcNumCyl(3:4,:)                 = NumCylSolution(R,'air1.m');              %_V
% VcNumCyl(5:6,:)                 = NumCylSolution(R,'exp');                 %_V
% VcNumSph(1:2,:)                 = NumSphSolution(R,'morrowair.m');         %_V
% VcNumSph(3:4,:)                 = NumSphSolution(R,'air1.m');              %_V
% VcNumSph(5:6,:)                 = NumSphSolution(R,'exp');                 %_V
% VcAnaSph                        = VcErf;                                   %_V

%%
% save VcNumCar.dat VcNumCar -ascii
% save VcNumCyl.dat VcNumCyl -ascii
% save VcNumSph.dat VcNumSph -ascii

%% Find location of the minimum
iCar(1) = 1;
for ii = 1:length(VcAnaCar)
    if (VcAnaCar(ii) < 0)
        VcAnaCar(ii) = Inf;
    end
end

iCar(1)  = find(VcAnaCar     ==min(VcAnaCar     ));
iCar(2)  = find(VcNumCar(1,:)==min(VcNumCar(1,:)));
iCyl(1)  = find(VcAnaCyl     ==min(VcAnaCyl     ));
iCyl(2)  = find(VcNumCyl(1,:)==min(VcNumCyl(1,:)));
iSph(1)  = find(VcAnaSph     ==min(VcAnaSph     ));
iSph(2)  = find(VcNumSph(1,:)==min(VcNumSph(1,:)));

%% Plot figure
% figure(2);
% loglog(...
%     delta*d,VcNumCar(2,:),'r', delta*d,VcNumCar(1,:),'r--',...
%     delta*d,VcNumCar(4,:),'b', delta*d,VcNumCar(3,:),'b--',...
%     delta*d,VcNumCar(6,:),'g', delta*d,VcNumCar(5,:),'g--',... 
%     delta*d,Raizer(d),'k' ...
% )
% box on
% set(gca,'Xscale','log','Yscale','log','TickDir','out')
% set(gca,'YMinorTick','on','XMinorTick','on')
% xlabel('${\rm\delta}\!R {\rm (m)} ; {\rm\delta}\!d {\rm (m)}$','Interpreter','Latex')
% ylabel('$V_{\rm c} {\rm (V)} ; \varepsilon_{\rm e} {\rm (eV)}$','Interpreter','Latex');
% legend('Num. w/o \nu_{a_2} & morrowair.m','Num. w/ \nu_{a_2} & morrowair.m','Num. w/o \nu_{a_2} & air1.m','Num. w/ \nu_{a_2} & air1.m','Num. w/ exp w/o \nu_{a_2}','Num. w/ exp w/ \nu_{a_2}','Analytical');
% legend('location','best')
% legend('boxoff')

%%
figure(1);
clf(1);
set(gcf,'Units','Normalized','OuterPosition',[0.1 0 .5 1],'Color',[1 1 1]);

subplot(3,1,1)
loglog(...
    delta*d*1e2*760,VcAnaCar,'r-', delta*d*1e2*760,VcNumCar(1,:),'r--',...
    delta*R*1e2*760,VcAnaCyl,'g-', delta*R*1e2*760,VcNumCyl(1,:),'g--',...
    delta*R*1e2*760,VcAnaSph,'b-', delta*R*1e2*760,VcNumSph(1,:),'b--',...
    [delta*d(iCar(1)) delta*d(iCar(2)) delta*R(iCyl(1)) delta*R(iCyl(2)) delta*R(iSph(1)) delta*R(iSph(2))]*1e2*760,[VcAnaCar(  iCar(1)) VcNumCar(1,iCar(2)) VcAnaCyl(  iCyl(1)) VcNumCyl(1,iCyl(2)) VcAnaSph(  iSph(1)) VcNumSph(1,iSph(2))],'kx')
axis([1e-2 1e3 1e2 3e4])
box on
set(gca,'Xscale','log','Yscale','log','TickDir','out')
set(gca,'YMinorTick','on','XMinorTick','on')
xlabel('pR--pd (cm Torr)')
ylabel('V_c (V)');
legend('// Th','// Num','C Th','C Num','O Th','O Num','minima');
legend('location','EastOutside')
legend('boxoff')

%%
subplot(3,1,2)
loglog(...
    delta*d,VcAnaCar     ,'r-',delta*d,VcNumCar(1,:),'r--',...
    delta*R,VcAnaCyl     ,'g-',delta*R,VcNumCyl(1,:),'g--',...
    delta*R,VcAnaSph     ,'b-',delta*R,VcNumSph(1,:),'b--',...
    [delta*d(iCar(1)) delta*d(iCar(2)) delta*R(iCyl(1)) delta*R(iCyl(2)) delta*R(iSph(1)) delta*R(iSph(2))],[VcAnaCar(  iCar(1)) VcNumCar(1,iCar(2)) VcAnaCyl(  iCyl(1)) VcNumCyl(1,iCyl(2)) VcAnaSph(  iSph(1)) VcNumSph(1,iSph(2))],'kx')
axis(    [min(delta*d) max(delta*d) 1e2 1e7])
box on
set(gca,'Xscale','log','Yscale','log','TickDir','out')
set(gca,'YMinorTick','on','XMinorTick','on')
xlabel('\deltaR--\deltad (m)')
ylabel('V_c (V)');
legend('// Th','// Num','C Th','C Num','O Th','O Num','minima');
legend('location','EastOutside')
legend('boxoff')

%%
subplot(3,1,3)
loglog(...
    [min(delta*d) max(delta*d)],[Ebd Ebd],'k:',...
    [min(delta*d) max(delta*d)],[Eo  Eo ],'k--',...
    delta*d,VcAnaCar./d              ,'r-',delta*d,VcNumCar(1,:)./d              ,'r--',...
    delta*R,VcAnaCyl./(R.*log(R2./R)),'g-',delta*d,VcNumCyl(1,:)./(R.*log(R2./R)),'g--',...
    delta*R,VcAnaSph./R              ,'b-',delta*R,VcNumSph(1,:)./R              ,'b--',...
    [delta*d(iCar(1)) delta*d(iCar(2)) delta*R(iCyl(1)) delta*R(iCyl(2)) delta*R(iSph(1)) delta*R(iSph(2))],[VcAnaCar(  iCar(1))./d(iCar(1)) VcNumCar(1,iCar(2))./d(iCar(2)) VcAnaCyl(iCyl(1))./(R(iCyl(1)).*log(R2./R(iCyl(1)))) VcNumCyl(1,(iCyl(2)))./(R(iCyl(2)).*log(R2./R(iCyl(2)))) VcAnaSph(  iSph(1))./R(iSph(1)) VcNumSph(1,iSph(2))./R(iSph(2))],'kx')
axis(    [min(delta*d) max(delta*d) .5*Eo 100*Eo])
box on
set(gca,'Xscale','log','Yscale','log','TickDir','out')
set(gca,'YMinorTick','on','XMinorTick','on')
xlabel('\deltaR--\deltad (m)')
ylabel('E_c (V/m)');
legend('E_{\rm bd}','E_{\rm k}','// Th','// Num','C Th','C Num','O Th','O Num','minima');
legend('location','EastOutside')
legend('boxoff')

%%
figure(2);
clf(2)
set(gcf,'Units','Normalized','OuterPosition',[0.5 0 .5 1],'Color',[1 1 1]);

ys      = (alpha>=eta).*(alpha-eta)./No*1e22;
x       = 1./(EE./No*1e21);
[A B u] = ExpFit(x(401:500),ys(401:500));
A       = A/1e22;
B       = abs(B)/1e21;
alphafit= No*A*exp(-B*No./EE);

B  = 2.08e12*1e4;
alphaLc = No*B*((EE/N).^2-(24.72e5/No).^2);
A    = .012;  
alphaLs = No*A*((EE/N)-(24.72e5/No));

alphaN  = ((EE*1e-5/delta<79.4)*delta.*(0.16053*(EE*1e-5/delta-21.65).^2-2.873) + (EE*1e-5/delta>=79.4)*delta.*(16.7766*EE*1e-5/delta-800.06))*1e2;

plot(EE/No*1e21,(alpha-eta)/No*1e22,'r',EE/No*1e21,alphafit/No*1e22,'g',EE/No*1e21,alphaLc/No*1e22,'b',EE/No*1e21,alphaLs/No*1e22,'b--',EE/No*1e21,alphaN/No*1e22,'k-')
box on
set(gca,'Xscale','linear','Yscale','linear','TickDir','out')
set(gca,'YMinorTick','on','XMinorTick','on')
xlabel('E/N (Td)')
ylabel('\alpha_{eff}/N 10^{-18} (cm^2) ');
legend('Numerical','exp fit','Lowke cylindrical','Lowke spherical','Naidis');
legend('location','best')
legend('boxoff')

%%
figure(3);
clf(3)
set(gcf,'Units','Normalized','OuterPosition',[0.1 0 .5 1],'Color',[1 1 1]);

subplot(2,1,1)
loglog(...
    delta*R,VcAnaCyl,'g-',delta*R,VcNumCyl(1,:),'g--',delta*R,VcLowCyl,'k-',...
    [delta*R(iCyl(1)) delta*R(iCyl(2))],[VcAnaCyl(iCyl(1)) VcNumCyl(1,iCyl(2))],'kx')
axis(    [min(delta*d) max(delta*d) 1e2 1e7])
box on
set(gca,'Xscale','log','Yscale','log','TickDir','out')
set(gca,'YMinorTick','on','XMinorTick','on')
xlabel('\deltaR--\deltad (m)')
ylabel('V_c (V)');
legend('cyl Ana','cyl Num','cyl low','minima');
legend('location','EastOutside')
legend('boxoff')

%%
subplot(2,1,2)
loglog(...
    delta*R,VcAnaSph,'g-',delta*R,VcNumSph(1,:),'g--',delta*R,VcLowSph,'k-',...
    [delta*R(iSph(1)) delta*R(iSph(2))],[VcAnaSph(iSph(1)) VcNumSph(1,iSph(2))],'kx')
axis(    [min(delta*d) max(delta*d) 1e2 1e7])
box on
set(gca,'Xscale','log','Yscale','log','TickDir','out')
set(gca,'YMinorTick','on','XMinorTick','on')
xlabel('\deltaR--\deltad (m)')
ylabel('V_c (V)');
legend('sph Ana','sphl Num','sph low','minima');
legend('location','EastOutside')
legend('boxoff')

%% Export
hgexport(1,'CarCylSph.eps');
hgexport(3,'CylSphLow.eps');

% clf; x = [0 1 1 0 ]; y=[0 0 1 1]; t = patch(x,y,[.5 .5 .5]);
% hatch(t,30,[1 0 0],'--',8,2)     
