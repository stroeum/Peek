function [Vc] = LowkeRiousset(R)
global Eo Q delta No
B    = 2.08e12*1e4;                                                        %_V^-2_m^-2
Vc = Eo*delta*R.*(1+sqrt(log(Q)./(Eo^2/No*delta*R*B))).^2;                 %_V
end