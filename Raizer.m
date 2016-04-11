function [Vc] = Raizer(d)
global p Q
A     = 11.2509252406;                                                     %_m^-1_Pa^-1
B     = 273.772514187;                                                     %_m^-1_Pa^-1

Vc   = B*p*d./(log(A/log(Q))+log(p*d));                           %_V
% 
% Vc    = B*p*d./(log(A/log(1e4))+log(p*d));                           %_V
% loglog(d,Vc1, 'r', d, Vc2,'b')
end

% function [Vc] = Raizer(d)
% global N Q alpha eta EE
% % A     = 11.2509252406;                                                     %_m^-1_Pa^-1
% % B     = 273.772514187;                                                     %V/_m^-1_Pa^-1
% 
% ys      = (alpha>=eta).*(alpha-eta)./N*1e22;
% x       = 1./(EE./N*1e21);
% [A B u] = ExpFit(x(401:500),ys(401:500));
% 
% yi      = A*exp(B*x);
% fprintf('fitting function     : y = %1.2e * exp(%1.2e*x)\nnorm of the residuals: r = %1.2e\n',A,B,u.normr);
% 
% A             = A/1e22;
% B             = abs(B)/1e21;
% 
% Vc    = B*N*d./(log(A/log(Q))+log(N*d));                           %_V
% 
% plot(1./x,ys,'r',1./x,yi,'b');
% end