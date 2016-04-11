function [Vc]=NumCarSolution(d,method)
global h N
h = 0;
Ec = d*0;
Vc = d*0;
% Ek = fminsearch(@(Ef)(air1(Ef, h, 10)-air1(Ef, h, 2)),28e5);
% for ii = 1:length(d)
%     fprintf('step = %3d.\n',ii);
%     [Ec(1,ii),fval,exitflag,output] = fminsearch(@(Ef)Eq1(d(ii),Ef,method),Ek*1.2);
%     fprintf('Algorithm used                                : %s\nNumber of function evaluation                 : %d\nNumber of iterations taken to find an interval: %d\nNumber of zero-finding iterations             : %d\nExit message                                  : %s\n\n',output.algorithm,output.funcCount,output.intervaliterations,output.iterations,output.message);
%     Vc(1,ii) = d(ii)*Ec(1,ii);
% 
%     [Ec(2,ii),fval,exitflag,output] = fminsearch(@(Ef)Eq2(d(ii),Ef,method),Ek*1.2);
%     fprintf('Algorithm used                                : %s\nNumber of function evaluation                 : %d\nNumber of iterations taken to find an interval: %d\nNumber of zero-finding iterations             : %d\nExit message                                  : %s\n\n',output.algorithm,output.funcCount,output.intervaliterations,output.iterations,output.message);
%     Vc(2,ii) = d(ii)*Ec(2,ii);
% end
options = optimset('Display','final','DiffMinChange',1e-20,'Diagnostics','off');%,'MaxIter',1e4);
for ii = 1:length(d)
    fprintf('step = %3d.\n',ii);
    Ec(1,ii) = fminsearch(@(Ef)Eq1(Ef,d(ii),method),31e5, options);
    Vc(1,ii) = d(ii)*Ec(1,ii);

    Ec(2,ii) = fminsearch(@(Ef)Eq2(Ef,d(ii),method),31e5, options);
    Vc(2,ii) = d(ii)*Ec(2,ii);
end
end

function [answer]=alpha1(E,NN)
global alpha eta EE

ys      = (alpha>=eta).*(alpha-eta)./NN*1e22;
x       = 1./(EE./NN*1e21);
[A B u] = ExpFit(x(401:500),ys(401:500));

A       = A/1e22;
B       = abs(B)/1e21;
% A = 21*kB*T/(1e-2*101325/760);
% B = 421*kB*T/(1e-2*101325/760);

answer = NN*A*exp(-B*NN./E);
end

function [answer]=Integrand1(Ef,method)
global h N
if strcmp(method,'morrowair.m')
    answer = (morrowair(Ef, h, 1) >= morrowair(Ef, h, 2)).*(morrowair(Ef, h, 1)-morrowair(Ef, h, 2))./(morrowair(Ef, h, 4).*Ef);
elseif strcmp(method,'air1.m')
    answer = (air1(Ef, h, 10) >= air1(Ef, h, 2)).*(air1(Ef, h, 10)-air1(Ef, h, 2))./(air1(Ef, h, 11).*Ef);
elseif strcmp(method,'exp')
    answer = alpha1(Ef,N);
else
    error('Wrong method')
end
end 

function [answer] = LHS1(Ef,d,method)
answer = Integrand1(Ef,method)*d;
end

function [answer] = Eq1(Ef,d,method)
global Q
answer = abs(LHS1(Ef,d,method)-log(Q));
end

function [answer]=alpha2(E,NN)
global alpha EE
ys      = alpha./NN*1e22;
x       = 1./(EE./NN*1e21);
[A B u] = ExpFit(x(401:500),ys(401:500));

A       = A/1e22;
B       = abs(B)/1e21;
answer  = NN*A*exp(-B*NN./E);
end

function [answer]=Integrand2(Ef,method)
global h N
if strcmp(method,'morrowair.m')
    answer = morrowair(Ef, h, 1)./(morrowair(Ef, h, 4).*Ef);
elseif strcmp(method,'air1.m')
    answer = air1(Ef, h, 10)./(air1(Ef, h, 11).*Ef);
elseif strcmp(method,'exp')
    answer = alpha2(Ef,N);
else
    error('Wrong method')
end
end

function [answer] = LHS2(Ef,d,method)
answer = Integrand2(Ef,method)*d;
end

function [answer] = Eq2(Ef,d,method)
global Q
answer = abs(LHS2(Ef,d,method)-log(Q));
end
