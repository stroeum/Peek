function [Vc]=NumSphSolution(R,method)
global h
h = 0;
CC = R*0;
Ec = R*0;
Vc = R*0;

for ii = 1:length(R)
    fprintf('step = %3d.\n',ii);
    [CC(1,ii),fval,exitflag,output] = fzero(@(c)Eq1(R(ii),c,method),.05);
    fprintf('Algorithm used                                : %s\nNumber of function evaluation                 : %d\nNumber of iterations taken to find an interval: %d\nNumber of zero-finding iterations             : %d\nExit message                                  : %s\n\n',output.algorithm,output.funcCount,output.intervaliterations,output.iterations,output.message);
    Ec(1,ii) = E(CC(1,ii),R(ii));
    Vc(1,ii) = R(ii)*Ec(1,ii);

    [CC(2,ii),fval,exitflag,output] = fzero(@(c)Eq2(R(ii),c,method),.05);
    fprintf('Algorithm used                                : %s\nNumber of function evaluation                 : %d\nNumber of iterations taken to find an interval: %d\nNumber of zero-finding iterations             : %d\nExit message                                  : %s\n\n',output.algorithm,output.funcCount,output.intervaliterations,output.iterations,output.message);
    Ec(2,ii) = E(CC(2,ii),R(ii));
    Vc(2,ii) = R(ii)*Ec(2,ii);
end
end

function [answer]=E(c,r)
global Eo delta
answer = Eo*delta*(c./r).^2;
end

function [answer]=alpha1(E,NN)
global alpha eta EE

ys      = (alpha>=eta).*(alpha-eta)./NN*1e22;
x       = 1./(EE./NN*1e21);
[A B u] = ExpFit(x(401:500),ys(401:500));

A             = A/1e22;
B             = abs(B)/1e21;
% A = 21*kB*T/(1e-2*101325/760);
% B = 421*kB*T/(1e-2*101325/760);;

answer = NN*A*exp(-B*NN./E);
end

function [answer]=Integrand1(c,r,method)
global h N
Ef     = E(c,r);
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

function [answer] = LHS1(R,c,method)
answer = quadl(@(r)Integrand1(c,r,method),R,c);%,1e-3);
end

function [answer] = Eq1(R,c,method)
global Q
answer = LHS1(R,c,method)-log(Q);
end

function [answer]=alpha2(E,NN)
global alpha EE
ys      = alpha./NN*1e22;
x       = 1./(EE./NN*1e21);
[A B u] = ExpFit(x(401:500),ys(401:500));

A       = A/1e22;
B       = abs(B)/1e21;
answer = NN*A*exp(-B*NN./E);
end

function [answer]=Integrand2(c,r,method)
global h N
Ef     = E(c,r);
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

function [answer] = LHS2(R,c,method)
answer = quadl(@(r)Integrand2(c,r,method),R,c);%,1e-3);
end

function [answer] = Eq2(R,c,method)
global Q
answer = LHS2(R,c,method)-log(Q);
end
