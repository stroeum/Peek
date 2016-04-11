function [VcPdLog VcApp A B x ys yi] = GibsonRioussetPasko(R)
global Eo Q delta N No alpha eta EE kB T R2

ys      = (alpha>=eta).*(alpha-eta)./N*1e22;
x       = 1./(EE./N*1e21);
[A B u] = ExpFit(x(401:500),ys(401:500));

yi      = A*exp(B*x);
fprintf('fitting function     : y = %1.2e * exp(%1.2e*x)\nnorm of the residuals: r = %1.2e\n',A,B,u.normr);

A             = A/1e22;
B             = abs(B)/1e21;

lambda = exp(-B*No/Eo);
mu     =  A*N/log(Q);
c      = R*log(lambda)./(mu*lambda*R-lambertw(mu*R.*exp(mu*lambda*R)));

VcPdLog= Eo*delta*c.*log(R2./R);

VcApp  = B*(A*N*R+log(Q)).*log(R2./R)/(A*(1-exp(-B*No/Eo)));

% Ec = B*(A*N*R+log(Q))./(R*A*(1-exp(-B*No/Eo)));
end