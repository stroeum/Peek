function [VcErf VcApp A B x ys yi] = RioussetPasko(R)
global Eo Q delta N No alpha eta EE kB T

ys      = (alpha>=eta).*(alpha-eta)./N*1e22;
x       = 1./(EE./N*1e21);
[A B u] = ExpFit(x(401:500),ys(401:500));

yi      = A*exp(B*x);
fprintf('fitting function     : y = %1.2e * exp(%1.2e*x)\nnorm of the residuals: r = %1.2e\n',A,B,u.normr);

A             = A/1e22;
B             = abs(B)/1e21;

% A = 21*kB*T/(1e-2*101325/760);
% B = 421*kB*T/(1e-2*101325/760);

VcApp(1,:) = Eo * delta .* R ./ R.^2 .* ...
    (...
        (log(Q)/(A*sqrt(Eo*No/B)*sqrt(pi/2)) + 2/sqrt(pi)*delta*R*sqrt(B/(Eo/No))) / ...
        (2/sqrt(pi)*(sqrt(B/(Eo/No))*delta)) ...
    ).^2;
VcApp(2,:) = Eo * delta * R .* (1 + log(Q)./(A*No*delta*R)).^2;

VcErf = Eo * delta .* R ./ R.^2 .* ... 
    (...
        2/sqrt(pi)*sqrt(B/(Eo/No))*(log(Q)./(A*N) + R) ...
    ).^2;
end