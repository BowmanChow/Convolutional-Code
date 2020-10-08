function y = channel(x, b, rho, NoiseSigma)

z = random('normal', 0, sqrt(0.5), size(x)) + 1i .* random('normal', 0, sqrt(0.5), size(x));

beta = filter(sqrt(1 - rho.^2), [1 -rho], z, random('normal', 0, sqrt(0.5)) + 1i *  random('normal', 0, sqrt(0.5)));

a = sqrt(1 - b.^2) + b .* beta;

n = random('normal', 0, NoiseSigma, size(x)) + 1i .* random('normal', 0, NoiseSigma, size(x));

y = a .* x + n;