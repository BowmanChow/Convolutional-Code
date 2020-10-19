function [y, noise, a_out, beta] = channel(x, b, rho, NoiseSigma, a)

z = random('normal', 0, sqrt(0.5), size(x)) + 1i .* random('normal', 0, sqrt(0.5), size(x));

beta = filter(sqrt(1 - rho.^2), [1 -rho], z, random('normal', 0, sqrt(0.5)) + 1i *  random('normal', 0, sqrt(0.5)));

if isempty(a)
    a_out = sqrt(1 - b.^2) + b .* beta;
else
    a_out = a;
end

noise = random('normal', 0, NoiseSigma, size(x)) + 1i .* random('normal', 0, NoiseSigma, size(x));

y = a_out .* x + noise;