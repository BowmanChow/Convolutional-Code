function y = DeComplexMapping(method, x, bit)

if method == 'circle'
    y = abs(x - exp(2i * pi * (0 : 2^bit-1)'/2^bit));
    y = -y;
elseif method == 'linear'
    y = abs(x - [0:7]');
    y = -y;
else
    error('Method not regconized')
end