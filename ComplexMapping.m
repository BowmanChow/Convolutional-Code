function y = ComplexMapping(method , x, bit)

if method == 'circle'
    y = exp(2i * pi * x / 2^bit);
elseif method == 'linear'
    y = x;
else
    error('Method not regconized')
end
