function y = DeComplexMapping(method, x, bit, HardOrSoft)

if strcmp(HardOrSoft , 'hard')
    func = @harddecode;
elseif strcmp(HardOrSoft,'soft')
    func = @abs;
else
    error('decode method unregconized')
end

if strcmp(method,'circle')
    y = func(x - exp(2i * pi * (0 : 2^bit-1)'/2^bit));
    y = -y;
elseif strcmp(method , 'linear')
    y = func(x - [0:7]');
    y = -y;
else
    error('Method not regconized')
end

function y = harddecode(dis)

y = zeros(size(dis))*Inf;
[~,j] = min(dis);
y(j + size(dis, 1) * (0 : size(dis, 2) - 1) ) = -1;