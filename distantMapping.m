function y = distantMapping(x, bit)
% 把Vertebi译码时相邻两个码的距离尽可能增大

if bit == 3
    n = [0 1 2 3 7 6 5 4];
    y = n(x + 1);
else
    error('Unregconized bit');
end

if size(x, 2) == 1
    y = reshape(y, [], 1);
elseif size(x, 1) == 1
    y = reshape(y, [], 1);
end