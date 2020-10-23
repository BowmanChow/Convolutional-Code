function [q,r]=deconvm2(b,a)
la = length(a);
lb = length(b);
if size(b,1)~=1
    b0 = b;
    b = b';
end
% q = zeros(1,lb-la+1);
r1 = [];
q = [];
for curr = 1:50:lb
    b1 = [r1(max(end-la+2,1):end) b(curr:min(curr+49,lb))];
    [q1,r1] = deconv(b1,a);
    q1 = mod(q1,2);
    r1 = mod(r1,2);
    q = [q q1];
end

r = r1;

if size(b0,1)~=1
    r = r';
    q = q';
end

end