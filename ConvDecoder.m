function info_out = ConvDecoder(est,nm,Poly,WindowSize)
% est 2^n by t, nm [n,m+1], Poly 10进制, WindowSize 滑动窗的尺寸
n = nm(1);
mp = nm(2);
m = mp-1;

if size(est,2)==2^n
    est = est';
end

if nargin == 3
    WindowSize = size(est,2);
end
    
Poly = reshape(Poly,[],1);
Poly = flip(de2bi(Poly,mp),2);
%Poly = logical(Poly>'0');

%低位均表示时刻靠前的，t=0, 一定在code(1)， 当前时刻在code(last)
% 十进制转二进制后，最高位对应a(1)，是时间最靠前的位。
StatePrev = [0:2^m-1;0:2^m-1]'; % 当前状态来自之前的哪个状态 行号表示状态序号+1，列表示m时刻之前是0还是1
CodeCurrent = StatePrev*2+[0,1]; % 该状态在当前(偏后)时刻看到的码 m+1位
StatePrev = mod(CodeCurrent,2^m);% 之前的状态为[0或1 State的高m-1位]

CodeCurrent_1 = de2bi(reshape(CodeCurrent,[],1),mp);
CodeIn = CodeCurrent_1*(Poly)';
CodeIn = mod(CodeIn,2);
CodeIn = bi2de(CodeIn);
CodeIn = reshape(CodeIn,[],2);

BestCode = zeros(2^m,size(est,2)+1); %当前状态认为的最优上一状态，0或1表示m时刻之前的输入为0或1
EstSum = zeros(2^m,size(est,2)+1);
EstSum(2:end,1) = -Inf;

for nn = 1:size(est,2)
    EstPrev = EstSum(:,nn);
    EstCurr = est(:,nn);
    [EstSum(:,nn+1),BestCode(:,nn+1)] = max(EstPrev(StatePrev+1)+EstCurr(CodeIn+1),[],2);
    
end

info_out = zeros(size(est,2),1);
State = 0;
for nn = size(est,2):-1:1
    info_out(nn) = BestCode(State+1,nn+1)-1;
    State = StatePrev(State+1,BestCode(State+1,nn+1));
    
end




end