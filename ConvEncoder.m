function out = ConvEncoder(info_in, nm, Poly)
% info_in 为输入的信息，double型列向量
% nm = [n,m+1]，n为一周期码长，m为寄存器数目,仅对k=1情形有效
% Poly 多项式形式，十进制
n = nm(1);
mp = nm(2);%m+1
if length(Poly)~=n
    msgID = 'myComponent:inputError';
    msgtext = 'Wrong size of Ploynomial.';
    ME = MException(msgID,msgtext);
    throw(ME)
elseif ~isempty(find(Poly>2^mp,1))
    msgID = 'myComponent:inputError';
    msgtext = 'Some element in Poly is larger than 2^(m+1).';
    ME = MException(msgID,msgtext);
    throw(ME)
end

Poly = reshape(Poly,[],1);
Poly = de2bi(Poly,mp);
% Poly = logical(Poly>'0');
%Poly = reshape(Poly',1,mp,n);
Poly = Poly';

cellsize = 1000;%1000位截断一次
info_in = [zeros(mp-1,1);info_in;zeros(mp-1,1)];
out = [];
info_in1 = info_in;
    for nn = 1:n
        tmp(:,nn) = conv(info_in1,(Poly(:,nn)),'Valid');
    end
    tmp = reshape(tmp.',[],1);
    out = mod(tmp,2);

% while length(info_in)>=mp
% idx = zeros(min(cellsize,size(info_in,1)+1-mp),mp);
% cellsize_1 = size(idx,1);
% for kk = 1:mp
%     idx(:,kk) = kk:kk+cellsize_1-1;
% end
%     cell = info_in(idx);
%     code = cell*Poly;
%     code = reshape(code',[],1);
%     code = mod(code,2);
%     out = [out;code];
%     info_in(1:cellsize_1) = [];
% %     tmp = zeros(min(cellsize+mp-1,length(info_in)),n);
% %     for nn = 1:n
% %         tmp(:,n) = conv(info_in(1:min(cellsize+mp-1,length(info_in))),flip(Poly(:,1)),'Valid');
% %     end
% %     tmp = reshape(tmp.',[],1);
% %     tmp = mod(tmp,2);
% %     out = [out;code];
% %     info_in(1:cellsize_1) = [];
%     
% end


end