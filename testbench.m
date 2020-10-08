clear all; clc;
%% 发送的数据
rng(1)
info = rand(10000,1)>0.5;
%% 编码参数
nm = [3,4];
Poly = [11,13,15];
n = nm(1);
m = nm(2)-1;
%% 编码
code = ConvEncoder(info,nm,Poly);
% %% 信道 理想
% est = zeros(2^nm(1),length(code)/nm(1));
% for nn = 1:length(code)/nm(1)
%     num = code(nn*nm(1)-nm(1)+1:nn*nm(1));
%     num = bi2de(char(num'+'0'));
%     est(num+1,nn) = 1;
% end

%% 信道 高斯，映射 角度等分
code1 = reshape(code,n,[]);
code1 = code1';
code1 = bi2de(code1);

vol = bin2gray(code1,'psk',2^n);
vol = ComplexMapping("circle", vol, n);

n00 = 0:0.1:3;
for k = 1:length(n00)
n0  = n00(k);

vol_out = channel(vol, 0.1, 0.1, n0 / 2);

vol1 = vol_out.';

est1 = abs(vol1-exp(2i*pi*(0:2^n-1)'/2^n));
est = -est1(bin2gray(0:2^n-1,'psk',2^n)+1,:);



%% 解码
info_out = ConvDecoder(est,nm,Poly);
info_out(1:nm(2)-1) = [];
%% 计算误码率
Error = sum(info_out~=info);
ErrorRate(k) = Error/length(info);
end
plot(n00,ErrorRate)
xlabel('n0')
ylabel('ErrorRate')
title('加性白噪声信道')