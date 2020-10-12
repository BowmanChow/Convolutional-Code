clear all; clc;
%% 发�?�的数据
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

%% 信道 高斯，映�? 角度等分
code1 = reshape(code,n,[]);
code1 = code1';
code1 = bi2de(code1);

gray = bin2gray(code1,'psk',2^n);
% gray = distantMapping(code1, n);
vol = ComplexMapping('circle', gray, n);

n00 = 1:0.1:5;
SNR = [];
for k = 1:length(n00)
n0  = n00(k);

[vol_out, noise] = channel(vol, 0.1, 0.1, n0 / 2);

vol1 = vol_out.';

est1 = DeComplexMapping('circle', vol1, n);
est = est1(bin2gray(0:2^n-1,'psk',2^n)+1,:);
%est = est1(distantMapping(0:2^n-1, n)+1,:);



%% 解码
info_out = ConvDecoder(est,nm,Poly);
info_out(1:nm(2)-1) = [];
%% 计算误码�?
Error = sum(info_out~=info);
ErrorRate(k) = Error/length(info);
%% SNR
SNR = [SNR; mean(abs(vol).^2) / mean(abs(noise).^2)];
end
plot(SNR,ErrorRate)
xlabel('snr')
ylabel('ErrorRate')
title('加�?�白噪声信道')