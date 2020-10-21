clear all; clc;
%% 发送的数据
% rng(1)
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

gray = bin2gray(code1,'psk',2^n);
% gray = distantMapping(code1, n);
vol = ComplexMapping('circle', gray, n);

%% scenes
scene = 3;

%% generate a
if scene == 4
    b = 1; rho = 1;
elseif scene == 3
    b = 0.5; rho = 0.95;
elseif scene == 2
    b = 0.1; rho = 0.1;
elseif scene == 1
    b = 0; rho = 0;
end
[~,~,a, beta] = channel(vol, b, rho, 1, []); 

%% ifKnowA == 0 means both dont know a, ifKnowA == 1 means reciever know a, ifKnowA == 2 means both know a
ifKnowA = 0;
if ifKnowA == 2
    vol = vol ./ (a ./ abs(a));
elseif ifKnowA == 0
    if scene == 4               %  b = 1, rho = 1 , need to transfer 1 before sending
        add_ones_num = 100;
        vol = [ones(add_ones_num,1);vol];
    elseif scene == 3           % b = 0.5, rho = 0.95, insert 1 cross
        tmp = zeros(2 * length(vol),1);
        tmp(1:2:end) = 1;
        tmp(2:2:end) = vol;
        vol = tmp;clear tmp;
    end
    [~,~,a, beta] = channel(vol, b, rho, 1, []); 
end

%%

n00 = 1 ./ sqrt([1.3:0.01:2.7]);
SNR = zeros(1,length(n00));ErrorRate = zeros(1,length(n00));
for k = 1:length(n00)
n0  = n00(k);

[vol_out, noise] = channel(vol, b, rho, n0 / 2, a);

%% ifKnowA == 0 means both dont know a, ifKnowA == 1 means reciever know a, ifKnowA == 2 means both know a
if ifKnowA == 2
    vol_out = vol_out ./ abs(a);
elseif ifKnowA == 1
    vol_out = vol_out ./ a;
else
    if scene == 3               % b = 0.5, rho = 0.95, insert 1 cross
        tmp = reshape(vol_out', 2, []);
        kalman_beta = KalmanFilter(tmp(1,:) - sqrt(1-b^2), rho^2, rho^2*(1-rho^2), b, 2*(n0 / 2)^2, beta(1));
        vol_out = tmp(2,:)' ./ (sqrt(1-b^2) + b * kalman_beta);
        clear tmp;
    elseif scene == 4           % b = 1, rho = 1 , need to decode 1 before recieving
        vol_out_mean = mean(vol_out(1:add_ones_num));
        vol_out_store = vol_out;
        vol_out = vol_out ./ vol_out_mean;
        vol_out(1:add_ones_num) = [];
    end
end

%%
vol_out = vol_out.';

est1 = DeComplexMapping('circle', vol_out, n, 'soft');
est = est1(bin2gray(0:2^n-1,'psk',2^n)+1,:);
%est = est1(distantMapping(0:2^n-1, n)+1,:);


%% 解码
info_out = ConvDecoder(est,nm,Poly);
info_out(1:nm(2)-1) = [];
%% 计算误码率
Error = sum(info_out~=info);
ErrorRate(k) = Error/length(info);
%% SNR
SNR(k) = mean(abs(vol).^2) / mean(abs(noise).^2);
end

semilogy(SNR,ErrorRate,'.');
hold on;
xq = SNR(1):0.01:SNR(end);
ModelFunc = @(p, x) exp(-p(1) .* x + p(2));
fitline = fitnlm(SNR,ErrorRate,ModelFunc,[1 1]);
semilogy(xq, predict(fitline, xq'));
xlabel('snr')
ylabel('ErrorRate')
title('加性白噪声信道')