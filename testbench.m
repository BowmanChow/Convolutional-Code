clear all; clc;close all
%% 参数
nm = [3,4];                 % n，m
Poly = [11,13,15];          % Polynomial
scene = 1;                  % scene; 1-4
ifKnowA = 2;                % ifKnowA == 0 means both dont know a, ifKnowA == 1 means reciever know a, ifKnowA == 2 means both know a
n00 = 1 ./ sqrt([1.3:0.01:2.7]);
n00 =  sqrt(0.2:0.01:1.3);    % noise
L0 = 10000;                 % length of infomation
rep = 4;                    % repeat times
CRC = 1;
CRC_Unit = 200;
CRC_RLength = 8;
CRC_Checks = [1 1 0 1 0 0 1 1 1];


SNR = zeros(rep,length(n00));
ErrorRate = zeros(rep,length(n00));
BlockErrorRate = zeros(rep,length(n00));
for r = 1:rep

%% 发送的数据
% rng(1)
info = rand(L0,1)>0.5;
%% 编码参数
n = nm(1);
m = nm(2)-1;
%% 编码
if CRC
infoCRC = zeros(ceil(length(info)/CRC_Unit)*CRC_RLength+length(info),1);

    for mm = 0:(length(info)-1)/CRC_Unit
        temp = info(mm*CRC_Unit+1:min((mm+1)*(CRC_Unit),length(info)));
        [~,r1] = deconvm2([temp;zeros(CRC_RLength,1)],CRC_Checks);
        infoCRC(mm*(CRC_Unit+CRC_RLength)+1:min((mm+1)*(CRC_Unit+CRC_RLength),length(infoCRC))) = ... 
            [temp;r1(max(end-CRC_RLength+1,0):end)];
    end
else
    infoCRC = info;
end
code = ConvEncoder(infoCRC,nm,Poly);

%% 信道 高斯，映射 角度等分
code1 = reshape(code,n,[]);
code1 = code1.';
code1 = bi2de(code1);

gray = bin2gray(code1,'psk',2^n);
% gray = distantMapping(code1, n);
vol = ComplexMapping('circle', gray, n);

%% scenes
if scene == 4 %  b = 1, rho = 1 , need to transfer 1 before sending
    add_ones_num = 100;
    vol = [ones(add_ones_num,1);vol];
elseif scene == 3 % b = 0.5, rho = 0.95, insert 1 cross
    tmp = zeros(2 * length(vol),1);
    tmp(1:2:end) = 1;
    tmp(2:2:end) = vol;
    vol = tmp;clear tmp;
end

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
if ifKnowA == 2
    vol = vol ./ a;
    vol = vol ./ abs(vol);
end   

%%
for k = 1:length(n00)
n0  = n00(k);

[vol_out, noise] = channel(vol, b, rho, n0 / 2, a);
%% ifKnowA == 0 means both dont know a, ifKnowA == 1 means reciever know a, ifKnowA == 2 means both know a
if ifKnowA == 0
    vol_out = vol_out ./ abs(a);
elseif ifKnowA == 1
    vol_out = vol_out ./ a;
end

%%
if scene == 3 % b = 0.5, rho = 0.95, insert 1 cross
    tmp = reshape(vol_out.', 2, []);
    kalman_beta = KalmanFilter(tmp(1,:) - sqrt(1-b^2), rho^2, rho^2*(1-rho^2), b, 2*(n0 / 2)^2, beta(1));
    vol_out = tmp(2,:).' ./ (sqrt(1-b^2) + b * kalman_beta);
    clear tmp;
elseif scene == 4 % b = 1, rho = 1 , need to decode 1 before recieving
    vol_out = vol_out ./ mean(vol_out(1:add_ones_num));
    vol_out(1:add_ones_num) = [];
end

%%
vol_out = vol_out.';

est1 = DeComplexMapping('circle', vol_out, n, 'soft');
est = est1(bin2gray(0:2^n-1,'psk',2^n)+1,:);
%est = est1(distantMapping(0:2^n-1, n)+1,:);



%% 解码
info_out = ConvDecoder(est,nm,Poly);
info_out(1:nm(2)-1) = [];
%% 解CRC
BlockError = 0;
if CRC
    info_out1 = zeros(size(info));
for mm = 0:(length(info_out)-1)/(CRC_Unit+CRC_RLength)
    temp = info_out(mm*(CRC_Unit+CRC_RLength)+1:min(length(info_out),(mm+1)*(CRC_Unit+CRC_RLength)));
    [~,r1] = deconvm2(temp,CRC_Checks);
    info_out1(mm*CRC_Unit+(1:(length(temp)-CRC_RLength)))=temp(1:(length(temp)-CRC_RLength));
    BlockError = BlockError+(sum(r1)>0);
end
else
    info_out1 = info_out;
    mm = 0;
end
BlockErrorRate(r,k) = BlockError/(mm+1);
%% 计算误码率
Error = sum(info_out1~=info);
ErrorRate(r,k) = Error/length(info);
%% SNR
SNR(r,k) = mean(abs(vol).^2) / mean(abs(noise).^2);
end
end
SNR = 10*log10(mean(SNR));
ErrorRate = mean(ErrorRate)+1e-16;
semilogy(SNR,ErrorRate)
hold on;
xq = SNR(1):0.01:SNR(end);
ModelFunc = @(p, x) exp(-p(1) .* x + p(2));
fitline = fitnlm(SNR,ErrorRate,ModelFunc,[1 1]);
semilogy(xq, predict(fitline, xq'));
xlabel('SNR/dB')
ylabel('ErrorRate')
title(sprintf('场景%d,条件%d',scene,ifKnowA+1))

figure;
BlockErrorRate = mean(BlockErrorRate)+1e-16;
semilogy(SNR,BlockErrorRate)
hold on;
xq = SNR(1):0.01:SNR(end);
ModelFunc = @(p, x) exp(-p(1) .* x + p(2));
fitline = fitnlm(SNR,BlockErrorRate,ModelFunc,[1 1]);
semilogy(xq, predict(fitline, xq'));
xlabel('SNR/dB')
ylabel('BlockErrorRate')
title(sprintf('场景%d,条件%d',scene,ifKnowA+1))
