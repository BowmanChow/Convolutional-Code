clear all; clc;
%% 发�?�的数据
rng(1)
info = rand(10000,1)>0.5;
%% 编码参数
% nm = [3,9];
% Poly = [367   435   457];
nm = [3,4];
Poly = [11 13 15];
n = nm(1);
m = nm(2)-1;

%% CRC 
CRC_Unit = 200;
% Poly1 = [1 zeros(1,CRC_Unit-1) 1];
% nn = 1;
% while length(Poly1)>1
%     Res1 = [1 1];
%     mm = 1;
%     Poly2 = Poly1;
%     while sum(Res1)>0 && length(Poly2)>length(Poly1)/2-2
%         mm=mm+1;
%         [Poly2,Res1] = deconvm2(Poly1,flip(de2bi(mm)));
%         Poly2 = mod(Poly2,2);
%         Res1 = mod(Res1,2);
% %         Res1 = bi2de(flip(Res1));
%     end   
% %    Poly2 = bi2de(flip(Poly2));
%     if length(Poly2)>length(Poly1)/2-2
%     Divisor(nn) = mm;
%     Poly1 = Poly2;
%     else 
%         Divisor(nn) = bi2de(flip(Poly1));
%         Poly1 = 0;
%     end
%     nn=nn+1;
%     
% end

CRC_RLength = 8;
CRC_Checks = [1 1 0 1 0 0 1 1 1];
infoCRC = zeros(ceil(length(info)/CRC_Unit)*CRC_RLength+length(info),1);

for mm = 0:(length(info)-1)/CRC_Unit
    temp = info(mm*CRC_Unit+1:min((mm+1)*(CRC_Unit),length(info)));
    [~,r1] = deconvm2([temp;zeros(CRC_RLength,1)],CRC_Checks);
    infoCRC(mm*(CRC_Unit+CRC_RLength)+1:min((mm+1)*(CRC_Unit+CRC_RLength),length(infoCRC))) = ... 
        [temp;r1(max(end-CRC_RLength+1,0):end)];
end

% info = infoCRC;
%% 编码
code = ConvEncoder(infoCRC,nm,Poly);


%% 信道 高斯，映�? 角度等分
code1 = reshape(code,n,[]);
code1 = code1';
code1 = bi2de(code1);

gray = bin2gray(code1,'psk',2^n);
% gray = distantMapping(code1, n);
vol = ComplexMapping('circle', gray, n);

n00 = 0.5:0.1:2;
SNR = [];
for k = 1:length(n00)
n0  = n00(k);

[vol_out, noise] = channel(vol, 0.1, 0.1, n0 / 2);

vol1 = vol_out.';

est1 = DeComplexMapping('circle', vol1, n, 'soft');
est = est1(bin2gray(0:2^n-1,'psk',2^n)+1,:);
%est = est1(distantMapping(0:2^n-1, n)+1,:);



%% 解码
info_out = ConvDecoder(est,nm,Poly);
info_out(1:nm(2)-1) = [];

%% 解CRC
BlockError = 0;
info_out1 = zeros(size(info));
for mm = 0:(length(info_out)-1)/(CRC_Unit+CRC_RLength)
    temp = info_out(mm*(CRC_Unit+CRC_RLength)+1:min(length(info_out),(mm+1)*(CRC_Unit+CRC_RLength)));
    [~,r1] = deconvm2(temp,CRC_Checks);
    info_out1(mm*CRC_Unit+(1:(length(temp)-CRC_RLength)))=temp(1:(length(temp)-CRC_RLength));
    BlockError = BlockError+(sum(r1)>0);
end
BlockErrorRate(k) = BlockError/(mm+1);
%% 计算误码�?
Error = sum(info_out1~=info);
ErrorRate(k) = Error/length(info);
%% SNR
SNR = [SNR; mean(abs(vol).^2) / mean(abs(noise).^2)];


end
plot(SNR,ErrorRate)
xlabel('snr')
ylabel('ErrorRate')
hold on
plot(SNR,BlockErrorRate);
legend('Bit Error Rate','Block Error Rate')
hold off

title('加性白噪声信道')