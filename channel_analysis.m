%% ���������Ĺ������֣���ƽӳ�䣬���͵����ݣ����к���ĳ���֮ǰ������������һ�ڣ�
clc;clear;close all;
number_of_channel_use = 1000000;   %�ŵ�ʹ�ô���
x = zeros(8,1);
for k = 1:8
    x(k) = 100*exp(1i*(k-1)*2*pi/8);  %ǰ���ϵ��Խ�󣬴����źŹ���Խ�󣨵���������������ڳ���һ�Ͷ���ø��õ�Ч����
end
index = randi([1 8],number_of_channel_use,1);
data = x(index);      %���͵����ݣ�Ҳ����xi


%% ����һ���൱������ŵ���
a1 = 1;
data_received1 = a1 * data + sqrt(0.5)*randn(number_of_channel_use,1) + 1i*sqrt(0.5)*randn(number_of_channel_use,1);

figure(1);
scatter(real(data_received1),imag(data_received1));
figure(2);
histogram(angle(data_received1),1000);

%% ��������Ч��������
a2 = sqrt(0.5)*randn(1,1) + 1i*sqrt(0.5)*randn(1,1);
data_received2 = a2 * data + sqrt(0.5)*randn(number_of_channel_use,1) + 1i*sqrt(0.5)*randn(number_of_channel_use,1);

figure(1);
scatter(real(data_received2),imag(data_received2));
figure(2);
histogram(angle(data_received2),1000);

%% ������(�ӵ��Ǹ���������Ҫ����Ҫ����ǰ���Ǹ�ϵ��ai������ʵ�����鲿���Ǹ�˹�ֲ���ģ�ͷ���Ҳ�Ǹ�˹��)
b = 0.5;rho = 0.95;
beta = zeros(number_of_channel_use,1);
beta(1) = sqrt(0.5)*randn(1,1) + 1i*sqrt(0.5)*randn(1,1);
for i = 2:number_of_channel_use
    beta(i) = rho*beta(i-1)+sqrt(1-rho^2)*(sqrt(0.5)*randn(1,1) + 1i*sqrt(0.5)*randn(1,1));
end
a3 = sqrt(1-b^2) + b*beta;     %������һ���ע��ȡ�������ܽ�a3�ķ��Ǳ�Ϊԭ����1/5������źžͺܿɹ���
%a3 = abs(a3).*exp(1i*1/5*angle(a3));   
data_received3 = a3 .* data + sqrt(0.5)*randn(number_of_channel_use,1) + 1i*sqrt(0.5)*randn(number_of_channel_use,1);

figure(1);
scatter(real(data_received3),imag(data_received3));
figure(2);
subplot(4,1,1);
histogram(angle(a3),1000);
title("a3�ķ���");
subplot(4,1,2);
histogram(angle(data),1000);
title("xi�ķ���");
subplot(4,1,3);
histogram(angle(a3 .* data),1000);
title("ai*xi�ķ���");
subplot(4,1,4);
histogram(angle(data_received3),1000);
title("datareceived3�ķ���");