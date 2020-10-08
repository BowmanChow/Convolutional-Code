%% 三个场景的公共部分：电平映射，发送的数据（运行后面的程序之前必须先运行这一节）
clc;clear;close all;
number_of_channel_use = 1000000;   %信道使用次数
x = zeros(8,1);
for k = 1:8
    x(k) = 100*exp(1i*(k-1)*2*pi/8);  %前面的系数越大，代表信号功率越大（调节这个参数可以在场景一和二获得更好的效果）
end
index = randi([1 8],number_of_channel_use,1);
data = x(index);      %发送的数据，也就是xi


%% 场景一（相当理想的信道）
a1 = 1;
data_received1 = a1 * data + sqrt(0.5)*randn(number_of_channel_use,1) + 1i*sqrt(0.5)*randn(number_of_channel_use,1);

figure(1);
scatter(real(data_received1),imag(data_received1));
figure(2);
histogram(angle(data_received1),1000);

%% 场景二（效果还不错）
a2 = sqrt(0.5)*randn(1,1) + 1i*sqrt(0.5)*randn(1,1);
data_received2 = a2 * data + sqrt(0.5)*randn(number_of_channel_use,1) + 1i*sqrt(0.5)*randn(number_of_channel_use,1);

figure(1);
scatter(real(data_received2),imag(data_received2));
figure(2);
histogram(angle(data_received2),1000);

%% 场景三(加的那个噪声不重要，重要的是前面那个系数ai，它的实部和虚部都是高斯分布，模和幅角也是高斯的)
b = 0.5;rho = 0.95;
beta = zeros(number_of_channel_use,1);
beta(1) = sqrt(0.5)*randn(1,1) + 1i*sqrt(0.5)*randn(1,1);
for i = 2:number_of_channel_use
    beta(i) = rho*beta(i-1)+sqrt(1-rho^2)*(sqrt(0.5)*randn(1,1) + 1i*sqrt(0.5)*randn(1,1));
end
a3 = sqrt(1-b^2) + b*beta;     %将下面一句的注释取消，就能将a3的幅角变为原来的1/5，输出信号就很可观了
%a3 = abs(a3).*exp(1i*1/5*angle(a3));   
data_received3 = a3 .* data + sqrt(0.5)*randn(number_of_channel_use,1) + 1i*sqrt(0.5)*randn(number_of_channel_use,1);

figure(1);
scatter(real(data_received3),imag(data_received3));
figure(2);
subplot(4,1,1);
histogram(angle(a3),1000);
title("a3的幅角");
subplot(4,1,2);
histogram(angle(data),1000);
title("xi的幅角");
subplot(4,1,3);
histogram(angle(a3 .* data),1000);
title("ai*xi的幅角");
subplot(4,1,4);
histogram(angle(data_received3),1000);
title("datareceived3的幅角");