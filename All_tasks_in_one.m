clc;
close all;
clear;

%%%Task2%%%
%%%%%%%%%%%

fs = 48000;
fc = 3000;
% N = 55;
N = 47;
m = (N-1)/2; % 延迟 m 使得 h_highpass[n] 是类型 I 线性相位滤波器
F = 12;

% 计算理想滤波器的冲击响应
n = -(N-1)/2 : (N-1)/2;
hd = (2 * fc / fs) * sinc(2 * fc / fs * n);

% 窗函数，chebyshev>hamming>blackman>hann
% w = blackman(N)';
% w = hamming(N)';
w = chebwin(N, 45)';
% 45 is recommended
h = hd .* w;

% 确保DC value为1
h = h / sum(h);

% 频率响应分析
[H_quant, f] = freqz(h, 1, 1024, fs);
normalized_f = f/fs;
H_db = 20*log10(abs(H_quant));

% 频率响应图
figure;
plot(normalized_f, H_db)
grid on;
title('Frequency Response of the Low-Pass Filter');
xlabel('Normalized Frequency f/f_s');
ylabel('Magnitude (dB)');
hold on;
index_1_8 = find(normalized_f == 1/8, 1);
index_1_16 = find(normalized_f == 1/16, 1);
plot(normalized_f(index_1_8), H_db(index_1_8), 'r*'); 
plot(normalized_f(index_1_16), H_db(index_1_16), 'g*'); 
text(normalized_f(index_1_8), H_db(index_1_8), [' \leftarrow (1/8, ' num2str(H_db(index_1_8)) ' dB)'], 'HorizontalAlignment', 'left');
text(normalized_f(index_1_16), H_db(index_1_16), [' \leftarrow (1/16, ' num2str(H_db(index_1_16)) ' dB)'], 'HorizontalAlignment', 'left');
hold off;

figure;
plot(normalized_f, abs(H_quant))
title('Frequency Response of the Low-Pass Filter');
xlabel('Normalized Frequency f/f_s');
ylabel('Magnitude'); 
grid on;

figure;
stem(0:N-1, h, "filled")
title('Impulse Response of the Low-Pass Filter');
xlabel('n');
ylabel('Magnitude');
grid on;

%%%Task3%%%
%%%%%%%%%%%
% 创建单位冲击函数
delta = [zeros(1, m) 1 zeros(1, m)]; 

h_highpass = delta - h;

[H_high, f_high] = freqz(h_highpass, 1, 1024, fs);
normalized_f_high = f_high/fs;
H_db_high = 20*log10(abs(H_high));

% 频率响应图
figure;
plot(normalized_f_high, H_db_high)
grid on;
title('Frequency Response of the High-Pass Filter');
ylim([-100,50])
line([1/16, 1/16], ylim, 'Color', 'cyan', 'LineStyle', '-.', 'DisplayName', 'v = 1/16');
line([1/8, 1/8], ylim, 'Color', 'cyan', 'LineStyle', '-', 'DisplayName', 'v = 1/8');
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');


% figure;
% plot(normalized_f_high, abs(H_high))
% grid on;
% title('Frequency Response of the High-Pass Filter');
% xlabel('Normalized Frequency');
% ylabel('Magnitude');

figure;
stem(0:N-1, h_highpass, "filled")
title('Impulse Response of the High-Pass Filter');
xlabel('n');
ylabel('Magnitude');
grid on;

%%%Task4%%%
%%%%%%%%%%%

h_quant = 2^(-F) * round(h * 2^F);
[H_quant, f] = freqz(h_quant, 1, 1024, fs);
normalized_f_quant = f/fs;
H_quant_db = 20*log10(abs(H_quant));

% figure;
% plot(normalized_f_quant, H_quant_db)
% grid on;
% title('Frequency Response of the Low-Pass Filter with Quantization');
% xlabel('Normalized Frequency f/f_s');
% ylabel('Magnitude (dB)');
% hold on;
% [~, index_1_8] = min(abs(normalized_f_quant - 1/8));
% [~, index_1_16] = min(abs(normalized_f_quant - 1/16));
% plot(normalized_f_quant(index_1_8), H_quant_db(index_1_8), 'r*'); % 红色圆圈
% plot(normalized_f_quant(index_1_16), H_quant_db(index_1_16), 'g*'); % 绿色星号
% text(normalized_f_quant(index_1_8), H_quant_db(index_1_8), [' \leftarrow (1/8, ' num2str(H_quant_db(index_1_8)) ' dB)'], 'HorizontalAlignment', 'left');
% text(normalized_f_quant(index_1_16), H_quant_db(index_1_16), [' \leftarrow (1/16, ' num2str(H_quant_db(index_1_16)) ' dB)'], 'HorizontalAlignment', 'left');
% hold off;

% Comparison
figure;
plot(normalized_f_quant, H_quant_db, 'r-', 'DisplayName', 'Quantized'); % 红色实线
grid on;
hold on;
plot(normalized_f, H_db, 'b--', 'DisplayName', 'Original'); % 蓝色虚线
line([1/16, 1/16], ylim, 'Color', 'cyan', 'LineStyle', '-.', 'DisplayName', 'v = 1/16'); % 画出 x = 1/16 的黑色实线
line([1/8, 1/8], ylim, 'Color', 'cyan', 'LineStyle', '-', 'DisplayName', 'v = 1/8'); % 画出 x = 1/8 的黑色实线
title('Frequency Response of the Low-Pass Filter');
xlabel('Normalized Frequency f/f_s');
ylabel('Magnitude (dB)');
legend show; % 显示图例
hold off;

%%%Task4%%%
%%%%%%%%%%%

% 从均匀分布区间采样信号
x_min = -2^11;
x_max = 2^11;
x = x_min + (x_max - x_min).*rand(1, 1000000);

xl = conv(x, h);
P_xl = mean(xl.^2);

xl_quant = round(conv(x, h_quant));
xl_quant = conv(x, h_quant);
P_quant_noise = mean((xl-xl_quant).^2);

SNQR = 10*log10(P_xl/P_quant_noise)





