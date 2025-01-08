% 圆极化波参数
E0 = 1;            % 振幅
omega = 2*pi;      % 角频率

% 时间范围
t = linspace(0, 1, 1000);  % 从0到1秒，取1000个时间点

% 左旋分量
E_left = E0 * exp(1j * omega * t);

% 右旋分量
E_right = E0 * exp(1j *(omega * t+pi/2));

% 圆极化波
E_circle = E_left + E_right;

% 计算功率
power_left = E_left*E_left'/length(t);
power_right = E_right*E_right'/length(t);
total_power = E_circle*E_circle'/length(t);

% 绘制波形
figure;
subplot(2, 1, 1);
plot(t, real(E_circle), 'b', t, imag(E_circle), 'r');
title('圆极化波波形');
xlabel('时间');
ylabel('电场强度');
legend('实部', '虚部');


