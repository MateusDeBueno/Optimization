clear; close all; clc;

x = 0:0.01:pi;
y = 5*sin(x).*cos(10*x) + rand(size(x));
figure; hold on;
plot(x, y);
xlabel('X');
ylabel('Y');

x_r = 0.95; y_r = -3.3; w_r = 0.2; h_r = 0.7;
rectangle('Position', [x_r-w_r/2, y_r-h_r/2, w_r, h_r], ...
'EdgeColor', [0.4, 0.1, 0.4], 'LineWidth',2);

x_a = 0.58; y_a = 0.18; w_a = 0.3; h_a = 0.3;
ax = axes('Units', 'Normalized', ...
'Position', [x_a, y_a, w_a, h_a], ...
'XTick', [], ...
'YTick', [], ...
'Box', 'on', ...
'LineWidth', 2, ...
'Color', [0.95, 0.99, 0.95]);
hold on;

plot(x, y);
xlabel('Detail at X==0.95');
axis([x_r-w_r/2, x_r+w_r/2, y_r-h_r/2, y_r+h_r/2]);
