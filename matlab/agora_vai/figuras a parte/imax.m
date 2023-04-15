clear; close all; clc;

% https://www.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end


Vbat = [400 600 400 396 405 800 360 450 350 360 400 320 360];
maxP = [300 300 250 150 150 270 250 150 50 100 50 100 50]*1e3;

Imax = maxP./Vbat;


figure
scatter(Imax,Vbat, 'filled')
xlabel('$I_{max}$[A]')
ylabel('$V_{bat}$[V]')
set(gca, 'FontSize', 20)
grid on
grid minor