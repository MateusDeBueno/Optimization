%%
% Demo of using HISTCOMPRESS function

%%
% Compressing sine wave datapoints

x = 0:0.05:10;
y = sin(x);
[xc, yc] = histcompress(x,y,0.02);
plot(x,y,'b.',xc,yc,'ro-');
legend('original data','compressed data');
title(['Data compression rate: ' num2str(length(x)/length(xc))]);

%%
% Compressing sawtooth wave datapoints

y = sawtooth(x);
[xc, yc] = histcompress(x,y,0.01);
plot(x,y,'b.',xc,yc,'ro-');
legend('original data','compressed data');
title(['Data compression rate: ' num2str(length(x)/length(xc))]);

%%
% Compressing square wave datapoints

y = square(x);
[xc, yc] = histcompress(x,y,0.01);
plot(x,y,'b.',xc,yc,'ro-');
legend('original data','compressed data');
title(['Data compression rate: ' num2str(length(x)/length(xc))]);
