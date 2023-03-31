function [xc, yc] = histcompress(x, y, eps)
% HISTCOMPRESS Compress two-dimensional data by removing redundant points.
%   [XC, YC] = HISTCOMPRESS(X, Y, EPS) removes redundant data points within
%   EPS tolerance interval using GE Historian Compression algorithm
%   (similar to Swinging Door algorithm) described here:
%   http://www.evsystems.net/files/GE_Historian_Compression_Overview.ppt

% Copyright (c) 2010, 2012 Yuriy Skalko <yuriy.skalko@gmail.com>

  N = length(x);
  xc = zeros(N, 1);
  yc = zeros(N, 1);

  j = 0; i = 2;
  slope_high = 0;
  slope_low = 0;
  slope = 0;
  while i < N
    if slope <= slope_low || slope >= slope_high % point is outside of the angle
      j = j + 1;
      xc(j) = x(i-1); % save previous point inside the angle
      yc(j) = y(i-1);
      slope_low = -Inf; % reset the angle
      slope_high = +Inf;
    end
    slope_high = min(slope_high,(y(i)+eps-yc(j))/(x(i)-xc(j))); % reduce the angle
    slope_low = max(slope_low,(y(i)-eps-yc(j))/(x(i)-xc(j)));
    i = i + 1;
    slope = (y(i)-yc(j))/(x(i)-xc(j)); % slope for next point
  end

  xc(j+1) = x(N); % save last point
  yc(j+1) = y(N);

  xc(j+2:end) = []; % delete the rest
  yc(j+2:end) = [];

end
