function fill_between_functions(f1, f2, fill_color, x_limits)
% fill_between_functions(f1, f2, fill_color, x_limits)
% Fills the area between two symbolic functions with the specified color.
%
% Inputs:
%   f1 - symbolic function 1
%   f2 - symbolic function 2
%   fill_color - color to fill area with, specified as a string or RGB triplet
%   x_limits - optional 2-element vector specifying the x-axis limits
%
% Example usage:
%   syms x
%   f1 = x;
%   f2 = x^2;
%   fill_between_functions(f1, f2, 'g', [-5, 5])

% Default x limits
if nargin < 4
    x_limits = [-5, 5];
end

% Evaluate functions at x points
x_points = linspace(x_limits(1), x_limits(2));
y1_points = double(subs(f1, x_points));
y2_points = double(subs(f2, x_points));

% Plot functions and fill between
fill([x_points fliplr(x_points)], [y1_points fliplr(y2_points)], fill_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
hold on
% plot(x_points, y1_points, fill_color, 'LineWidth', 0.2)
% plot(x_points, y2_points, fill_color, 'LineWidth', 0.2)
xlim(x_limits)

end