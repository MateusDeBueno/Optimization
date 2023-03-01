function get_packing_plot(d_wire,d_litz,x_vec,y_vec,fill,n)
% Plot the geometry of a litz wire.
%
%    Plot the following parameters:
%        - the discrete strands
%        - the equivalent solid wire
%
%    Parameters:
%        wire (struct): struct with the litz wire parameters
%
%    (c) 2016-2020, ETH Zurich, Power Electronic Systems Laboratory, T. Guillod

% plot init
figure('name', 'packing')

% angle vector
phi_vec = linspace(0, 2*pi);

% equivalent solid wire
x_wire = (d_wire./2).*cos(phi_vec);
y_wire = (d_wire./2).*sin(phi_vec);

% discrete strand
x_litz = (d_litz./2).*cos(phi_vec);
y_litz = (d_litz./2).*sin(phi_vec);

% plot the geometry
subplot(1,1,1)
plot(x_wire, y_wire, 'b')
hold('on');
for i=1:n
    x_litz_tmp = x_litz+x_vec(i);
    y_litz_tmp = y_litz+y_vec(i);
    plot(x_litz_tmp, y_litz_tmp, 'r')
end
axis('equal');
xlabel('x [mm]')
ylabel('y [mm]')
title('wire geometry', 'interpreter', 'none')

% plot title
d_wire_str = sprintf('d_wire = %.2f mm', d_wire);
d_litz_str = sprintf('d_litz = %.2f mm', d_litz);
fill_str = sprintf('fill = %.2f %%', 100.*fill);
n_str = sprintf('n = %d', n);
str = [d_wire_str ' / ' d_litz_str ' / ' fill_str ' / ' n_str];
sgtitle(str, 'interpreter', 'none')

end