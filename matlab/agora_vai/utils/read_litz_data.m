function [position_strand,radii] = read_litz_data(strands)
    % http://hydra.nat.uni-magdeburg.de/packing/cci/
    rad = readtable('radius.txt');

    radii = rad.Var2(strands); %get radius
    pos = readtable(strcat('cci',num2str(strands),'.txt')); %get positions
    for k=1:strands
        positions_x(k) = pos.Var2(k);
        positions_y(k) = pos.Var3(k);
    end
    position_strand = [positions_x;positions_y]';

end

