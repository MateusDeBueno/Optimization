function [position_litz,group] = layer_position_litz(d_w,d_litz,position,fill)
    
    d_wire = d_w;
    d_wire_target = d_wire;
%     fill = 0.55;

    [n, x_vec, y_vec] = get_packing_pattern(d_litz, d_wire_target, fill);

    position_strand = [x_vec;y_vec]';

    get_packing_plot(d_wire,d_litz,x_vec,y_vec,fill,n)

    all_position_x = zeros(1,n*length(position));
    all_position_y = zeros(1,n*length(position));
    group = zeros(1,n*length(position));
    j = 0;
    for k=1:length(position)
        
        for strand=1:n
            j = j + 1;
            all_position_x(j) = position_strand(strand,1) + position(k,1);
            all_position_y(j) = position_strand(strand,2) + position(k,2);
            group(j) = k;
        end
    end
    
    group = group';
    position_litz = [all_position_x' all_position_y'];
end
