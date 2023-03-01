function [position] = layer_position2(N,max_spires_by_layer)
    
%     positions = zeros(layer,ceil(N/layer));
    layer = ceil(N/max_spires_by_layer);

    m = max_spires_by_layer;
    n = layer;
    
    positions_layer = zeros(1,N);
    positions_turn = zeros(1,N);
    
    k = 0;
    
    for y = 1:n %layer
        for x = 1:m %spiral
            k = k + 1;
            if (k>N)
                break
            end
            positions_layer(k) = y;
            if (floor(y/2)==y/2) % code for even
                positions_turn(k) = x - 0.25;
            else % code for odd
                positions_turn(k) = x + 0.25;
            end
        end
    end
    
%     positions_turn = positions_turn - sum(positions_turn)/(m*n);
    positions_turn = positions_turn -(max(positions_turn)+min(positions_turn))/2;
    positions_layer = positions_layer - 1;
    position = [positions_layer' positions_turn'];
end

