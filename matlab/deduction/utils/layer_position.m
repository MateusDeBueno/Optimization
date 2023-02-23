function [position] = layer_position(N,layer)
    
%     positions = zeros(layer,ceil(N/layer));
    
    m = ceil(N/layer);
    n = layer;
    
    positions_layer = zeros(1,m*n);
    positions_turn = zeros(1,m*n);
    
    k = 0;
    
    for y = 1:n
        for x = 1:m
            k = k + 1;
            positions_layer(k) = y;
            if (floor(y/2)==y/2) % code for even
                positions_turn(k) = x - 0.25;
            else % code for odd
                positions_turn(k) = x + 0.25;
            end
%             if (k>N)
%                 break
%             end
        end
    end
    
    positions_turn = positions_turn - sum(positions_turn)/(m*n);
    positions_layer = positions_layer - 1;
    position = [positions_layer' positions_turn'];
end

