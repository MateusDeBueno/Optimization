function [position] = vertical_multiply(position,s_turn)
%     position = position.*(([ones(1,size(position,2))', s_turn*ones(1,size(position,2))']))';
    position = position.*(([ones(1,size(position,1))', s_turn*ones(1,size(position,1))']));

end

