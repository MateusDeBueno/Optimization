function [position] = horizontal_multiply(position,s_layer)
    position = position.*(([s_layer*ones(1,size(position,1))', ones(1,size(position,1))']));
end

