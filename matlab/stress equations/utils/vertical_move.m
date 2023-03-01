function [p] = vertical_move(p,x)
    p = p + (ones(length(p),2).*[0 x])';
end