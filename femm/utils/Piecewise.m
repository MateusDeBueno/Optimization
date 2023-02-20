% function [answ] = Piecewise(cond1,ans1,cond2,ans2)
%     if (cond1)
%         answ = ans1;
%     elseif (cond2)
%         answ = ans2;
%     end
% end

function [answ] = Piecewise(varargin)
    switch nargin
        case 4 %quatro entradas
            cond1 = varargin{1};
            ans1 = varargin{2};
            cond2 = varargin{3};
            ans2 = varargin{4};

            if (cond1)
                answ = ans1;
            elseif (cond2)
                answ = ans2;
            end
        case 3 %3 entradas
            cond1 = varargin{1};
            ans1 = varargin{2};
            ans2 = varargin{3};

            if (cond1)
                answ = ans1;
            else
                answ = ans2;
            end
        otherwise   
    end
end
