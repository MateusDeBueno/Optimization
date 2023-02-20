function Rbasic = f_Rbasic(w,d,h,u0)
    % https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9332553
    
    Rbasic = 1/(u0*(w/d + 2/pi*(1+log(pi*h/(2*d)))));
end

