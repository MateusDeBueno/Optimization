function [Fr] = get_Fr(f,d_str,k,Nl,p)
    % Winding Resistance and Power Loss for Inductors With Litz and Solid-Round Wires
    
    global cu;
    global u0;
    
    qsi_w = sqrt(cu/(pi*f*u0)); % eq 8
    A = (pi/4)^(.75)*d_str/qsi_w*sqrt(d_str/p*k); % eq 13
    Fr = A*((sinh(2*A)+sin(2*A))/(cosh(2*A)-cos(2*A))) ...
                        + A*((2*(Nl*Nl*k-1)/3)*(sinh(A)-sin(A))/(cosh(A)+cos(A))); % eq 2
end

