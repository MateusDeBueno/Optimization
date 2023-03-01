function Fr = f_Fr_litz2(f,Nl)
    u0 = 4*pi*1e-7;
    cu = 1.72e-8; % [ohm*m]

    %o fio litz la de baixo
    k = 180;
    d_l = 101.6e-6;

    delta_w = sqrt(cu/(pi*f*u0));
    n = 0.8;
    A = (pi/4)^(.75)*d_l/delta_w*sqrt(n); %as vezes chamado de Delta
    Fs = (sinh(2*A)+sin(2*A))/(cosh(2*A)-cos(2*A));
    Fp = (2*(Nl*Nl*k-1)/3)*(sinh(A)-sin(A))/(cosh(A)+cos(A));
    
    Fr = A*(Fs + Fp)*2;
end