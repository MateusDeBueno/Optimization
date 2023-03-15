function Fr = f_Fr_litz(h,Nl,dl,strand,n,l)
    u0 = 4*pi*1e-7;
    cu = 1.72e-8; % [ohm*m]

    delta_w = sqrt(cu/(pi*h*l.pr.fs*u0));
%     n = 0.8;
    A = (pi/4)^(.75)*dl/delta_w*sqrt(n); %as vezes chamado de Delta
    Fs = (sinh(2*A)+sin(2*A))/(cosh(2*A)-cos(2*A));
    Fp = (2*(Nl*Nl*strand-1)/3)*(sinh(A)-sin(A))/(cosh(A)+cos(A));
    
    Fr = A*(Fs + Fp)*2;
end