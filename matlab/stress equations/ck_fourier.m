function c_k_rms = ck_fourier(ts,derivadas,pts_inics)
    syms t k fs real positive
    Ts = 1/fs;
    c_k = 0;
    for i=1:length(ts)-1
        c_k = c_k + int(((derivadas(i)*(t-ts(i)) + pts_inics(i))*exp(1i*2*pi*k*t/Ts)), [ts(i) ts(i+1)]);
    end
    c_k = simplify(abs(2*c_k/Ts));
    c_k_rms = c_k/sqrt(2);
end