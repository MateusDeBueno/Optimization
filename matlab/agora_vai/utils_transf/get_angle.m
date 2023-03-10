function [k] = get_angle(Pout_med,Ptarget)
    k = ceil(length(Pout_med)/2);
    Px = abs(Pout_med-Ptarget);
    Pdx = diff(Px);
    for i=1:100 %Least mean squares filter
        k = round(k - 1e-1*Px(k)/Pdx(k), 0);
    end
end