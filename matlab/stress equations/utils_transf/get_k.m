function [k] = get_k(vP,Ptarget)
    k = ceil(length(vP)/2);
    Px = abs(vP-Ptarget);
    Pdx = diff(Px);
    for i=1:100 %Least mean squares filter
        k = round(k - 1e-1*Px(k)/Pdx(k), 0);
    end
end