function d_target = f_get_d_target(strands,d_litz,sbw)
    r_litz = d_litz/2;
    [position_strand2,radii] = read_litz_data(strands);
    k_litz_radii = r_litz/radii;
    d_target = 2*k_litz_radii*sbw; %diametro do envelope
end

