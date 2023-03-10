function Rm = f_Rm_air_gap(d,w_x,w_y,R_x,R_y,u0)

    sigma_x = R_x/(d/(u0*w_y));
    sigma_y = R_y/(d/(u0*w_x));
    sigma = sigma_y*sigma_x;

    Rm = sigma*d/(u0*w_x*w_y);

end
