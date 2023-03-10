function cmap = f_create_cmap(n_unique_color, left_color, right_color)

    cols = n_unique_color;
    cmap = interp1([0, 1], [left_color; right_color], linspace(0, 1, cols));

end