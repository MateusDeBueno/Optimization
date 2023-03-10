function gap = f_get_gap(core_geometry,N,ur,L_goal)

    k = 0;
    for gap=1e-2:1e-4:5

        L = f_get_indutance(core_geometry,N,gap*1e-3,ur);
        if (real(L)>0 && imag(L)==0)
            k = k + 1;
            vector_L(k) = L;
            vector_gap(k) = gap;
        end
    end

    [ d, ix ] = min( abs( vector_L-L_goal ) );

    gap = vector_gap(ix);
end

