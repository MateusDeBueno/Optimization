function [out] = f_equations_plot(l,trafoo)
    if (trafoo == "YY")
        out = fplot_equationsYY(l);
    elseif (trafoo == "YD")
        out = fplot_equationsYD(l);
    elseif (trafoo == "DfD")
        out = fplot_equationsDfD(l);
    elseif (trafoo == "DiD")
        out = fplot_equationsDiD(l);
    elseif (trafoo == "DiY")
        out = fplot_equationsDiY(l);
    elseif (trafoo == "DfY")
        out = fplot_equationsDfY(l);
    end
end