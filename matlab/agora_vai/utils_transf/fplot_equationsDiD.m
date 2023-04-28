function out = fplot_equationsDiD(l)
    if (-pi <= l.pr.phi && l.pr.phi <= -2*pi/3)
        out = l.eq.DiD.eq1.f_plot(l.pr.L1,l.pr.L2,l.pr.Ldab,l.pr.M,l.pr.Vi,l.pr.d,l.pr.fs,l.pr.phi);
    elseif (-2*pi/3 <= l.pr.phi && l.pr.phi <= -pi/3)
        out = l.eq.DiD.eq2.f_plot(l.pr.L1,l.pr.L2,l.pr.Ldab,l.pr.M,l.pr.Vi,l.pr.d,l.pr.fs,l.pr.phi);
    elseif (-pi/3 <= l.pr.phi && l.pr.phi <= 0)
        out = l.eq.DiD.eq3.f_plot(l.pr.L1,l.pr.L2,l.pr.Ldab,l.pr.M,l.pr.Vi,l.pr.d,l.pr.fs,l.pr.phi);
    elseif (0 <= l.pr.phi && l.pr.phi <= pi/3)
        out = l.eq.DiD.eq4.f_plot(l.pr.L1,l.pr.L2,l.pr.Ldab,l.pr.M,l.pr.Vi,l.pr.d,l.pr.fs,l.pr.phi);
    elseif (pi/3 <= l.pr.phi && l.pr.phi <= 2*pi/3)
        out = l.eq.DiD.eq5.f_plot(l.pr.L1,l.pr.L2,l.pr.Ldab,l.pr.M,l.pr.Vi,l.pr.d,l.pr.fs,l.pr.phi);
    elseif (2*pi/3 <= l.pr.phi && l.pr.phi <= pi)
        out = l.eq.DiD.eq6.f_plot(l.pr.L1,l.pr.L2,l.pr.Ldab,l.pr.M,l.pr.Vi,l.pr.d,l.pr.fs,l.pr.phi);
    end
end