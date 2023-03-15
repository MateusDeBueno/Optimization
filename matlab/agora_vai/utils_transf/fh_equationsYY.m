function out = fh_equationsYY(l,nn)
    if (-pi <= l.pr.phi && l.pr.phi <= -2*pi/3)
        out = l.eq.YY.eq1.fh(l.pr.L1,l.pr.L2,l.pr.Ldab,l.pr.M,l.pr.Vi,l.pr.d,l.pr.fs,l.pr.phi,nn);
    elseif (-2*pi/3 <= l.pr.phi && l.pr.phi <= -pi/3)
        out = l.eq.YY.eq2.fh(l.pr.L1,l.pr.L2,l.pr.Ldab,l.pr.M,l.pr.Vi,l.pr.d,l.pr.fs,l.pr.phi,nn);
    elseif (-pi/3 <= l.pr.phi && l.pr.phi <= 0)
        out = l.eq.YY.eq3.fh(l.pr.L1,l.pr.L2,l.pr.Ldab,l.pr.M,l.pr.Vi,l.pr.d,l.pr.fs,l.pr.phi,nn);
    elseif (0 <= l.pr.phi && l.pr.phi <= pi/3)
        out = l.eq.YY.eq4.fh(l.pr.L1,l.pr.L2,l.pr.Ldab,l.pr.M,l.pr.Vi,l.pr.d,l.pr.fs,l.pr.phi,nn);
    elseif (pi/3 <= l.pr.phi && l.pr.phi <= 2*pi/3)
        out = l.eq.YY.eq5.fh(l.pr.L1,l.pr.L2,l.pr.Ldab,l.pr.M,l.pr.Vi,l.pr.d,l.pr.fs,l.pr.phi,nn);
    elseif (2*pi/3 <= l.pr.phi && l.pr.phi <= pi)
        out = l.eq.YY.eq6.fh(l.pr.L1,l.pr.L2,l.pr.Ldab,l.pr.M,l.pr.Vi,l.pr.d,l.pr.fs,l.pr.phi,nn);
    end
end