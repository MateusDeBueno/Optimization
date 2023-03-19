function out = f_equationsDfD_dt(l)
    out = f_equationsDfD(l);
    C = num2cell(out);
    [~,~,Ip,Is,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = C{:};
    
    if (Ip>0)
        l.pr.phi = l.pr.phi - (l.pr.dt/2)*l.pr.fs*2*pi;
        out = f_equationsDfD(l);
    end
    if (Is>0)
        l.pr.phi = l.pr.phi + (l.pr.dt/2)*l.pr.fs/(2*pi);
        out = f_equationsDfD(l);
    end
end