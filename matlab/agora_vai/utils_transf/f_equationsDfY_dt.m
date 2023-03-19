function out = f_equationsDfY_dt(l)
    out = f_equationsDfY(l);
    C = num2cell(out);
    [~,~,Ip,Is,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = C{:};
    
    if (Ip>0)
        l.pr.phi = l.pr.phi - (l.pr.dt/2)*l.pr.fs*2*pi;
        out = f_equationsDfY(l);
    end
    if (Is>0)
        l.pr.phi = l.pr.phi + (l.pr.dt/2)*l.pr.fs/(2*pi);
        out = f_equationsDfY(l);
    end
end