function air_draw(x2,x4)
    P_air = [x2+x4,0];
    mi_addblocklabel(P_air);
    mi_selectlabel(P_air);
    mi_setblockprop('Air',0,0,'',0,0,0);
%     all_center = [all_center; P_air];
    mi_clearselected();
end

