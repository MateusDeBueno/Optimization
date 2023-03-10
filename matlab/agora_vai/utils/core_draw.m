function core_draw(x2,x3,x4,y1,y2,y3)
    p_core_u_right = [x2,y1
                    x2,y2
                    x3,y2
                    x3,y1
                    x4,y1
                    x4,y3];
    p_core_u = [p_core_u_right; flip(horizontal_flip(p_core_u_right))];
    DrawClosedPolygon(p_core_u);
    P_core_u = [0,y2];
    mi_addblocklabel(P_core_u);
    mi_selectlabel(P_core_u);
    mi_setblockprop('Ferrite',0,0,'',0,0,0);
%     all_center = [all_center; P_core_u];
    mi_clearselected()

    p_core_d = vertical_flip(p_core_u);
    DrawClosedPolygon(p_core_d);
    P_core_d = [0,-y2];
    mi_addblocklabel(P_core_d);
    mi_selectlabel(P_core_d);
    mi_setblockprop('Ferrite',0,0,'',0,0,0);
%     all_center = [all_center; P_core_d];
    mi_clearselected()
end

