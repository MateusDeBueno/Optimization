% function put_wires(wires_r,d_w,turn,n,name)
function put_wires(wires_r,d_w,turn,group,name)
    for k=1:length(wires_r)
        point = wires_r(k,:);
        DrawClosedCirc(point,d_w)
        mi_addblocklabel(point);
        mi_selectlabel(point);
        mi_setblockprop('Copper',0,0,strcat('Coil_',num2str(group(k)),name),0,0,turn);
        mi_clearselected()
    end
end