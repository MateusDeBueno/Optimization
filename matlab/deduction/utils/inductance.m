function y=inductance(n,ri,ro,z)
    openfemm;
    newdocument(0);
    mi_probdef(0,'inches','axi',1e-8,0,30);
    mi_drawrectangle(ri,-z/2,ro,z/2);
    r=2*max([ro,ri,z]);
    mi_drawarc(0,-r,0,r,180,5);
    mi_drawline(0,-r,0,r);
    mi_addcircprop('icoil',1,1);
    mi_addblocklabel((ri+ro)/2,0);
    mi_addblocklabel(0.75*r,0);
    mi_addmaterial('coil',1,1,0,0,0,0,0,1,0,0,0);
    mi_addmaterial('air' ,1,1,0,0,0,0,0,1,0,0,0);
    mi_addboundprop('abc',0,0,0,0,0,0,1/(r*0.0254*pi*4.e-7),0,2);
    mi_selectlabel((ri+ro)/2,0);
    mi_setblockprop('coil',0,r/20,'icoil',0,0,n);
    mi_clearselected;
    mi_selectlabel(0.75*r,0);
    mi_setblockprop('air',0,r/100,'<None>',0,0,0);
    mi_clearselected;
    mi_selectarcsegment(r,0);
    mi_setarcsegmentprop(5,'abc',0,0);
    mi_saveas('c:\\femm42\\examples\\tmp.fem');
    mi_analyze;
    mi_loadsolution;
    c=mo_getcircuitproperties('icoil');
    y=c(3);
%     closefemm;