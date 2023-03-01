function DrawClosedCirc(p,dm)
    point_up = [p(1,1), p(1,2)+dm/2];
    point_down = [p(1,1), p(1,2)-dm/2];

    mi_addnode(point_up);
    mi_addnode(point_down);

    mi_addarc(point_up(1,1),point_up(1,2),point_down(1,1),point_down(1,2),180,2)
    mi_addarc(point_down(1,1),point_down(1,2),point_up(1,1),point_up(1,2),180,2)
end