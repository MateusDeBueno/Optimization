function DrawClosedPolygon(p)
    for i = 1:size(p,1)
        mi_addnode(p(i,1),p(i,2));
    end

    for i = 1:size(p,1)-1
        mi_addsegment(p(i,1),p(i,2),p(i+1,1),p(i+1,2));
    end
    mi_addsegment(p(end,1),p(end,2),p(1,1),p(1,2));
end