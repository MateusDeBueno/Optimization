function [legen_name] = create_legend_contourf(contourLevels, colors)
% function [legen_name] = create_legend_contourf(contourLevels)

    clear legen_name;
    
    n_contour = numel(contourLevels);
    legen_name(1) = append(' $<$ ',string(contourLevels(1)));
    legen_name(2:n_contour) = append(string(contourLevels(1:n_contour-1)),' to ',string(contourLevels(2:n_contour)));
    legen_name(end+1) = append(' $>$ ',string(contourLevels(end)));
    
    for ii=1:n_contour+1
        fill([1 1], [1 1],colors(ii,:),'Edgecolor', 'none')
    end
end