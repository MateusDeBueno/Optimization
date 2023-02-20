function [position_litz,position,group,max_spires_by_layer] = layer_position_litz2(d_litz,strands,N,sbw,G_carretel,x2,s)

    %deixa proporcional o tamanho do envelope em relacao ao tamanho de cada
    %radio dos strings
    
    
    
    
    
    
    if strands==1
        
        d_target = d_litz*sbw; %diametro do envelope
    %     s_turn = d_target*1.1;
        max_spires_by_layer = floor(G_carretel/d_target); %numero maximo de espiras
        s_layer = d_target*0.88;
        s_turn = d_target;
    %     s_layer = d_target;
    
        position = layer_position2(N,max_spires_by_layer);
        position = horizontal_multiply(position,s_layer);
        position = vertical_multiply(position,s_turn);
        position = horizontal_move(position,x2+s+d_target/2);
        
        
        group = 1:N;
        position_litz = position;
    else
        r_litz = d_litz/2;
        [position_strand2,radii] = read_litz_data(strands);
        k_litz_radii = r_litz/radii;
        position_strand = position_strand2*k_litz_radii*sbw; %posicao individual de cada string
        d_target = 2*k_litz_radii*sbw; %diametro do envelope
        max_spires_by_layer = floor(G_carretel/d_target); %numero maximo de espiras
    %     s_turn = d_target*1.1;
        s_layer = d_target*0.88;
        s_turn = d_target;
    %     s_layer = d_target;
    
        position = layer_position2(N,max_spires_by_layer);
        position = horizontal_multiply(position,s_layer);
        position = vertical_multiply(position,s_turn);
        position = horizontal_move(position,x2+s+d_target/2);

        all_position_x = zeros(1,strands*N);
        all_position_y = zeros(1,strands*N);
        group = zeros(1,strands*N);

        if (strands==2 || strands==4)
            theta = 90;
            R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
            position_strand = (R*[position_strand(:,1),position_strand(:,2)]')'; %rotation 90 degree
        end

        j = 0;
        for k=1:N

            if (strands ~= 2 && strands ~= 4)
                theta = rand(1)*360; % to rotate 90 counterclockwise
                R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
                position_strand = (R*[position_strand(:,1),position_strand(:,2)]')'; %randon rotation
            end

            for str=1:strands
                j = j + 1;
                all_position_x(j) = position_strand(str,1) + position(k,1);
                all_position_y(j) = position_strand(str,2) + position(k,2);
                group(j) = k;
            end
        end
        
        group = group';
        position_litz = [all_position_x' all_position_y'];
    end
end

