function [ coords ] = nodes2calibrationpoints(corner,filename,num,distrib,type,minelem,Zmax,Zmin)
% Part of FAST Calibration v2.4 - GPLv3
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% DOI:      http://doi.org/10.5880/wsm.2023.002
% Manual:  	http://doi.org/10.48440/wsm.2023.002
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates the coordinates for calibration points based on
% the nodes provided in the input file (filename). The coordinate system in
% the local and regional model have to be identical.
%
% Corners: X and Y coordinates of the corners of the local model.
% filename: Name of the file that contains the nodes.
% Num: Number of desired calibration points.
% distrib: This is a factor for the distribution of calibration points. The
% larger the number the more evenly distributed are the calibration points.
% This variable only applies to type corner. It is ignored (but
% has to be present) for type random and border. A value between 10 and 30
% is recommended.
% type: 'border', corner', or 'random'.
% minelem: The distance between the border of the model and the closest
% calibration point (only for type 'border').
%

if (isequal(type,'corner')) || (isequal(type,'border')) || (isequal(type,'random'))
    % Reduce the regional nodes to the aproximate area of the local model 
    Xmax = max(corner(:,1));
    Xmin = min(corner(:,1));

    Ymax = max(corner(:,2));
    Ymin = min(corner(:,2));

    data = csvread(filename);
    data(:,1) = [];

    i =  (data(:,1) > Xmin) & (data(:,1) < Xmax)...
        & (data(:,2) > Ymin) & (data(:,2) < Ymax) & (data(:,3) < Zmax);

    nodes_reg = data(i,:);

    % Find datapoints in polygon
    polyx = [ corner(:,1) ; corner(1,1) ];
    polyy = [ corner(:,2) ; corner(1,2) ];
    polygon = [ polyx polyy];

    % Compute new corner points outside of zone prone to boundary effects
    poly_small = polygon;
    for i = 1:4
        % Compute change in corner nodes
        diag = sqrt( 2 * minelem^2 );
        x = polygon(i,1) - polygon((i+1),1);
        y = polygon(i,2) - polygon((i+1),2);
        h = sqrt( x^2 + y^2 );
        beta = atand(x/y) + 45;
        dx = abs(diag * sind(beta));
        dy = abs(diag * cosd(beta));

        % Compute new corners
        if (x ~= 0) && (y ~= 0)
            % X Values:
            if x < 0
                poly_small(i,1) = polygon(i,1) + dx;
            elseif x > 0
                poly_small(i,1) = polygon(i,1) - dx;
            end
            %Y Values:
            if y < 0
                poly_small(i,2) = polygon(i,2) + dy;
            elseif y > 0
                poly_small(i,2) = polygon(i,2) - dy;
            end
        else
            if y == 0 && x < 0
                poly_small(i,1) = polygon(i,1) + minelem;
                poly_small(i,2) = polygon(i,2) - minelem;
            elseif x == 0 && y > 0
                poly_small(i,1) = polygon(i,1) - minelem;
                poly_small(i,2) = polygon(i,2) - minelem;
            elseif y == 0 && x > 0
                poly_small(i,1) = polygon(i,1) - minelem;
                poly_small(i,2) = polygon(i,2) + minelem;
            elseif x == 0 && y < 0
                poly_small(i,1) = polygon(i,1) + minelem;
                poly_small(i,2) = polygon(i,2) + minelem;
            end
        end
    end
    poly_small(5,:) = poly_small(1,:);

    [in , ~ ] = inpolygon(nodes_reg(:,1),nodes_reg(:,2),poly_small(:,1),poly_small(:,2));
    inside_nodes = nodes_reg(in,:);
    i = inside_nodes(:,3) < ( Zmin + 50);
    inside_nodes(i,:) = [];

    if num > ( length(inside_nodes) * 0.8 )
        disp('Too many calibration points - select a smaller ''num''')
        return
    end

    clear in Xmax Xmin Ymax Ymin Zmax Zmin beta diag dx dy h data i nodes_reg
    clear polyx polyy x y

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isequal(type,'corner')
        % Corner
        coords = zeros(1,3);
        for i = 1:4
            dist = bsxfun(@hypot, inside_nodes(:,1) - corner(i,1), inside_nodes(:,2) - corner(i,2));
            [~,list] = sort(dist);
            prel = inside_nodes(list(1:((num/4)*distrib)),:);
            prel = sortrows(prel,3);
            dist = zeros(1,(length(prel)-1));
            for k = 2:length(prel)
                dist(k-1) = prel((k-1),3) - prel(k,3);
            end
            [~,list] = sort(dist,'descend');    
            coords = [ coords; prel(list(1:(num/4)),:)];   
        end

    coords(1,:) = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif isequal(type,'border')
    % Border
    num4 = round(num/4);
    pd = makedist('Normal','mu',0.5,'sigma',0.4);
    coords = zeros(round(num*1.5),3);
    c = 1;

    for i = 1:4
        % Create random normal distribution between 0 and 1
        bp = random(pd,1,round(num4*1.05));
        bp(bp > 1 | bp < 0) = [];
        % Compute border length
        x = polygon(i,1) - polygon((i+1),1);
        y = polygon(i,2) - polygon((i+1),2);
        h = sqrt( x^2 + y^2 );
        alpha = atand(x/y);
        % Compute positions along border
        h_ = bp * h;
        if x < 0 && y > 0
            x_ = -h_ .* sind(alpha) + corner(i,1);
            y_ = -h_ .* cosd(alpha) + corner(i,2);
        elseif x > 0 && y > 0
            x_ = -h_ .* sind(alpha) + corner(i,1);
            y_ = -h_ .* cosd(alpha) + corner(i,2);
        elseif x > 0 && y < 0
            x_ = h_ .* sind(alpha) + corner(i,1);
            y_ = h_ .* cosd(alpha) + corner(i,2);
        elseif x < 0 && y < 0
            x_ = h_ .* sind(alpha) + corner(i,1);
            y_ = h_ .* cosd(alpha) + corner(i,2);

        elseif x == 0 && y > 0
            x_ = ones(length(h_)) .* corner(i,1);
            y_ = -h_ + corner(i,2);
        elseif x == 0 && y < 0
            x_ = ones(length(h_)) .* corner(i,1);
            y_ = h_ + corner(i,2);
        elseif x > 0 && y == 0
            x_ = -h_ + corner(i,1);
            y_ = ones(length(h_)) .* corner(i,2);
        elseif x < 0 && y == 0
            x_ = h_ + corner(i,1);
            y_ = ones(length(h_)) .* corner(i,2);
        end
        % Find closest nodes
        for j = 1:length(h_)
            distance = sqrt( (inside_nodes(:,1) - x_(j)).^2 + (inside_nodes(:,2) - y_(j)).^2 );
            k = find(distance == min(distance));
            if length(k) > 1
                cali = round(rand(1)*(length(k) - 1)) + 1;
            else
                cali = 1;
            end
            coords(c,:) = inside_nodes(k(cali),:);
            c = c + 1;
            inside_nodes(k(cali),:) =[];
        end
    end
    % Remove Zero values and surplus calibration points
    coords(c:length(coords),:) = [];
    k = round(rand((length(coords) - num),1) * length(coords));
    coords(k,:) = [];

    clear k j i h_ h distance bp alpha c num4 pd x y x_ y_

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif isequal(type,'random')
        % Random
        i = rand(1,num);
        i = round(i*length(inside_nodes));
        coords = inside_nodes(i,:);
    end
    
    csvwrite('branch_calibration_nodes.csv',coords);

elseif isequal(type,'user')
    disp(['User defined calibration nodes are used.'])
    coords = csvread('user_defined_branch_calibration_nodes.csv');
    
    % Find datapoints in polygon
    polyx = [ corner(:,1) ; corner(1,1) ];
    polyy = [ corner(:,2) ; corner(1,2) ];
    polygon = [ polyx polyy];
    
else
    disp(['Type of calibration procedure not specified/unknown.']);

end

scatter3(coords(:,1),coords(:,2),coords(:,3));
hold on
scatter3(corner(:,1),corner(:,2),[0,0,0,0],50,'black','LineWidth',1.5);
plot3(polygon(:,1),polygon(:,2),[0,0,0,0,0],'k','LineWidth',1.5)
hold off

end