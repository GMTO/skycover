function gmacs_probes()

    starfield_dir = '../starfiles';
    filename = '../plot_gmacs.out';
    
    %Beware there are a few limitations of the matlab scripts that I haven't found a good way to abstract away. 
    % If you would like to simulate a dgnf test, you must go into the file, either iterprobes.m or iterprobes_tracking.m, 
    % and change the 'dims' variable at the top of the script from [4 4 4 4] to [4 4 4]. 
    % This is because the dgnf doesn't have an obscuration polygon to read out of the file. 
    % The 'dims' variable gives the number of point on the polygons described in the file.

    %Another thing you must do if you would like to run a dgnf visual simulation is go into the file, 
    % iterprobes.m or iterprobes_tracking.m, and change the 'npolygons' variable from 13 to 12. 
    % This tells the script to display 12 polygons (3 for each probe) per frame, 
    % instead of 13 (3 for each probe plus the obscuration).

    dims = [4 183 4];
    probes = read_polygons(filename, dims);
    
    figure; hold on;
    ylim([-2000, 2000]);
    xlim([-2000, 2000]);
    axis equal;
    axis square;
    xlabel('mm')
    ylabel('mm')

    axes = [ 0,  1400, ...
            -1400,  0, ...
             0, -1400, ...
             1400,  0];
        
    maxradius = 65.02345*7 + 2.49e-3*7^3;
    %circle(0, 0, 0.1 * 3600, 'k');
    circle(0, 0, maxradius, 'k');
        
    for i=1:4
        plot([0, axes(i*2-1)], [0, axes(i*2)], 'k');
    end

    darkgrey    = [0.662745, 0.662745, 0.662745];
    dimgray     = [0.411765, 0.411765, 0.411765];
    sgidarkgray = [0.333333, 0.333333, 0.333333];
    red = [1.0, 0.0, 0.0];
    grn = [0.0, 1.0, 0.0];
    blu = [0.0, 0.0, 1.0];
    
   colors = [sgidarkgray; dimgray; darkgrey];

    npolygons = 13;
    nconfigurations = size(probes, 2) / npolygons;
    probe_handles = double.empty(1, npolygons, 0);
    
    for i=1:nconfigurations
        starfile = sprintf('%s/starfield%d.cat', starfield_dir, i);
        fprintf('Starfile: %s\n',starfile)
        [starsx, starsy] = readstars(starfile);
        
        
        offset = (i-1)*npolygons;
        probe_handles = [];
        for j=1:npolygons
            xcoords = lfilter(@(x) ~isnan(x), probes(offset+j).xs);
            ycoords = lfilter(@(x) ~isnan(x), probes(offset+j).ys);
            probe_handles(j) = fill(xcoords, ycoords, colors(mod(j, size(colors,1))+1, :));
        end
        
        starplot = plot(starsx, starsy, 'b.', 'markersize', 12);
    
    pause(1.5);
        
        if i < nconfigurations
            for j=1:npolygons
                delete(probe_handles(j));
            end
            delete(starplot);
        end
    end
end
