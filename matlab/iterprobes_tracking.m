function M = iterprobes_tracking(starfield_dir, filename)

    #Beware there are a few limitations of the matlab scripts that I haven't found a good way to abstract away. If you would like to simulate a dgnf test, you must go into the file, either iterprobes.m or iterprobes_tracking.m, and change the 'dims' variable at the top of the script from [4 4 4 4] to [4 4 4]. This is because the dgnf doesn't have an obscuration polygon to read out of the file. The 'dims' variable gives the number of point on the polygons described in the file.

    #Another thing you must do if you would like to run a dgnf visual simulation is go into the file, iterprobes.m or iterprobes_tracking.m, and change the 'npolygons' variable from 13 to 12. This tells the script to display 12 polygons (3 for each probe) per frame, instead of 13 (3 for each probe plus the obscuration).

    dims = [4 4 4 4];
    probes = read_polygons(filename, dims);
    
    figure; hold on;
    ylim([-2000, 2000]);
    xlim([-2000, 2000]);
    axis equal;
    plot(0, 0);
    grid on;
    
    axes = [ 0,  1400, ...
            -1400,  0, ...
             0, -1400, ...
             1400,  0];
         
    circle(0, 0, 0.1 * 3600, 'k');
    circle(0, 0, -.167 * 3600, 'k');
        
    for i=1:4
        plot([0, 3600*axes(i*2-1)], [0, 3600*axes(i*2)], 'k');
    end
    
    darkgrey    = [0.662745, 0.662745, 0.662745];
    dimgray     = [0.411765, 0.411765, 0.411765];
    sgidarkgray = [0.333333, 0.333333, 0.333333];
        
    colors = [sgidarkgray; darkgrey; dimgray];
    
    npolygons = 13;
    nconfigurations = 60;
    probe_handles = double.empty(1, npolygons, 0);
    num_valid_configs = size(probes, 2) / (nconfigurations * npolygons);
    
    M = [];
    for i=1:num_valid_configs
        starfile = sprintf('%s/starfield%d.cat', starfield_dir, i);
        [starsx, starsy] = readstars(starfile);
        
        % leg = legend(starfile);
        % set(leg, 'FontSize', 16);
        
        for j=1:nconfigurations
            
            offset = ((i-1) * nconfigurations * npolygons) + (j-1)*npolygons;
            
            probe_handles = [];
            for k=1:npolygons
                xcoords = lfilter(@(x) ~isnan(x), probes(offset+k).xs);
                ycoords = lfilter(@(x) ~isnan(x), probes(offset+k).ys);
                probe_handles(k) = fill(xcoords, ycoords, colors(mod(k, size(colors,1))+1, :));
            end
            
            [rotatedx, rotatedy] = rotate2dcoord(starsx, starsy, (j-1) * (pi / 180));
            starplot = plot(rotatedx, rotatedy, 'b.', 'markersize', 12);
           
                    
            f = getframe(gcf);
            M = [M f f f f];
           
            pause(0.05);
        
            if j < nconfigurations
                for k=1:npolygons
                    delete(probe_handles(k));
                end 
                delete(starplot);
            elseif i < num_valid_configs
                pause(1.5);
                for k=1:npolygons
                    delete(probe_handles(k));
                end 
                delete(starplot);
            end
        end

                            
            f = getframe(gcf);
            M = [M f f f f f f f f f f f f f f f f f f f];
        
        pause(1);
    end
end
