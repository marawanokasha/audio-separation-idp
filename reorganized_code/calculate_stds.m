function [ stds1, stds2 ] = calculate_stds( scatts_diretory )
    % Compute standard deviations for all speakers

    allScattFiles = dir(scatts_diretory);
    X1All = [];
    X2All = [];

    for i=1:length(allScattFiles) % male or female level
        % make sure it's a directory but not the . or .. directory
        if allScattFiles(i).isdir == false && strncmp(allScattFiles(i).name,'.',1) == 0
            
            load(strcat(scatts_diretory, allScattFiles(i).name));
            X1All = [X1All X1];
            X2All = [X2All X2];
        end
    end
            
    % get the stds of the combined X1 of all speakers and combined X2 of all
    % speakers
    
    stds1 = std(X1All,0,2);
    stds2 = std(X2All,0,2);

end

