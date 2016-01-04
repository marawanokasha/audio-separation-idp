function [ Dnmf1, Dnmf2 ] = get_speaker_dictionaries( speakerName, dictionaryDirectory )
% Gets the dictionaries (1st and 2nd level) of a specific speaker

    allDictionaryFiles = dir(dictionaryDirectory);

    % Create a dictionary from the scattering representation of each speaker
    for i=1:length(allDictionaryFiles)
        % make sure it's a directory but not the . or .. directory
        if allDictionaryFiles(i).isdir == false && strncmp(allDictionaryFiles(i).name,'.',1) == 0

            fileName = allDictionaryFiles(i).name;
            if strfind(fileName, speakerName) >= 0

                load(strcat(dictionaryDirectory, fileName));
                Dnmf1 = Dnmf1;
                Dnmf2 = Dnmf2;
            end
        end
    end

end