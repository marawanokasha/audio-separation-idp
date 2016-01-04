function speakers = load_speakers(dataDirectory, trainingCount, testCount, resample_freq, Npad)
    randomSeed = 1234;
    rng(randomSeed,'twister'); % seeding for repeatability

    allDirectories = dir(dataDirectory);
    speaker_index = 1;
    speakers = [];

    for i=1:length(allDirectories) % male or female level
       if allDirectories(i).isdir && strncmp(allDirectories(i).name,'.',1) == 0
           type = allDirectories(i).name;
           display(strcat(dataDirectory,'/' , type));
           allSubDirectories = dir(strcat(dataDirectory,'/' , type));

            for j=1:length(allSubDirectories) % speaker level
               display(allSubDirectories(j).name)
               if allSubDirectories(j).isdir && strncmp(allSubDirectories(j).name,'.',1) == 0
                   speaker_name = allSubDirectories(j).name
                   speakers(speaker_index).type = type;
                   speakers(speaker_index).name = speaker_name;
                   speakers(speaker_index).trainingSignals = {};
                   speakers(speaker_index).testSignals = {};
                   recordings_paths = strcat(dataDirectory,'/',type ,'/' ,speaker_name,'/*.wav');
                   allSubSubFiles = dir(recordings_paths);

                   % taking {trainingCount + testCount} random files, trainingCount for training and
                   % testCount for testing
                   n = length(allSubSubFiles);
                   display(n)
                   chosen = trainingCount + testCount;
                   chosen_indices = randperm(n,min(chosen,n));
                   
                   for k=1:length(chosen_indices)
                           file = allSubSubFiles(chosen_indices(k)).name;
                           %display(file);
                           file_path = strcat(dataDirectory ,'/', type, '/', speaker_name, '/',file);
                           [file_content, original_freq] = audioread(file_path);
                           size(file_content)
                           file_content = resample(file_content,resample_freq,original_freq);
                           file_content = pad_mirror(file_content,Npad);
                       if k <= trainingCount
                           speakers(speaker_index).trainingSignals{k} = file_content;
                           
                       else
                           speakers(speaker_index).testSignals{(length(speakers(speaker_index).testSignals)+1)} = file_content;
                       end
                   end
                   speakers(speaker_index).fullTraining = cell2mat(speakers(speaker_index).trainingSignals);

                   speaker_index = speaker_index+1;
               end
           end
       end
    end
end