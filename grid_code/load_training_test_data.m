function speakers = load_training_test_data(dataDirectory, trainingCount, testCount, resample_freq, Npad)
    randomSeed = 1234;
    rng(randomSeed,'twister'); % seeding for repeatabilitys

    allDirectories = dir(dataDirectory);
    speaker_index = 1;
    
    % define the data structure fields
    speakers = struct('type','','name','','trainingFileNames',{},'testingFileNames',{},'fullTraining',[],'fullTesting',[]);

    for i=1:length(allDirectories) % male or female level
        % make sure it's a directory but not the . or .. directory
        if allDirectories(i).isdir && strncmp(allDirectories(i).name,'.',1) == 0
            type = allDirectories(i).name;
            display(strcat(dataDirectory,'/' , type));
            allSubDirectories = dir(strcat(dataDirectory,'/' , type));

            for j=1:length(allSubDirectories) % speaker level
                display(allSubDirectories(j).name)
                if allSubDirectories(j).isdir && strncmp(allSubDirectories(j).name,'.',1) == 0
                    speaker_name = allSubDirectories(j).name
                    new_speaker = struct;
                    new_speaker.type = type;
                    new_speaker.name = speaker_name;
                    new_speaker.trainingSignals = {};
                    new_speaker.testingSignals = {};
                    new_speaker.trainingFileNames = {};
                    new_speaker.testingFileNames = {};
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
                        
                        % add file content to either training signals or
                        % test signal
                        if k <= trainingCount
                            new_speaker.trainingSignals{k} = file_content;
                            new_speaker.trainingFileNames{k} = file;
                        else
                            new_speaker.testingSignals{k-trainingCount} = file_content;
                            new_speaker.testingFileNames{k-trainingCount} = file;
                        end
                   end
                   new_speaker.fullTraining = cell2mat(new_speaker.trainingSignals);
                   new_speaker.fullTesting = cell2mat(new_speaker.testingSignals);
                   
                   % remove fields we won't be using
                   new_speaker = rmfield(new_speaker,'trainingSignals');
                   new_speaker = rmfield(new_speaker,'testingSignals');
                   
                   speakers(speaker_index) = new_speaker;

                   speaker_index = speaker_index+1;
               end
           end
       end
    end
end