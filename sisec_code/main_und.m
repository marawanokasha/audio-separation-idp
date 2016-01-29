
%%

%%% Training Params
trainingCount = 100;
testingCount = 0;


% General params
eps=1e-3;
fs = 16000; % resample freq
Npad = 2^17; % max size of audio clip


%%% Scattering Params
T = 2048;
Q = 32;
scparam.N = Npad;
scparam.T = T;
scparam.Q = Q;
filts = create_scattfilters( scparam );

options.renorm=1;
options.parallel = 1;


%%% NMF Params
dictionaryAtomsX1 = 200;
dictionaryAtomsX2 = 800;
lambda = 0.1;
numberOfIterations = 500;

% Params for X1: 1st level scattering
param1.K = [dictionaryAtomsX1];
param1.posAlpha = 1; % dictionaries must be positive
param1.posD = 1; % activations must be positive
param1.pos = 1;
param1.lambda = [lambda];
param1.iter = numberOfIterations;
param1.numThreads=16;
param1.batchsize=512;

% Params for X2: 2nd level scattering
param2.K = dictionaryAtomsX2;
param2.posAlpha = 1;
param2.posD = 1;
param2.pos = 1;
param2.lambda = lambda;
param2.iter = numberOfIterations;
param2.numThreads=16;
param2.batchsize=512;


% setting up output directories
dataDirectory = '../data/sisec/und/training';
saveDirectory = strcat('final/sisec/und/','T_',int2str(T),'_D1_',int2str(dictionaryAtomsX1),'_D2_',int2str(dictionaryAtomsX2),'_nIter_',int2str(numberOfIterations),'_Npad_',int2str(Npad),'_fs_',int2str(fs),'/');
speakersSaveDirectory = strcat(saveDirectory, 'speakers/');
scattsSaveDirectory = strcat(saveDirectory, 'saved_scatts/');
paramsSaveDirectory = strcat(saveDirectory, 'params/');
dictionariesSaveDirectory = strcat(saveDirectory, 'saved_dicts/');
resultsSaveDirectory = strcat(saveDirectory, 'results/');

mkdir(saveDirectory);
mkdir(speakersSaveDirectory);
mkdir(scattsSaveDirectory);
mkdir(paramsSaveDirectory);
mkdir(dictionariesSaveDirectory);
mkdir(resultsSaveDirectory);


%save(strcat(paramsSaveDirectory, 'general_params.mat'), 'trainingCount', 'testingCount','eps','Npad','fs','param1','param2','scparam','options')

%% Generate the speakers matrix

speakers = load_speakers_mono(dataDirectory, trainingCount, testingCount, fs, Npad);

%%
save(strcat(speakersSaveDirectory, 'speakers_resampled_','training_',int2str(trainingCount),'_testing_',int2str(testingCount),'.mat'), 'speakers')

%%
%load(strcat(saveDirectory, 'speakers/speakers_resampled_','training_',trainingCount,'_testing_',testCount,'_.mat'))

%% Creating the scatt matrices

for i=1:size(speakers,2)
    
    X1 = [];
    X2 = [];
    display(strcat('Started Speaker: ', speakers(i).name))
    [X2, X1] = audioscatt_fwd_haar(speakers(i).fullTraining, filts, options);
    
    X1 = reshape(X1, [size(X1,1), size(X1,2) * size(X1,3)]);
    X2 = reshape(X2, [size(X2,1), size(X2,2) * size(X2,3)]);
    
    display(strcat('Finished Speaker: ', speakers(i).name))
    save(strcat(scattsSaveDirectory, 'scatt_n', int2str(trainingCount),'_T_',int2str(T),'_Q_',int2str(scparam.Q),'_Npad_',int2str(Npad), '__',speakers(i).type, '_speaker_', speakers(i).name),'X1', 'X2')
    clear X1 X2
end


%% Calculating the standard deviations for later normalization

%speakers_to_consider = {'f1','f2','f3'};
%[stds1, stds2] = calculate_stds(scattsSaveDirectory, speakers_to_consider);
%save(strcat(paramsSaveDirectory, 'renorm_params.mat'), 'stds1','stds2');

[stds1, stds2] = calculate_stds(scattsSaveDirectory);

%save(strcat(paramsSaveDirectory, 'renorm_params.mat'), 'stds1','stds2');


%% Creating the NMF dictionaries

allScattFiles = dir(scattsSaveDirectory);

% Create a dictionary from the scattering representation of each speaker
for i=1:length(allScattFiles)
    % make sure it's a directory but not the . or .. directory
    if allScattFiles(i).isdir == false && strncmp(allScattFiles(i).name,'.',1) == 0
        
        % make sure it's one of the speakers we want to use
        if length(strmatch(get_speaker_name_from_file(allScattFiles(i).name), speakers_to_consider)) == 0
            continue;
        end
        
        fileName = allScattFiles(i).name;
        load(strcat(scattsSaveDirectory, fileName));
        data.X1 = X1;
        data.X2 = X2;
        if options.renorm
           %renormalize data: whiten each frequency component
           data.X1 = renorm_spect_data(data.X1, stds1, eps);
           data.X2 = renorm_spect_data(data.X2, stds2, eps);
        end
        
        % create the dictionaries
        Dnmf1 = mexTrainDL(abs(data.X1),param1);
        Dnmf2 = mexTrainDL(abs(data.X2),param2);
        
        save(strcat(dictionariesSaveDirectory,'dict_n', int2str(trainingCount), '_', fileName(strfind(fileName,'__')+2:end)),'Dnmf1','Dnmf2')
        
    end
end

%%



%% Testing
%%% Using predefined dictionaries for 3 speakers

testDirectory = '../data/sisec/und/test/test/';

testFileName = 'test_female3_liverec_130ms_1m_mix.wav';
file_path = strcat(testDirectory,testFileName); 
test_file_content = read_audio_processed(file_path,fs,Npad);

%%
Dnmf1_list = {};
Dnmf2_list = {};

[Dnmf11, Dnmf21] = get_dictionaries('f1', dictionariesSaveDirectory);
[Dnmf12, Dnmf22] = get_dictionaries('f2', dictionariesSaveDirectory);
[Dnmf13, Dnmf23] = get_dictionaries('f3', dictionariesSaveDirectory);
 
Dnmf1_list = [Dnmf1_list Dnmf11 Dnmf12 Dnmf13];
Dnmf2_list = [Dnmf2_list Dnmf21 Dnmf22 Dnmf23];           

[lvl2_speech_list, lvl1_speech_list] = demix_scatt2top_multi(test_file_content, Dnmf1_list, Dnmf2_list, stds1, stds2, eps, filts, scparam, param1, param2, Npad);

%%
soundsc(test_file_content, fs);

%%
soundsc(lvl2_speech_list{1}, fs)

%%
soundsc(lvl2_speech_list{2}, fs)

%%
soundsc(lvl2_speech_list{3}, fs)

%%

testSourcesDirectory = '../data/sisec/und/test/test-speaker-sources/';
source_mat = []; % used in bss
allSourceFiles = dir(testSourcesDirectory);

ind = strfind(testFileName, 'mix.wav');
prefix = testFileName(1:ind-1);

% Getting the speakers to be used for this mix
for j=1:length(allSourceFiles)
    source_file_name = allSourceFiles(j).name;
    if strfind(source_file_name, prefix)
        source_file_path = strcat(testSourcesDirectory, source_file_name);
        source_file_content = read_audio_processed(source_file_path,fs,Npad);
        source_mat = [source_mat; source_file_content'];
    end
end


lvl1_speech_mat = [];
lvl2_speech_mat = [];
for k=1:length(lvl2_speech_list)
    lvl1_speech_mat = [lvl1_speech_mat; lvl1_speech_list{k}];
    lvl2_speech_mat = [lvl2_speech_mat; lvl2_speech_list{k}];
end
% measuring signal separation
[SDRi_lvl1,ISRi_lvl1,SIRi_lvl1,SARi_lvl1] = bss_eval_images_nosort(lvl1_speech_mat,source_mat);
[SDRi_lvl2,ISRi_lvl2,SIRi_lvl2,SARi_lvl2] = bss_eval_images_nosort(lvl2_speech_mat,source_mat);



%% Testing
%%% Using predefined dictionaries for 3 speakers and on doing manual
%%% mix

sourcePath1 = '../data/sisec/und/training/female/f1/dev1_female3_inst_sim_1.wav';
sourcePath2 = '../data/sisec/und/training/female/f2/dev1_female3_inst_sim_2.wav';
sourcePath3 = '../data/sisec/und/training/female/f3/dev1_female3_inst_sim_3.wav';

source1 = read_audio_processed(sourcePath1,fs,Npad);
source2 = read_audio_processed(sourcePath2,fs,Npad);
source3 = read_audio_processed(sourcePath3,fs,Npad);


% create the mixture
mix = (source1+source2+source3);
%%
soundsc(mix,fs);

%%
soundsc(source1, fs);
%%
soundsc(source2, fs);
%%
soundsc(source3, fs);
%%
Dnmf1_list = {};
Dnmf2_list = {};

[Dnmf11, Dnmf21] = get_dictionaries('f1', dictionariesSaveDirectory);
[Dnmf12, Dnmf22] = get_dictionaries('f2', dictionariesSaveDirectory);
[Dnmf13, Dnmf23] = get_dictionaries('f3', dictionariesSaveDirectory);
 
Dnmf1_list = [Dnmf1_list Dnmf11 Dnmf12 Dnmf13];
Dnmf2_list = [Dnmf2_list Dnmf21 Dnmf22 Dnmf23];           

[lvl2_speech_list, lvl1_speech_list] = demix_scatt2top_multi(mix, Dnmf1_list, Dnmf2_list, stds1, stds2, eps, filts, scparam, param1, param2, Npad);

display('Finished demix')

source_mat = [source1';source2';source3'];

lvl1_speech_mat = [];
lvl2_speech_mat = [];
for k=1:length(lvl2_speech_list)
    lvl1_speech_mat = [lvl1_speech_mat; lvl1_speech_list{k}];
    lvl2_speech_mat = [lvl2_speech_mat; lvl2_speech_list{k}];
end

% measuring signal separation
[SDRi_lvl1,ISRi_lvl1,SIRi_lvl1,SARi_lvl1] = bss_eval_images_nosort(lvl1_speech_mat,source_mat);
[SDRi_lvl2,ISRi_lvl2,SIRi_lvl2,SARi_lvl2] = bss_eval_images_nosort(lvl2_speech_mat,source_mat);

%%
soundsc(mix, fs);

%%
soundsc(lvl2_speech_list{1}, fs)

%%
soundsc(lvl2_speech_list{2}, fs)

%%
soundsc(lvl2_speech_list{3}, fs)

%%

testSourcesDirectory = '../data/sisec/und/test/test-speaker-sources/';
source_mat = []; % used in bss
allSourceFiles = dir(testSourcesDirectory);

ind = strfind(testFileName, 'mix.wav');
prefix = testFileName(1:ind-1);

% Getting the speakers to be used for this mix
for j=1:length(allSourceFiles)
    source_file_name = allSourceFiles(j).name;
    if strfind(source_file_name, prefix)
        source_file_path = strcat(testSourcesDirectory, source_file_name);
        source_file_content = read_audio_processed(source_file_path,fs,Npad);
        source_mat = [source_mat; source_file_content'];
    end
end


lvl1_speech_mat = [];
lvl2_speech_mat = [];
for k=1:length(lvl2_speech_list)
    lvl1_speech_mat = [lvl1_speech_mat; lvl1_speech_list{k}];
    lvl2_speech_mat = [lvl2_speech_mat; lvl2_speech_list{k}];
end
% measuring signal separation
[SDRi_lvl1,ISRi_lvl1,SIRi_lvl1,SARi_lvl1] = bss_eval_images_nosort(lvl1_speech_mat,source_mat);
[SDRi_lvl2,ISRi_lvl2,SIRi_lvl2,SARi_lvl2] = bss_eval_images_nosort(lvl2_speech_mat,source_mat);



        
%% Testing
%%% Using all dictionaries


testDirectory = '../data/sisec/und/test/';

testFileName = 'test_Ca1_Co_A_mix.wav';
file_path = strcat(testDirectory,testFileName); 
[test_file_content, original_freq] = audioread(file_path);
test_file_content = resample(test_file_content,fs,original_freq);
test_file_content = pad_mirror(test_file_content,Npad)';

%%

% load dictionaries
Dnmf1_list = {};
Dnmf2_list = {};
speaker_name_list = {};
j = 1;

allDictionaryFiles = dir(dictionariesSaveDirectory);

% Create a dictionary from the scattering representation of each speaker
for i=1:length(allDictionaryFiles)
    % make sure it's a directory but not the . or .. directory
    if allDictionaryFiles(i).isdir == false && strncmp(allDictionaryFiles(i).name,'.',1) == 0

        fileName = allDictionaryFiles(i).name;

        load(strcat(dictionariesSaveDirectory, fileName));
        
        Dnmf1_list = [Dnmf1_list Dnmf1];
        Dnmf2_list = [Dnmf2_list Dnmf2];
        speaker_name_list{j} = fileName(strfind(fileName,'__')+2:strfind(fileName,'.mat')-1);
        j = j +1;

    end
end

%%
[lvl2_speech_list, lvl1_speech_list] = demix_scatt2top_multi(test_file_content, Dnmf1_list, Dnmf2_list, stds1, stds2, eps, filts, scparam, param1, param2, Npad);


%%
soundsc(test_file_content, fs);

%%
soundsc(lvl2_speech_list{1}, fs)

%%
soundsc(lvl2_speech_list{2}, fs)

%%
soundsc(lvl2_speech_list{3}, fs)

%% Testing
%%% testing with all test files and using BSS

% mode = 'dev';
mode = 'test';
% testDirectory = '../data/sisec/und/test/dev-mix/';
% testSourcesDirectory = '../data/sisec/und/test/dev-mix-speaker-sources/';
testDirectory = '../data/sisec/und/test/test/';
testSourcesDirectory = '../data/sisec/und/test/test-speaker-sources/';

allTestFiles = dir(testDirectory);
results = [];
result_index = 1;

% Create a dictionary from the scattering representation of each speaker
for i=1:length(allTestFiles)
    % make sure it's a directory but not the . or .. directory
    if allTestFiles(i).isdir == false && strncmp(allTestFiles(i).name,'.',1) == 0
        file_name = allTestFiles(i).name
        file_path = strcat(testDirectory, file_name);
        test_file_content = read_audio_processed(file_path,fs,Npad);
        
        ind = strfind(file_name, 'mix.wav');
        prefix = file_name(1:ind-1);
        
        allSourceFiles = dir(testSourcesDirectory);

        Dnmf1_list = {};
        Dnmf2_list = {};
        source_mat = []; % used in bss

        % Getting the speakers to be used for this mix
        for j=1:length(allSourceFiles)
            source_file_name = allSourceFiles(j).name;
            if strfind(source_file_name, prefix)
                if strfind(source_file_name, 'female') gender = 'f'; else gender= 'm'; end;
                number = source_file_name(strfind(source_file_name, '.wav')-1);
                speaker_name = strcat(gender, number)
                source_file_path = strcat(testSourcesDirectory, source_file_name);
                source_file_content = read_audio_processed(source_file_path,fs,Npad);
                source_mat = [source_mat; source_file_content'];
                [Dnmf1, Dnmf2] = get_dictionaries(speaker_name, dictionariesSaveDirectory);

                Dnmf1_list = [Dnmf1_list Dnmf1];
                Dnmf2_list = [Dnmf2_list Dnmf2];
            end
        end
        
        [lvl2_speech_list, lvl1_speech_list] = demix_scatt2top_multi(test_file_content, Dnmf1_list, Dnmf2_list, stds1, stds2, eps, filts, scparam, param1, param2, Npad);
        
        lvl1_speech_mat = [];
        lvl2_speech_mat = [];
        for k=1:length(lvl2_speech_list)
            lvl1_speech_mat = [lvl1_speech_mat; lvl1_speech_list{k}];
            lvl2_speech_mat = [lvl2_speech_mat; lvl2_speech_list{k}];
        end
        % measuring signal separation
        [SDRi_lvl1,ISRi_lvl1,SIRi_lvl1,SARi_lvl1] = bss_eval_images_nosort(lvl1_speech_mat,source_mat);
        [SDRi_lvl2,ISRi_lvl2,SIRi_lvl2,SARi_lvl2] = bss_eval_images_nosort(lvl2_speech_mat,source_mat);
        
        bss_sig_results_lvl1 = struct; 
        bss_sig_results_lvl2 = struct;
        bss_sig_results_lvl1.SDRi = SDRi_lvl1;
        bss_sig_results_lvl1.ISRi = ISRi_lvl1;
        bss_sig_results_lvl1.SIRi = SIRi_lvl1;
        bss_sig_results_lvl1.SARi = SARi_lvl1;
        bss_sig_results_lvl2.SDRi = SDRi_lvl2;
        bss_sig_results_lvl2.ISRi = ISRi_lvl2;
        bss_sig_results_lvl2.SIRi = SIRi_lvl2;
        bss_sig_results_lvl2.SARi = SARi_lvl2;
        
        bss_sig_results_lvl1
        bss_sig_results_lvl2
        
        results(result_index).file_name = file_name;
        results(result_index).speakers_original = source_mat;
        
        results(result_index).speakers_reconstructed_lvl1 = lvl1_speech_list;
        results(result_index).speakers_reconstructed_lvl2 = lvl2_speech_list;
        
        results(result_index).bss_sig_results_lvl1 = bss_sig_results_lvl1;
        results(result_index).bss_sig_results_lvl2 = bss_sig_results_lvl2;
        
        display(strcat('Finished test file: ', int2str(result_index)));
        result_index = result_index+1;
    end
end

save(strcat(resultsSaveDirectory,'results_', mode,'.mat'), 'results','-v7.3')

%%
test_index = 1;

%%
soundsc(results(test_index).speakers_reconstructed_lvl2{1}, fs);
%%
soundsc(results(test_index).speakers_reconstructed_lvl2{2}, fs);
%% Writing out sample results

resultsDirectory = '../results_audio/sisec/und/';

wavwrite(results(test_index).speakers_original(1,:), fs, strcat(resultsDirectory,'speaker1_original.wav'));
wavwrite(results(test_index).speakers_original(2,:), fs, strcat(resultsDirectory,'speaker2_original.wav'));
wavwrite(results(test_index).speakers_original(3,:), fs, strcat(resultsDirectory,'speaker3_original.wav'));
wavwrite(results(test_index).speakers_reconstructed_lvl2{1}, fs, strcat(resultsDirectory,'speaker1_reconstructed.wav'));
wavwrite(results(test_index).speakers_reconstructed_lvl2{2}, fs, strcat(resultsDirectory,'speaker2_reconstructed.wav'));
wavwrite(results(test_index).speakers_reconstructed_lvl2{3}, fs, strcat(resultsDirectory,'speaker3_reconstructed.wav'));

%% Writing out sample results

resultsDirectory = '../results_audio/sisec/und/test/';

wavwrite(results(test_index).speakers_original(1,:), fs, strcat(resultsDirectory,'speaker1_original.wav'));
wavwrite(results(test_index).speakers_original(2,:), fs, strcat(resultsDirectory,'speaker2_original.wav'));
wavwrite(results(test_index).speakers_original(3,:), fs, strcat(resultsDirectory,'speaker3_original.wav'));
wavwrite(results(test_index).speakers_reconstructed_lvl2{1}, fs, strcat(resultsDirectory,'speaker1_reconstructed.wav'));
wavwrite(results(test_index).speakers_reconstructed_lvl2{2}, fs, strcat(resultsDirectory,'speaker2_reconstructed.wav'));
wavwrite(results(test_index).speakers_reconstructed_lvl2{3}, fs, strcat(resultsDirectory,'speaker3_reconstructed.wav'));


%%
% mix_name='male3_liverec_130ms_1m'; % covers both male and female
% mix_name='male3_liverec_130ms_5cm'; % covers both male and female
% mix_name='male3_liverec_250ms_1m'; % covers both male and female
mix_name='male3_liverec_250ms_5cm'; % covers both male and female

avg_bss = struct;
avg_bss.lvl1_SDRi = 0;
avg_bss.lvl1_ISRi = 0;
avg_bss.lvl1_SIRi = 0;
avg_bss.lvl1_SARi = 0;
avg_bss.lvl2_SDRi = 0;
avg_bss.lvl2_ISRi = 0;
avg_bss.lvl2_SIRi = 0;
avg_bss.lvl2_SARi = 0;
count = 0;
for i=1:length(results)
    if strfind(results(i).file_name, mix_name)
        bss_sig_results_lvl1.SDRi
        bss_sig_results_lvl1 = results(i).bss_sig_results_lvl1;
        bss_sig_results_lvl2 = results(i).bss_sig_results_lvl2;
        
        avg_bss.lvl1_SDRi = avg_bss.lvl1_SDRi + bss_sig_results_lvl1.SDRi;
        avg_bss.lvl1_ISRi = avg_bss.lvl1_ISRi + bss_sig_results_lvl1.ISRi;
        avg_bss.lvl1_SIRi = avg_bss.lvl1_SIRi + bss_sig_results_lvl1.SIRi;
        avg_bss.lvl1_SARi = avg_bss.lvl1_SARi + bss_sig_results_lvl1.SARi;
        
        avg_bss.lvl2_SDRi = avg_bss.lvl2_SDRi + bss_sig_results_lvl2.SDRi;
        avg_bss.lvl2_ISRi = avg_bss.lvl2_ISRi + bss_sig_results_lvl2.ISRi;
        avg_bss.lvl2_SIRi = avg_bss.lvl2_SIRi + bss_sig_results_lvl2.SIRi;
        avg_bss.lvl2_SARi = avg_bss.lvl2_SARi + bss_sig_results_lvl2.SARi;
        count = count +1;
    end
end

avg_bss.lvl1_SDRi = mean(avg_bss.lvl1_SDRi / count);
avg_bss.lvl1_ISRi = mean(avg_bss.lvl1_ISRi / count);
avg_bss.lvl1_SIRi = mean(avg_bss.lvl1_SIRi / count);
avg_bss.lvl1_SARi = mean(avg_bss.lvl1_SARi / count);

avg_bss.lvl2_SDRi = mean(avg_bss.lvl2_SDRi / count);
avg_bss.lvl2_ISRi = mean(avg_bss.lvl2_ISRi / count);
avg_bss.lvl2_SIRi = mean(avg_bss.lvl2_SIRi / count);
avg_bss.lvl2_SARi = mean(avg_bss.lvl2_SARi / count);

save(strcat(resultsSaveDirectory, strcat('avg_results_', mode,'_', mix_name, '.mat')), 'avg_bss','-v7.3');