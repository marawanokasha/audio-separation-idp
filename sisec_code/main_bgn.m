
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
scparam.N = Npad;
scparam.T = T;
scparam.Q = 32;
filts = create_scattfilters( scparam );

options.renorm=1;
options.parallel = 1;


%%% NMF Params
dictionaryAtomsX1 = 200;
dictionaryAtomsX2 = 800;
lambda = 0.1;
numberOfIterations = 1000;

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
dataDirectory = '../data/sisec/bgn';
saveDirectory = strcat('final/sisec/bgn/','T_',int2str(T),'_D1_',int2str(dictionaryAtomsX1),'_D2_',int2str(dictionaryAtomsX2),'_Npad_',int2str(Npad),'_fs_',int2str(fs),'/');
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


save(strcat(paramsSaveDirectory, 'general_params.mat'), 'trainingCount', 'testingCount','eps','Npad','fs','param1','param2','scparam','options')

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

[stds1, stds2] = calculate_stds(scattsSaveDirectory);

save(strcat(paramsSaveDirectory, 'renorm_params.mat'), 'stds1','stds2');


%% Creating the NMF dictionaries

allScattFiles = dir(scattsSaveDirectory);

% Create a dictionary from the scattering representation of each speaker
for i=1:length(allScattFiles)
    % make sure it's a directory but not the . or .. directory
    if allScattFiles(i).isdir == false && strncmp(allScattFiles(i).name,'.',1) == 0
        
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

testDirectory = '../data/sisec/bgn/test/';

testFileName = 'test_Ca1_Co_A_mix.wav';
file_path = strcat(testDirectory,testFileName); 
[fc, original_freq] = audioread(file_path);
fc = resample(fc,fs,original_freq);
fc = pad_mirror(fc,Npad)';

%%


%% Testing
%%% Using predefined dictionaries -> Result is great

testDirectory = '../data/sisec/bgn/test/';

testFileName = 'test_Ca1_Co_A_mix.wav';
file_path = strcat(testDirectory,testFileName); 
test_file_content = read_audio_processed(file_path,fs,Npad);

[Dnmf11, Dnmf21] = get_dictionaries('f1', dictionariesSaveDirectory);
[Dnmf12, Dnmf22] = get_dictionaries('n1', dictionariesSaveDirectory);
            
[lvl2_speaker1, lvl2_speaker2, lvl1_speaker1, lvl1_speaker2] = demix_scatt2top(test_file_content, Dnmf11, Dnmf12, Dnmf21, Dnmf22, stds1, stds2, eps, filts, scparam, param1, param2, Npad);

%%
soundsc(test_file_content, fs);

%%
soundsc(lvl2_speaker1, fs)

%%
soundsc(lvl2_speaker2, fs)

%% Testing
%%% Using all dictionaries


testDirectory = '../data/sisec/bgn/test/';

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

testDirectory = '../data/sisec/bgn/test/';
testSourcesDirectory = '../data/sisec/bgn/test-speaker-sources/';
testNoisesDirectory = '../data/sisec/bgn/test-noise-sources/';

allTestFiles = dir(testDirectory);
results = [];
result_index = 1;

% Create a dictionary from the scattering representation of each speaker
for i=1:length(allTestFiles)
    % make sure it's a directory but not the . or .. directory
    if allTestFiles(i).isdir == false && strncmp(allTestFiles(i).name,'.',1) == 0
        file_name = allTestFiles(i).name;
        ind = strfind(file_name, '_');
        
        file_path = strcat(testDirectory, file_name);
        test_source_file_path = strcat(testSourcesDirectory, file_name(1:ind(end)), 'sim.wav');
        test_noise_file_path = strcat(testNoisesDirectory, file_name(1:ind(end)), 'noi.wav');
        
        test_file_content = read_audio_processed(file_path,fs,Npad);
        test_source_file_content = read_audio_processed(test_source_file_path,fs,Npad);
        test_noise_file_content = read_audio_processed(test_noise_file_path,fs,Npad);
        
        if strfind(file_name, 'A')
            speaker_name = 'f1';
            [Dnmf11, Dnmf21] = get_dictionaries('f1', dictionariesSaveDirectory);
        else
            speaker_name = 'm1';
            [Dnmf11, Dnmf21] = get_dictionaries('m1', dictionariesSaveDirectory);
        end
        
        [Dnmf12, Dnmf22] = get_dictionaries('n1', dictionariesSaveDirectory);

        [lvl2_speaker1, lvl2_speaker2, lvl1_speaker1, lvl1_speaker2] = demix_scatt2top(test_file_content, Dnmf11, Dnmf12, Dnmf21, Dnmf22, stds1, stds2, eps, filts, scparam, param1, param2, Npad);
        
        % measuring signal DOA
        [SDR_lvl1,SIR_lvl1,SAR_lvl1] = bss_eval_source_denoising(lvl1_speaker1, test_source_file_content', test_noise_file_content');
        [SDR_lvl2,SIR_lvl2,SAR_lvl2] = bss_eval_source_denoising(lvl2_speaker1, test_source_file_content', test_noise_file_content');
        
        % measuring signal separation
        [SDRi_lvl1,ISRi_lvl1,SIRi_lvl1,SARi_lvl1] = bss_eval_images_nosort([lvl1_speaker1;lvl1_speaker2],[test_source_file_content';test_noise_file_content']);
        [SDRi_lvl2,ISRi_lvl2,SIRi_lvl2,SARi_lvl2] = bss_eval_images_nosort([lvl2_speaker1;lvl2_speaker2],[test_source_file_content';test_noise_file_content']);
        
        bss_DOA_results_lvl1 = struct; bss_DOA_results_lvl2 = struct; 
        bss_DOA_results_lvl1.SDR = SDR_lvl1;
        bss_DOA_results_lvl1.SIR = SIR_lvl1;
        bss_DOA_results_lvl1.SAR = SAR_lvl1;
        bss_DOA_results_lvl2.SDR = SDR_lvl2;
        bss_DOA_results_lvl2.SIR = SIR_lvl2;
        bss_DOA_results_lvl2.SAR = SAR_lvl2;
        
        bss_sig_results_lvl1.SDRi = SDRi_lvl1;
        bss_sig_results_lvl1.ISRi = ISRi_lvl1;
        bss_sig_results_lvl1.SIRi = SIRi_lvl1;
        bss_sig_results_lvl1.SARi = SARi_lvl1;
        bss_sig_results_lvl2.SDRi = SDRi_lvl2;
        bss_sig_results_lvl2.ISRi = ISRi_lvl2;
        bss_sig_results_lvl2.SIRi = SIRi_lvl2;
        bss_sig_results_lvl2.SARi = SARi_lvl2;
        
        bss_DOA_results_lvl1
        bss_DOA_results_lvl2
        bss_sig_results_lvl1
        bss_sig_results_lvl2
        
        results(result_index).speaker_name = speaker_name;
        results(result_index).file_name = file_name;
        results(result_index).speaker_original = test_source_file_content;
        results(result_index).noise_original = test_noise_file_content;
        
        results(result_index).speaker_reconstructed_lvl1 = lvl1_speaker1;
        results(result_index).noise_reconstructed_lvl1 = lvl1_speaker2;
        
        results(result_index).speaker_reconstructed_lvl2 = lvl2_speaker1;
        results(result_index).noise_reconstructed_lvl2 = lvl2_speaker2;
        
        results(result_index).bss_DAO_results_lvl1 = bss_DOA_results_lvl1;
        results(result_index).bss_DAO_results_lvl2 = bss_DOA_results_lvl2;
        results(result_index).bss_sig_results_lvl1 = bss_sig_results_lvl1;
        results(result_index).bss_sig_results_lvl2 = bss_sig_results_lvl2;
        
        result_index = result_index+1;
        display(strcat('Finished speaker: ', int2str(result_index)));
    end
end

save(strcat(resultsSaveDirectory,'results.mat'), 'results','-v7.3')

%%
soundsc(test_file_content, fs);

%%
soundsc(lvl2_speaker1, fs)

%%
soundsc(lvl2_speaker2, fs)
%%

soundsc(results(test_index).speaker2_original, fs);

%%

soundsc(results(test_index).speaker1_reconstructed_lvl2, fs);

%%

soundsc(results(test_index).speaker2_reconstructed_lvl2, fs);
        
%%

results(test_index).bss_results_lvl1


results(test_index).bss_results_lvl2