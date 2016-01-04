
%%

%%% Training Params
trainingCount = 1;
testingCount = 1;


% General params
eps=1e-3;
fs = 16000; % resample freq
Npad = 2^18; % max size of audio clip


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
numberOfIterations = 2000;

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
dataDirectory = '../data/MASS/tamy-que_pena_tanto_faz (bossanova)/tamy-que_pena_tanto_faz_6-19';
saveDirectory = strcat('final/mass/tamy_que_pena/','T_',int2str(T),'_D1_',int2str(dictionaryAtomsX1),'_D2_',int2str(dictionaryAtomsX2),'_Npad_',int2str(Npad),'_fs_',int2str(fs),'/');
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

[sourcesL,sourcesR,sourceNameList,mixL,mixR,trainingFs,nbits] = loadSources(dataDirectory,0);

sources = [];

for i=1:size(sourceNameList,1)
    file_content = ((sourcesL(:,i) + sourcesR(:,i))/2);
    size(file_content)
    file_content = resample(file_content,fs, trainingFs);
    size(file_content)
    file_content = pad_mirror(file_content,Npad);
    sources = [sources file_content];
end
                            

%% Creating the scatt matrices

for i=1:size(sourceNameList,1)
    
    X1 = [];
    X2 = [];
    display(strcat('Started Speaker: ', sourceNameList(i).name))
    [X2, X1] = audioscatt_fwd_haar(sources(:,i), filts, options);
    
    X1 = reshape(X1, [size(X1,1), size(X1,2) * size(X1,3)]);
    X2 = reshape(X2, [size(X2,1), size(X2,2) * size(X2,3)]);
    
    display(strcat('Finished Speaker: ', sourceNameList(i).name))
    save(strcat(scattsSaveDirectory, 'scatt_n', int2str(trainingCount),'_T_',int2str(T),'_Q_',int2str(scparam.Q),'_Npad_',int2str(Npad), '__', 'speaker_', sourceNameList(i).name),'X1', 'X2')
    clear X1 X2
end


%%% Calculating the standard deviations for later normalization

[stds1, stds2] = calculate_stds(scattsSaveDirectory);

save(strcat(paramsSaveDirectory, 'renorm_params.mat'), 'stds1','stds2');


%%% Creating the NMF dictionaries

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
        
        save(strcat(dictionariesSaveDirectory,'dict_n', int2str(trainingCount), '__', fileName(strfind(fileName,'__')+2:end)),'Dnmf1','Dnmf2')
        
    end
end


%%% Testing

% load testing file
testDirectory = '../data/MASS/tamy-que_pena_tanto_faz (bossanova)/tamy-que_pena_tanto_faz_46-57';

[testSourcesL,testSourcesR,testSourceNameList,testMixL,testMixR,testingFs,nbits] = loadSources(testDirectory,0);

testSources = [];

for i=1:size(testSourceNameList,1)
    file_content = ((testSourcesL(:,i) + testSourcesR(:,i))/2);
    file_content = resample(file_content,fs, testingFs);
    file_content = pad_mirror(file_content,Npad);
    testSources = [testSources file_content];
end

testMix = ((testMixL + testMixR) / 2)';

%%%

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

[lvl_2_speech_list, lvl_1_speech_list] = demix_scatt2top_multi(testMix', Dnmf1_list, Dnmf2_list, stds1, stds2, eps, filts, scparam, param1, param2, Npad);


%%

result = struct;
result.lvl_1_speech_list = lvl_1_speech_list;
result.lvl_2_speech_list = lvl_2_speech_list;
result.mix = testMix;
result.sources = testSources;
result.sourceNames = testSourceNameList;

save(strcat(resultsSaveDirectory, 'results.mat'), 'result');
%%

soundsc(testMix,testingFs);

%%

soundsc(lvl_2_speech_list{1},testingFs);

%%

soundsc(lvl_2_speech_list{2},testingFs);