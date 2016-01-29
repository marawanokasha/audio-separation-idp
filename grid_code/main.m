%% For GRID Data
%%

%%% Training Params
trainingCount = 100;
testingCount = 5;


% General params
eps=1e-3;
fs = 16000; % resample freq
Npad = 2^15; % max size of audio clip


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
dataDirectory = '../data/grid';
saveDirectory = 'final/grid/';
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

speakers = load_training_test_data(dataDirectory, trainingCount, testingCount, fs, Npad);

%%
save(strcat(speakersSaveDirectory, 'speakers_resampled_','training_',int2str(trainingCount),'_testing_',int2str(testingCount),'.mat'), 'speakers')

%%
%load(strcat(saveDirectory, 'speakers/speakers_resampled_','training_',trainingCount,'_testing_',testCount,'_.mat'))

%% Creating the scatt matrices

for i=1:size(speakers,2)
    
    X1 = [];
    X2 = [];
    display(strcat('Started Speaker ', speakers(i).name))
    [X2, X1] = audioscatt_fwd_haar(speakers(i).fullTraining, filts, options);
    
    X1 = reshape(X1, [size(X1,1), size(X1,2) * size(X1,3)]);
    X2 = reshape(X2, [size(X2,1), size(X2,2) * size(X2,3)]);
    
    display(strcat('Finished Speaker ', speakers(i).name))
    save(strcat(scattsSaveDirectory, 'scatt_n', int2str(trainingCount),'_T_',int2str(T),'_Q_',int2str(scparam.Q),'_Npad_',int2str(Npad), '__',speakers(i).type, '_speaker_', speakers(i).name),'X1', 'X2')
    clear X1 X2
end


%% Calculating the standard deviations for later normalization

[stds1, stds2] = calculate_stds(scattsSaveDirectory);

%save(strcat(paramsSaveDirectory, 'renorm_params.mat'), 'stds1','stds2');


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


%% Testing

results = [];
resultsIndex = 1;
for i=1:size(speakers,2)
    
    % Load dictionaries of the first speaker
    [Dnmf11, Dnmf21] = get_dictionaries(speakers(i).name, dictionariesSaveDirectory);
    
    for j=1:size(speakers(i).fullTesting,2)
        
        x1 = speakers(i).fullTesting(:,j)'; T1 = length(x1 );
        
        for k=i+1:size(speakers,2)

            % Load dictionaries of the second speaker
            [Dnmf12, Dnmf22] = get_dictionaries(speakers(k).name, dictionariesSaveDirectory);
            
            for l=1:size(speakers(k).fullTesting,2)
                
                x2 = speakers(k).fullTesting(:,l)'; T2 = length(x2);

                % create the mixture
                T = min([T1,T2,Npad]);
                x1 = x1(1:T);
                x2 = x2(1:T);
                mix = (x1+x2);
                
                [speech1, speech2, xest1, xest2] = demix_scatt2top(mix, Dnmf11, Dnmf12, Dnmf21, Dnmf22, stds1, stds2, eps, filts, scparam, param1, param2, Npad);

                % evaluate results
                bss_scatt2  =  BSS_EVAL(x1', x2', speech1(1:T)', speech2(1:T)', mix');
                bss_scatt1 =  BSS_EVAL(x1', x2', xest1(1:T)', xest2(1:T)', mix');
                
                results(resultsIndex).speaker1 = speakers(i).name;
                results(resultsIndex).speaker2 = speakers(k).name;
                
                results(resultsIndex).speaker1_file = speakers(i).testingFileNames{j} ;
                results(resultsIndex).speaker2_file = speakers(k).testingFileNames{l} ;
                
                results(resultsIndex).speaker1_original = x1 ;
                results(resultsIndex).speaker2_original = x2 ;
                
                results(resultsIndex).speaker1_reconstructed_lvl2 = speech1 ;
                results(resultsIndex).speaker2_reconstructed_lvl2 = speech2 ;
                
                results(resultsIndex).speaker1_reconstructed_lvl1 = xest1;
                results(resultsIndex).speaker2_reconstructed_lvl1 = xest2;
                
                results(resultsIndex).bss_results_lvl1 = bss_scatt1;
                results(resultsIndex).bss_results_lvl2 = bss_scatt2;
                
                display(strcat('finished result: ', int2str(resultsIndex)));
                resultsIndex = resultsIndex + 1;
                break;

            end
        end

    end
end

%%
save(strcat(resultsSaveDirectory,'results.mat'), 'results','-v7.3')


%%
test_index = 3;

%%
soundsc(results(test_index).speaker1_original, fs);

%%

soundsc(results(test_index).speaker2_original, fs);

%%

soundsc(results(test_index).speaker1_reconstructed_lvl2, fs);

%%

soundsc(results(test_index).speaker2_reconstructed_lvl2, fs);
        
%%

results(test_index).bss_results_lvl1


results(test_index).bss_results_lvl2

%% Writing out sample results

resultsDirectory = '../results_audio/grid/';

wavwrite(results(test_index).speaker1_original, fs, strcat(resultsDirectory,'speaker1_original.wav'));
wavwrite(results(test_index).speaker2_original, fs, strcat(resultsDirectory,'speaker2_original.wav'));
wavwrite(results(test_index).speaker1_reconstructed_lvl2, fs, strcat(resultsDirectory,'speaker1_reconstructed.wav'));
wavwrite(results(test_index).speaker2_reconstructed_lvl2, fs, strcat(resultsDirectory,'speaker2_reconstructed.wav'));

mix = results(test_index).speaker1_original + results(test_index).speaker2_original;
wavwrite(mix, fs, strcat(resultsDirectory,'mix.wav'));



%% Testing
%%% For 3 speakers

results = [];
resultsIndex = 1;
for i=1:size(speakers,2)
    
    % Load dictionaries of the first speaker
    [Dnmf11, Dnmf21] = get_dictionaries(speakers(i).name, dictionariesSaveDirectory);
    
    for j=1:size(speakers(i).fullTesting,2)
        
        x1 = speakers(i).fullTesting(:,j)'; T1 = length(x1 );
        
        for k=i+1:size(speakers,2)

            % Load dictionaries of the second speaker
            [Dnmf12, Dnmf22] = get_dictionaries(speakers(k).name, dictionariesSaveDirectory);
            
            for l=1:size(speakers(k).fullTesting,2)
                
                x2 = speakers(k).fullTesting(:,l)'; T2 = length(x2);
                
                for m=k+1:size(speakers,2)

                    % Load dictionaries of the third speaker
                    [Dnmf13, Dnmf23] = get_dictionaries(speakers(m).name, dictionariesSaveDirectory);

                    for n=1:size(speakers(m).fullTesting,2)

                        x3 = speakers(m).fullTesting(:,n)'; T3 = length(x3);

                        % create the mixture
                        T = min([T1,T2,T3,Npad]);
                        x1 = x1(1:T);
                        x2 = x2(1:T);
                        x3 = x3(1:T);
                        mix = (x1+x2+x3);
                        
                        Dnmf1_list = {};
                        Dnmf2_list = {};

                        Dnmf1_list = [Dnmf1_list Dnmf11 Dnmf12 Dnmf13];
                        Dnmf2_list = [Dnmf2_list Dnmf21 Dnmf22 Dnmf23];           

                        [lvl2_speech_list, lvl1_speech_list] = demix_scatt2top_multi(mix, Dnmf1_list, Dnmf2_list, stds1, stds2, eps, filts, scparam, param1, param2, Npad);

                        results(resultsIndex).speaker1 = speakers(i).name;
                        results(resultsIndex).speaker2 = speakers(k).name;
                        results(resultsIndex).speaker3 = speakers(m).name;

                        results(resultsIndex).speaker1_file = speakers(i).testingFileNames{j} ;
                        results(resultsIndex).speaker2_file = speakers(k).testingFileNames{l} ;
                        results(resultsIndex).speaker3_file = speakers(m).testingFileNames{n} ;
                        
                        results(resultsIndex).mix = mix;
                        results(resultsIndex).speaker1_original = x1 ;
                        results(resultsIndex).speaker2_original = x2 ;
                        results(resultsIndex).speaker3_original = x3 ;
            
                        results(resultsIndex).speaker1_reconstructed_lvl1 = lvl1_speech_list{1};
                        results(resultsIndex).speaker2_reconstructed_lvl1 = lvl1_speech_list{2};
                        results(resultsIndex).speaker3_reconstructed_lvl1 = lvl1_speech_list{3};

                        results(resultsIndex).speaker1_reconstructed_lvl2 = lvl2_speech_list{1};
                        results(resultsIndex).speaker2_reconstructed_lvl2 = lvl1_speech_list{2};
                        results(resultsIndex).speaker3_reconstructed_lvl2 = lvl1_speech_list{3};

                        display(strcat('finished result: ', int2str(resultsIndex)));
                        resultsIndex = resultsIndex + 1;
                        break;
                    end
                    break;
                end
                break;
            end
            break;
        end
        break;

    end
    break;
end

%save(strcat(resultsSaveDirectory,'results_3speakers.mat'), 'results','-v7.3')

%%
soundsc(mix,fs);

%%
soundsc(results(1).speaker3_reconstructed_lvl2,fs);
%%
soundsc(results(7).speaker2_reconstructed_lvl2, fs)


%% Testing
%%% For 3 speakers (selected speakers)

results = [];
resultsIndex = 1;


% Load dictionaries of the first speaker
[Dnmf11, Dnmf21] = get_dictionaries('s1', dictionariesSaveDirectory);
speaker1 = get_speaker(speakers, 's1');

% Load dictionaries of the second speaker
[Dnmf12, Dnmf22] = get_dictionaries('s2', dictionariesSaveDirectory);
speaker2 = get_speaker(speakers, 's2');

% Load dictionaries of the third speaker
[Dnmf13, Dnmf23] = get_dictionaries('s7', dictionariesSaveDirectory);
speaker3 = get_speaker(speakers, 's7');

Dnmf1_list = {};
Dnmf2_list = {};

Dnmf1_list = [Dnmf1_list Dnmf11 Dnmf12 Dnmf13];
Dnmf2_list = [Dnmf2_list Dnmf21 Dnmf22 Dnmf23];           

x1 = speaker1.fullTesting(:,1)'; T1 = length(x1 );
x2 = speaker2.fullTesting(:,1)'; T2 = length(x2);
x3 = speaker3.fullTesting(:,1)'; T3 = length(x3);

% create the mixture
T = min([T1,T2,T3,Npad]);
x1 = x1(1:T);
x2 = x2(1:T);
x3 = x3(1:T);
mix = (x1+x2+x3);

soundsc(mix,fs);

%%
soundsc(x1,fs);
%%
soundsc(x2,fs);
%%
soundsc(x3,fs);
%%
[lvl2_speech_list, lvl1_speech_list] = demix_scatt2top_multi(mix, Dnmf1_list, Dnmf2_list, stds1, stds2, eps, filts, scparam, param1, param2, Npad);

%save(strcat(resultsSaveDirectory,'results_3speakers.mat'), 'results','-v7.3')

%%
soundsc(lvl2_speech_list{3}, fs)