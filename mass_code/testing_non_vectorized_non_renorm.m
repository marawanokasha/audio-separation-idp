

%% General parameters

% define scattering parameters and construct filters
eps=1e-3;
fs = 16000;
Npad = 2^17;
T = 2048;
scparam.N = Npad;
scparam.T = T;
scparam.Q = 32;
filts = create_scattfilters( scparam );


%% NMF Params

KK1 = [160];
LL1 = [0.1];
param1.K = KK1;
param1.posAlpha = 1;
param1.posD = 1;
param1.pos = 1;
param1.lambda = LL1;
param1.iter = 4000;
param1.numThreads=16;
param1.batchsize=512;

KK2 = [768];
LL2 = [0.1];
param2.K = KK2;
param2.posAlpha = 1;
param2.posD = 1;
param2.pos = 1;
param2.lambda = LL2;
param2.iter = 4000;
param2.numThreads=16;
param2.batchsize=512;
    


%% Loading the renorm parameters (for stds1 and stds2)

load('sisec_code/saved_dicts/und/renorm_params_all.mat')


%% Loading the dicts

dataDirectory = 'sisec_code/saved_dicts/und';

dictionaries = {
%     'dict_non_vectorized_non_renorm_female_speaker_f1.mat'
%     'dict_non_vectorized_non_renorm_female_speaker_f2.mat'
%     'dict_non_vectorized_non_renorm_female_speaker_f3.mat'
%     'dict_non_vectorized_non_renorm_female_speaker_f4.mat'
    'dict_all_non_vectorized_non_renorm_male_speaker_m1.mat'
    'dict_all_non_vectorized_non_renorm_male_speaker_m2.mat'
    'dict_all_non_vectorized_non_renorm_male_speaker_m3.mat'
    'dict_all_non_vectorized_non_renorm_male_speaker_m4.mat'
%     'dict_non_vectorized_non_renorm_music_speaker_drums.mat'
%     'dict_non_vectorized_non_renorm_music_speaker_flute1.mat'
%     'dict_non_vectorized_non_renorm_music_speaker_flute2.mat'
%     'dict_non_vectorized_non_renorm_music_speaker_guitar1.mat'
%     'dict_non_vectorized_non_renorm_music_speaker_guitar2.mat'
%     'dict_non_vectorized_non_renorm_music_speaker_guitar3.mat'
};

Dnmf1_list = {};
Dnmf2_list = {};
for i=1:length(dictionaries)
    load(strcat(dataDirectory,'/' , dictionaries{i}));
    Dnmf1_list = [Dnmf1_list Dnmf1];
    Dnmf2_list = [Dnmf2_list Dnmf2];
end

clear Dnmf1 Dnmf2


%% Load the mix file from the dataset

%mix_file = '../data/sisec/und/test/test_female3_inst_mix.wav';
%mix_file = '../data/sisec/und/test/dev1_female3_inst_mix.wav';
mix_file = '../data/sisec/und/test/test_male4_inst_mix.wav';
[mix, Fs] = audioread(mix_file);
mix = resample(mix,fs,Fs);
mix = ((mix(:,1) + mix(:,2))/2)';
mix_len = length(mix);
mix = mix(1:min([mix_len,Npad]));

%% Creating my own mix

file_1 = '../data/sisec/und/female/f1/dev1_female4_src_1.wav';
file_2 = '../data/sisec/und/female/f2/dev1_female4_src_2.wav';
file_3 = '../data/sisec/und/female/f3/dev1_female4_src_3.wav';

% use this pair of files for a more challenging case
%file_1 = 'audio_files/bgat5a.wav'
%file_2 = 'audio_files/bgbh5s.wav'

% Load files and resample to fs = 16kHz
[x1, Fs] = audioread( file_1  );
x1 = resample(x1,fs,Fs);
x1 = x1(:)'; T1 = length(x1 );

[x2, Fs] = audioread( file_2 );
x2 = resample(x2,fs,Fs);
x2 = x2(:)'; T2 = length(x2);


[x3, Fs] = audioread( file_3 );
x3 = resample(x3,fs,Fs);
x3 = x3(:)'; T3 = length(x3);

% create the mixture
T = min([T1,T2,T3,Npad]);
x1 = x1(1:T);
x2 = x2(1:T);
x3 = x3(1:T);
mix = (x1+x2+x3);

%%
soundsc(x1,fs)
%%
soundsc(x2,fs)
%%
soundsc(x3,fs)
%%
soundsc(mix,fs)
%%
% perform separation

[speech_list, xest_list] = demix_scatt2top_multi(mix, Dnmf1_list, Dnmf2_list, stds1, stds2, eps, filts, scparam, param1, param2, Npad);


%% evaluate results
Parms_scatt2  =  BSS_EVAL(x1', x2', speech_list{1}', speech_list{1}', mix');
Parms_scatt1 =  BSS_EVAL(x1', x2', xest_list{1}', xest_list{2}', mix');

Parms_scatt2
Parms_scatt1

%%
soundsc(speech_list{1}, fs)
%%
soundsc(speech_list{2}, fs)
%%
soundsc(speech_list{3}, fs)
%%
soundsc(xest_list{1}, fs)
%%
soundsc(xest_list{2}, fs)
%%
soundsc(xest_list{3}, fs)
