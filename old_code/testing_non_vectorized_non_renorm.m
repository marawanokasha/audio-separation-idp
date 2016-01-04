

%% General parameters

% define scattering parameters and construct filters
eps=1e-3;
fs = 16000;
Npad = 2^15;
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

load('saved_dicts/renorm_params_s1_s7.mat')


%% Loading the dicts

load('saved_dicts/dict_n100_non_vectorized_non_renorm_male_speaker_s1.mat')

Dnmf11 = Dnmf1;
Dnmf21 = Dnmf2;
clear Dnmf1 Dnmf2

load('saved_dicts/dict_n100_non_vectorized_non_renorm_female_speaker_s7.mat')

Dnmf12 = Dnmf1;
Dnmf22 = Dnmf2;
clear Dnmf1 Dnmf2


%% Load the files to separate and create the mix


% use this pair of files for an example of good performance
file_1 = 'audio_files/sbba2n.wav';
file_2 = 'audio_files/sbbh2n.wav';

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


% create the mixture
T = min([T1,T2,Npad]);
x1 = x1(1:T);
x2 = x2(1:T);
mix = (x1+x2);

%%
soundsc(x1,fs)
%%
soundsc(x2,fs)
%%
soundsc(mix,fs)
%%
% perform separation
[speech1, speech2, xest1, xest2] = demix_scatt2top(mix, Dnmf11, Dnmf12, Dnmf21, Dnmf22, stds1, stds2, eps, filts, scparam, param1, param2, Npad);


%% evaluate results
Parms_scatt2  =  BSS_EVAL(x1', x2', speech1(1:T)', speech2(1:T)', mix');
Parms_scatt1 =  BSS_EVAL(x1', x2', xest1(1:T)', xest2(1:T)', mix');

Parms_scatt2
Parms_scatt1

%%
soundsc(xest2, fs)
%%
soundsc(speech2, fs)
%%
soundsc(speech1, fs)
%%
soundsc(xest1, fs)
