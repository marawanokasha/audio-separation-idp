

%% read data

load('saved_scatts/scatt_n100_non_vectorized_non_renorm_X1_X2_male_speaker_s1.mat')
X11 = X1;
X21 = X2;
clear X1 X2

load('saved_scatts/scatt_n100_non_vectorized_non_renorm_X1_X2_female_speaker_s7.mat')
X12 = X1;
X22 = X2;
clear X1 X2


%% Compute stds

% get the stds of the combined X1 of all speakers and combined X2 of all
% speakers
X1 = [reshape(X11, [size(X11,1), size(X11,2) * size(X11,3)]) reshape(X12, [size(X12,1), size(X12,2) * size(X12,3)]) ];
X2 = [reshape(X21, [size(X21,1), size(X21,2) * size(X21,3)]) reshape(X22, [size(X22,1), size(X22,2) * size(X22,3)]) ];

stds1 = std(X1,0,2);
stds2 = std(X2,0,2);

%%
save('saved_dicts/renorm_params_s1_s7.mat', 'stds1', 'stds2')

%% Running dictionary learning for speaker 1

data.X1 = reshape(X11, [size(X11,1), size(X11,2) * size(X11,3)]);
data.X2 = reshape(X21, [size(X21,1), size(X21,2) * size(X21,3)]);

%%

Npad = 2^15;
    
options.renorm=1;
if options.renorm
   %renormalize data: whiten each frequency component.
   eps=2e-3;
   data.X1 = renorm_spect_data(data.X1, stds1, eps);
        
   eps=1e-3;
   data.X2 = renorm_spect_data(data.X2, stds2, eps);
end

%% train modelss
model = 'NMF-scatt2';
    
    
%%%%Plain NMF%%%%%%%
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

% need spams to compute this dictionary, but any algorithm for NMF would do     
Dnmf1 = mexTrainDL(abs(data.X1),param1);
    
%%
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
    
% need spams to compute this dictionary, but any algorithm for NMF would do     
Dnmf2 = mexTrainDL(abs(data.X2),param2);
   
%%
save('saved_dicts/dict_n100_non_vectorized_non_renorm_male_speaker_s1.mat','Dnmf1','Dnmf2')

%% Running dictionary learning for speaker 1

data.X1 = reshape(X12, [size(X12,1), size(X12,2) * size(X12,3)]);
data.X2 = reshape(X22, [size(X22,1), size(X22,2) * size(X22,3)]);

%%

Npad = 2^15;
    
options.renorm=1;
if options.renorm
   %renormalize data: whiten each frequency component.
   eps=2e-3;
   data.X1 = renorm_spect_data(data.X1, stds1, eps);
        
   eps=1e-3;
   data.X2 = renorm_spect_data(data.X2, stds2, eps);
end

%% train modelss
model = 'NMF-scatt2';
    
    
%%%%Plain NMF%%%%%%%
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

% need spams to compute this dictionary, but any algorithm for NMF would do     
Dnmf1 = mexTrainDL(abs(data.X1),param1);
    
%%
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
    
% need spams to compute this dictionary, but any algorithm for NMF would do     
Dnmf2 = mexTrainDL(abs(data.X2),param2);
   
%%
save('saved_dicts/dict_n100_non_vectorized_non_renorm_female_speaker_s7.mat','Dnmf1','Dnmf2')



