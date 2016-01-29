

%% read data

load('saved_scatts/scatt_n500_non_vectorized_X1_X2_female_speaker_s7.mat')

%%
data.X1 = X1;
data.X2 = X2;


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
save('saved_dicts/dict_n500_non_vectorized_female_speaker_s7.mat','Dnmf1','Dnmf2')

