
%%

% define scattering parameters and construct filters
eps=1e-3;
fs = 16000;
Npad = 2^15;
T = 2048;

scparam.N = Npad;
scparam.T = T;
scparam.Q = 32;
filts = create_scattfilters( scparam );

options.renorm=1;
options.parallel = 1;

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
dataDirectory = '/media/Work/workspace/matlab/idp/data';
%speakers = read_audio(dataDirectory, 500, 200, Npad);

%%
%save('saved_dicts/speakers.mat', 'speakers')

%%
load('saved_dicts/speakers.mat')


%%
parpool('local',2); 

%%

for i=1:size(speakers,2)
    
    X1 = [];
    X2 = [];
    stds1all = [];
    stds2all = [];
    display(strcat('Started Speaker ', speakers(i).name))
    for j=1: size(speakers(i).fullTraining,2)
        trainingSample = speakers(i).fullTraining(:,j);
        [X2Sample, X1Sample] = audioscatt_fwd_haar(trainingSample, filts, options);
        stds1 = std(X1Sample,0,2);
        stds2 = std(X2Sample,0,2);
        
        options.renorm=1;
        if options.renorm
           %renormalize data
           eps=2e-3;
           X1Sample = renorm_spect_data(X1Sample, stds1, eps);

           eps=1e-3;
           X2Sample = renorm_spect_data(X2Sample, stds2, eps);
        end

        X1Sample = reshape(X1Sample, [size(X1Sample,1) * size(X1Sample,2),1]);
        X2Sample = reshape(X2Sample, [size(X2Sample,1) * size(X2Sample,2),1]);
        
        X1(:,j) = X1Sample;
        X2(:,j) = X2Sample;

        stds1all(:,j) = stds1;
        stds2all(:,j) = stds2;
        
    end
    display(strcat('Finished Speaker ', speakers(i).name))
    save(strcat('saved_scatts/scatt_n500_X1_X2_', speakers(i).type, '_speaker_', speakers(i).name),'X1', 'X2', 'stds1all', 'stds2all')

end
    
%%
%speakers2 = read_audio(dataDirectory, 2, 10, Npad);
[XX2, XX1] = audioscatt_fwd_haar(speakers2(1).fullTraining, filts, options);

%%
% parpool('local',2); 
[X2, X1] = audioscatt_fwd_haar(speakers(1).fullTraining, filts, options);    

%%
X1 = reshape(X1, [size(X1,1), size(X1,2)* size(X1,3)]);
X2 = reshape(X2, [size(X2,1), size(X2,2)* size(X2,3)]);



%%
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
%compute renormalization parameters    
if ~exist('stds1','var')
stds1 = std(X1,0,2);
stds2 = std(X2,0,2);
end


%%

data.X1 = X1;
data.X2 = X2;
% this code assumes that data contains in the fields X1 and X2 the first and second level scattering coefficients for some training data.  
% these coefficients can be computed as:
% [X2, X1] = audioscatt_fwd_haar(pad_mirror(x',Npad), filts, options);
% where x is the audio signal in the time domain

    
options.renorm=1;
if options.renorm
   %renormalize data: whiten each frequency component.
   eps=2e-3;
   data.X1 = renorm_spect_data(data.X1, stds1, eps);
        
   eps=1e-3;
   data.X2 = renorm_spect_data(data.X2, stds2, eps);
end
   
    
%% train models
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
save('saved_dicts/matlab_first_trial.mat','Dnmf1','Dnmf2')


%%
XS1 = data.X1;
XS2 = data.X2;
save('saved_scatts/scatt_dim23_n30.mat','XS1','XS2')
