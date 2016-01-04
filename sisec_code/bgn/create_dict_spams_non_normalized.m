

%% read data


dataDirectory = 'sisec_code/saved_scatts/bgn/';

bigX1 = [];
bigX2 = [];
allFiles = dir(dataDirectory);
scatts = [];
for i=1:length(allFiles)
    if allFiles(i).isdir == false
        load(strcat(dataDirectory,'/' , allFiles(i).name))
        reshapedX1 = reshape(X1, [size(X1,1), size(X1,2) * size(X1,3)]);
        reshapedX2 = reshape(X2, [size(X2,1), size(X2,2) * size(X2,3)]);
        bigX1 = [bigX1 reshapedX1];
        bigX2 = [bigX2 reshapedX2];
        fileName = allFiles(i).name;
        identifier = fileName(size('scatt_all_non_vectorized_non_renorm_X1_X2_',2)+1:strfind(fileName,'.mat')-1)
        scatt.identifier = identifier;
        scatt.X1 = reshapedX1;
        scatt.X2 = reshapedX2;
        scatts = [scatts scatt];
    end
end


%% Compute stds

% get the stds of the combined X1 of all speakers and combined X2 of all
% speakers

stds1 = std(bigX1,0,2);
stds2 = std(bigX2,0,2);


%%
save('sisec_code/saved_dicts/bgn/renorm_params_all.mat', 'stds1', 'stds2')

%%
load('sisec_code/saved_dicts/bgn/renorm_params_all.mat')


%%

scatts = [];
load('sisec_code/saved_scatts/bgn/scatt_all_non_vectorized_non_renorm_X1_X2_male_speaker_m1.mat')
reshapedX1 = reshape(X1, [size(X1,1), size(X1,2) * size(X1,3)]);
reshapedX2 = reshape(X2, [size(X2,1), size(X2,2) * size(X2,3)]);
scatt.X1 = reshapedX1;
scatt.X2 = reshapedX2;
scatt.identifier = 'male_speaker_m1';
scatts = [scatts scatt];
load('sisec_code/saved_scatts/bgn/scatt_all_non_vectorized_non_renorm_X1_X2_male_speaker_m2.mat')
reshapedX1 = reshape(X1, [size(X1,1), size(X1,2) * size(X1,3)]);
reshapedX2 = reshape(X2, [size(X2,1), size(X2,2) * size(X2,3)]);
scatt.X1 = reshapedX1;
scatt.X2 = reshapedX2;
scatt.identifier = 'male_speaker_m2';
scatts = [scatts scatt];
load('sisec_code/saved_scatts/bgn/scatt_all_non_vectorized_non_renorm_X1_X2_male_speaker_m3.mat')
reshapedX1 = reshape(X1, [size(X1,1), size(X1,2) * size(X1,3)]);
reshapedX2 = reshape(X2, [size(X2,1), size(X2,2) * size(X2,3)]);
scatt.X1 = reshapedX1;
scatt.X2 = reshapedX2;
scatt.identifier = 'male_speaker_m3';
scatts = [scatts scatt];

%% Running dictionary learning for every speaker

for i=1:length(scatts)

    data.X1 = scatts(i).X1;
    data.X2 = scatts(i).X2;
    
    Npad = 2^17;

    options.renorm=1;
    if options.renorm
       %renormalize data: whiten each frequency component.
       eps=2e-3;
       data.X1 = renorm_spect_data(data.X1, stds1, eps);

       eps=1e-3;
       data.X2 = renorm_spect_data(data.X2, stds2, eps);
    end
    

    %%% NMF for X1
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

    %%% NMF for X2
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

    save(strcat('sisec_code/saved_dicts/bgn/dict_all_non_vectorized_non_renorm_', scatts(i).identifier, '.mat'),'Dnmf1','Dnmf2')

end
