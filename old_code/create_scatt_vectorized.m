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