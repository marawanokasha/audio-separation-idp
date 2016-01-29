
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
speakers = load_speakers(dataDirectory, 500, 200, fs, Npad);

%%
%save('saved_dicts/speakers_resampled.mat', 'speakers')

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
    
        X1 = [X1 X1Sample];
        X2 = [X2 X2Sample];
        
        if mod(j, 10) == 0
            display(j)
        end
%         if j > 3
%             break
%         end
    end
    stds1 = std(X1,0,2);
    stds2 = std(X2,0,2);

    if options.renorm
       %renormalize data
       eps=2e-3;
       X1 = renorm_spect_data(X1, stds1, eps);

       eps=1e-3;
       X2 = renorm_spect_data(X2, stds2, eps);
    end
    

    display(strcat('Finished Speaker ', speakers(i).name))
%     save(strcat('saved_scatts/scatt_n500_X1_X2_', speakers(i).type, '_speaker_', speakers(i).name),'X1', 'X2', 'stds1all', 'stds2all')
    save(strcat('saved_scatts/scatt_n500_non_vectorized_X1_X2_', speakers(i).type, '_speaker_', speakers(i).name),'X1', 'X2', 'stds1', 'stds2')
%     break
end
    