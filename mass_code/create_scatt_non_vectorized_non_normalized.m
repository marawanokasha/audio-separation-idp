
%%

datasetFolder = '/home/marawan/workspace/matlab/idp/data/MASS/bearlin-roads (pop-rock)/bearlin-roads_0-28';

[sourcesL,sourcesR,sourceNameList,mixL,mixR,fs,nbits] = loadSources(datasetFolder,0);

% define scattering parameters and construct filters
eps=1e-3;
Npad = 1234801;
T = 2048;

scparam.N = Npad;
scparam.T = T;
scparam.Q = 32;
filts = create_scattfilters( scparam );

options.renorm=1;
options.parallel = 1;


dataDirectory = '/home/marawan/workspace/matlab/idp/data/sisec/und/';

%%

for i=1:size(speakers,2)
    
    X1 = [];
    X2 = [];
    stds1all = [];
    stds2all = [];
    display(strcat('Started Speaker: ', speakers(i).name))
    [X2, X1] = audioscatt_fwd_haar(speakers(i).fullTraining, filts, options);

    display(strcat('Finished Speaker: ', speakers(i).name))
%     save(strcat('saved_scatts/scatt_n500_X1_X2_', speakers(i).type, '_speaker_', speakers(i).name),'X1', 'X2', 'stds1all', 'stds2all')
    save(strcat('sisec_code/saved_scatts/und/scatt_all_non_vectorized_non_renorm_X1_X2_', speakers(i).type, '_speaker_', speakers(i).name),'X1', 'X2')
%     break
end
    