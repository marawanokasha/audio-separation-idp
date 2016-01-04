%%

datasetFolder = '/home/marawan/workspace/matlab/idp/data/MASS/bearlin-roads (pop-rock)/bearlin-roads_0-28';

[sourcesL,sourcesR,sourceNameList,mixL,mixR,fs,nbits] = loadSources(datasetFolder,0);

%%
soundsc(sourcesL(:,2),fs)
%%
soundsc(mixL,fs)