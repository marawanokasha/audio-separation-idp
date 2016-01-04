function [sourcesL,sourcesR,sourceNameList,mixL,mixR,fs,nbits] = loadSources(datasetFolder,effects,timeRange)

% Loads a data set from the MTG MASS Resources (http://www.mtg.upf.edu/static/mass/resources)
% The left channel and right channel tracks of a stereo song as well as their mixture are 
% stored in different matlab matrices.
%
% Usage: [sourcesL,sourcesR,sourceNameList,mixL,mixR,fs,nbits]
%        = loadSources(datasetFolder,effects)
%        [sourcesL,sourcesR,sourceNameList,mixL,mixR,fs,nbits] 
%        = loadSources(datasetFolder,effects,timeRange)
%
% Input:
%   - datasetFolder: folder where the data set to be loaded was uncompressed.
%   - effects: 0 (load sources without effects), 1 (load sources with effects).
%   - timeRange: if set, only the timeRange in seconds (e.g. [0 4]) from each track is loaded.
%
% Ouput:
%   - sourcesL: matrix with the sources of the left channel as column vectors.
%   - sourcesR: matrix with the sources of the right channel as column vectors.
%   - sourceNameList, list of structs where "sourceNameList(i).name" stores the name of source i.
%   - mixL: column vector of the left channel of the mixture.
%   - mixR: column vector of the right channel of the mixture.
%   - fs: sampling rate for all sources. If it is not constant among all sources, an error is thrown.
%   - nbits: number of bits per sample for all sources. If it is not constant among all sources, an error is thrown.
%
% Author: MarC Vinyes Raso (Music Technology Group - Universitat Pompeu Fabra) contact: m.drwx.org
% Last update: 29/11/2008

sourcesL=[];
sourcesR=[];
sourceNameList=[];
mixL=[];
mixR=[];

if (nargin<2) error('Missing input parameters. Usage: [sourcesL,sourcesR,sourceNameList,mixL,mixR] = loadSources(datasetFolder,effects)'); end

if (nargin==3 && ~(size(timeRange,1)==1 && size(timeRange,2)==2))
    error('Third input argument (timeRange) should be a row vector of two columns.')
end

if (nargin==3 && timeRange(1)>timeRange(2))
    error('Third input argument (timeRange)''s "start time" should be lesser than the "end time".')
end

if (~isdir(datasetFolder))
	error(['Folder "' datasetFolder '" doesn''t exist.']);
end

if (effects)
    effects=1;
	subFolder='tracks_with_effects';
else
    effects=0;
	subFolder='tracks_without_effects';
end

sourcesFolder=fullfile(datasetFolder,subFolder);

disp(['Loading data from folder: "' sourcesFolder '"...']);

if (~isdir(sourcesFolder))
	warning(['Subfolder "' subFolder '" is empty. Trying with effects=' num2str(1-effects) '.']);  
  effects=1-effects;
  if (effects)
    effects=1;
    subFolder='tracks_with_effects';
  else
    effects=0;
    subFolder='tracks_without_effects';
  end
  
  sourcesFolder=fullfile(datasetFolder,subFolder);

  disp(['Loading data from folder: "' sourcesFolder '"...']);

  if (~isdir(sourcesFolder))
    error('Wrong database folder.');
  end
end

sourceNameList=dir(fullfile(sourcesFolder,'*.wav'));
fs=0;
nbits=0;
numChannels=0;

if (size(sourceNameList,1)==0)
	error(['No .wav files were found in folder "' sourcesFolder '". Please check that the data set was uncompressed successfully.']);
end

for i=1:size(sourceNameList,1)
    %strip the .wav file extension
	filename=sourceNameList(i).name;
	l=size(filename,2);
	while (filename(l)~='.')
		l=l-1;
    end
    filename=filename(1:l-1);
    
    disp(['Loading source ' num2str(i) ': ' filename '...']);
    
	[source,fsnow,nbitsnow]=wavread(fullfile(sourcesFolder,sourceNameList(i).name));
	if (fs==0)
		fs=fsnow;
	else
		if (fs~=fsnow) 
            error('Sources have different sample rate.');
        end
	end
	
	if (nbits==0)
		nbits=nbitsnow;
    else
		if (nbits~=nbitsnow) 
            error('Sources have different number of bits per sample.');
        end
	end
	
	
	sourcesL=[sourcesL,source(:,1)];
    
    if (numChannels==0)
        numChannels=size(source,2);
    elseif (numChannels~=size(source,2));
            error(['Some sources of the song ' sourcesFolder ' are stereo and others mono. Please make them all stereo or mono.']);
    end
        
	if (numChannels>1) 
        sourcesR=[sourcesR,source(:,2)];
    end
    
	sourceNameList(i).name=filename;
    
end

disp('Creating mix...');

mixL=sum(sourcesL,2);
mixR=sum(sourcesR,2);

if (nargin<3)
    timeRange=[1 size(mixL,1)];
else
    timeRange=round(timeRange*fs+1);
    if (timeRange(2)>size(mixL,1))
        timeRange(2)=size(mixL,1);
        warning('End time, greater than song length, automatically corrected to the end of the song.');
    end
        if (timeRange(1)<0)
        timeRange(1)=1;
        warning('Negative start time automatically corrected to 0.');
    end
end

if (size(sourceNameList,1)>0)
    sourcesL=sourcesL(timeRange(1):timeRange(2),:);
    mixL=mixL(timeRange(1):timeRange(2));
    
    if (numChannels>1)
        sourcesR=sourcesR(timeRange(1):timeRange(2),:);
        mixR=mixR(timeRange(1):timeRange(2));  
    end
end