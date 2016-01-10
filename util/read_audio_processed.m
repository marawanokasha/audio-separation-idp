function [ file_content ] = read_audio_processed( file_path, resample_freq, Npad )
%READ_AUDIO_PROCESSED reads audio and applies resampling, mono-channeling
%and padding as necessary and returns a column vector

[file_content, original_freq] = audioread(file_path);
if size(file_content,2) > 1
    file_content = ((file_content(:,1) + file_content(:,2))/2);
end
file_content = resample(file_content,resample_freq,original_freq);
file_content = pad_mirror(file_content,Npad);

end

