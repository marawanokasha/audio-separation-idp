function [ speaker_name ] = get_speaker_name_from_file( file_name )
%GET_SPEAKER_NAME_FROM_FILE Does string processing and gets the speaker
%name from the input file name

ind = strfind(file_name,'_');
speaker_name = file_name(ind(end)+1:strfind(file_name,'.mat')-1);

end

