function [ speaker ] = get_speaker( speakers, speaker_name )
%GET_SPEAKER Get Speaker from speaker struct

    for i=1:length(speakers)
        if strcmp(speakers(i).name, speaker_name) == 1
            speaker = speakers(i);
            break;
        end
    end

end

