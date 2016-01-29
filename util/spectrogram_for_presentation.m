%%
x = speakers(1).fullTraining(:,1);

%%
soundsc(x, 16000)
%% This is what was used for producing the figures in the presentation
figure;
spectrogram(x,100,0,1000,16000,'yaxis')

%% Wavelet Transform example
figure;
coefs = cwt(x,1:32, 'cgau4');
sc = wscalogram('image',coefs);


%% STFT Magnitude Spectrum

figure;
sc = spectrogram(x,1024,512,1024,16000,'yaxis');


A = abs(sc);

plot(A)
