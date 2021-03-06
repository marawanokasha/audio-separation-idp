function [out] = audioreconstruct(in, data ,phase)

options = data.scparam;
os = getoptions(options,'os',1);
Q = getoptions(options,'Q',32);
T = getoptions(options,'T',4096);

[N,M, L]=size(in);
[N,M0, L]=size(phase);

%1) upscale modulus 
filts = data.filts;
J=size(filts.psi,2);
out=zeros(M0,L);
dse = getoptions(options,'dse',round(T*2^(-floor(J/Q)-os)));
%xx=repmat([1:N]',1,M);
%yy=repmat([1:dse:M0],N,1);
%yy0=repmat([1:M0],N,1);

[xx0,yy0]=meshgrid(1:dse:M0,1:N);
[xx,yy]=meshgrid(1:M0,1:N);

for l=1:L
inpu=interp2(xx0,yy0,in(:,:,l),xx,yy,'spline');

%2) use mix phase
coeffs=inpu.*phase(:,:,l);

%3) pseudoinverse
coeffs=transpose(coeffs);
for j=1:J
out(:,l) = out(:,l) + real(ifft(fft(coeffs(:,j)).*filts.dpsi{j}));
end
out(:,l) = out(:,l) + real(ifft(fft(coeffs(:,J+1)).*filts.dphi));

end



