% % Part 2 Comfortably Numb Bass
clear all; clc;
[y2, Fs2] = audioread('Floyd.m4a');
time2 = length(y2)/Fs2; % record time in seconds
L=15; n=length(y2)/time2*L;
t2=linspace(0,L,n+1); t=t2(1:n);
k=(1/L)*[0:n/2-1 -n/2:-1]; ks=fftshift(k);
% Bandpass filter
y_filter = bandpass(y2,[75,175],Fs2);
ygt_spec2=[];
tslide=0:0.03:t(end);
for j=1:length(tslide)
g=0.07*exp(-300*(t-tslide(j)).^2)'; % Gabor window
yg2=g.*y_filter(1:n); %GNR in Gabor Window time domain
ygt2=fft(yg2)'; %fft of gabor window GNR
ygt_spec2=[ygt_spec2; fftshift(abs(ygt2))];
end

figure(1)
pcolor(tslide,ks(330751:334976),log(abs(ygt_spec2(:,330751:334976).'+1))), shading interp
set(gca,'Fontsize',[14]);
colormap(hot)
axis([0 L 75 200])
xlabel('Time (s)','Fontsize',14), ylabel('Frequency (Hz)','Fontsize',14);
title('Comfortably Numb Bass','FontSize',18);
 hold on
 yline(82.31,'w','E2','FontSize',16),hold on,yline(92.5,'w','F#2','FontSize',16),hold on,
 yline(98.00,'w','G2','FontSize',16), hold on, yline(110,'w','A2','FontSize',16),hold on, 
 yline(123.47,'w','B2','FontSize',16), hold on, 
 yline(146.83,'w','D3','FontSize',16), hold on, yline(164.81,'w','E3','FontSize',16);