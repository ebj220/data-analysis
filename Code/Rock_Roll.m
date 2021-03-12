%Rock&Roll

clc; clear all;

%%Part 1 Sweet Child'O Mine
[y1, Fs1] = audioread('GNR.m4a');
time1 = length(y1)/Fs1; % record time in seconds

L=5; %domain size (could be 10 seconds, hours, etc)
n=length(y1)/time1*L;
t2=linspace(0,L,n+1); 
t=t2(1:n);
k=(1/L)*[0:n/2-1 -n/2:-1]; %frequency Hz
ks=fftshift(k);

%% Gabor Filter
ygt_spec1=[]; 
tslide=0:0.05:t(end);
for j=1:length(tslide)
g=.05*exp(-5000*(t-tslide(j)).^2); % Gabor window size
yg1=g.*transpose(y1(1:n)); 
ygt1=fft(yg1); 
ygt_spec1=[ygt_spec1; fftshift(abs(ygt1))]; 
end

figure(1)
pcolor(tslide,ks,log(abs(ygt_spec1).'+1));
shading interp
set(gca,'Fontsize',[20]);
xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14),title('Sweet Child Of Mine Guitar','FontSize',18)
axis([0 L 0 1000]);
colormap(hot);
hold on, yline(277.18,'w','C#4','Fontsize',16), hold on, yline(311.13,'w','D#4','Fontsize',16),hold on
yline(369.99,'w','F#4','Fontsize',16), hold on, yline(415.3,'w','G#4','Fontsize',16), hold on, yline(554.37,'w','C#5','Fontsize',16);
hold on, yline(698.46,'w','F5','Fontsize',16), hold on, yline(739.99,'w','F#5','Fontsize',16);

% % % Part 2 Comfortably Numb Bass
% clear all; clc;
% [y2, Fs2] = audioread('Floyd.m4a');
% time2 = length(y2)/Fs2; % record time in seconds
% %
% %p8 = audioplayer(y2,Fs2); playblocking(p8);
% L=15; n=length(y2)/time2*L;
% t2=linspace(0,L,n+1); t=t2(1:n);
% k=(1/L)*[0:n/2-1 -n/2:-1]; ks=fftshift(k);
% % Bandpass filter
% y_filter = bandpass(y2,[75,175],Fs2);
% ygt_spec2=[];
% tslide=0:0.03:t(end);
% for j=1:length(tslide)
% g=0.07*exp(-300*(t-tslide(j)).^2)'; % Gabor window
% yg2=g.*y_filter(1:n); %GNR in Gabor Window time domain
% ygt2=fft(yg2)'; %fft of gabor window GNR
% ygt_spec2=[ygt_spec2; fftshift(abs(ygt2))];
% % subplot(3,1,1), plot(t,y_filter(1:n),'k',t,g,'r')
% % subplot(3,1,2), plot(t,yg2,'k')
% % subplot(3,1,3), plot(ks(n/2:n),fftshift(abs(ygt2(n/2:n)))/max(abs(ygt2(n/2:n))))
% % xlim([0 250])
% % drawnow
% % pause(0.01)
% end
% %close all
% 
% figure(2)
%  pcolor(tslide,ks(330751:334976),log(abs(ygt_spec2(:,330751:334976).'+1))), shading interp
% % pcolor(tslide,ks,log(abs(ygt_spec2).'+1));
% %shading interp
% set(gca,'Fontsize',[14]);
% colormap(hot)
% axis([0 L 75 200])
% xlabel('Time (s)','Fontsize',14), ylabel('Frequency (Hz)','Fontsize',14);
% title('Comfortably Numb Bass','FontSize',18);
%  hold on
%  yline(82.31,'w','E2','FontSize',16),hold on,yline(92.5,'w','F#2','FontSize',16),hold on,
%  yline(98.00,'w','G2','FontSize',16), hold on, yline(110,'w','A2','FontSize',16),hold on, 
%  yline(123.47,'w','B2','FontSize',16), hold on, 
%  yline(146.83,'w','D3','FontSize',16), hold on, yline(164.81,'w','E3','FontSize',16);
% % 

% % Part 3 Comfortably Numb Guitar
% clear all; clc
% [y2, Fs2] = audioread('Floyd.m4a');
% time2 = length(y2)/Fs2; % record time in seconds
% %
% %p8 = audioplayer(y2,Fs2); playblocking(p8);
% L=10; n=length(y2)/time2*L;
% t2=linspace(0,L,n+1); t=t2(1:n);
% k=(1/L)*[0:n/2-1 -n/2:-1]; ks=fftshift(k);
% % Bandpass filter
% y_filter = bandpass(y2,[200,700],Fs2);
% ygt_spec2=[];
% tslide=0:0.03:t(end);
% for j=1:length(tslide)
% g=0.1*exp(-450*(t-tslide(j)).^2)'; % Gabor window
% yg2=g.*y_filter(1:n); %GNR in Gabor Window time domain
% ygt2=fft(yg2)'; %fft of gabor window GNR
% ygt_spec2=[ygt_spec2; fftshift(abs(ygt2))];
% % subplot(3,1,1), plot(t,y_filter(1:n),'k',t,g,'r')
% % subplot(3,1,2), plot(t,yg2,'k')
% % subplot(3,1,3), plot(ks(n/2:n),fftshift(abs(ygt2(n/2:n)))/max(abs(ygt2(n/2:n))))
% % xlim([0 250])
% % drawnow
% % pause(0.01)
% end
% %close all
% 
% figure(3)
% %pcolor(tslide,ks(330751:334976),log(abs(ygt_spec2(:,330751:334976).'+1))), shading interp
% pcolor(tslide,ks,log(abs(ygt_spec2).'+1));
% shading interp
% set(gca,'Fontsize',[14]);
% colormap(hot)
% axis([0 L 200 750])
% xlabel('Time (t)','Fontsize',14), ylabel('frequency (Hz)','Fontsize',14);
% title('Comfortably Numb Guitar','FontSize',18);
% hold on, yline(293.66,'w','D4','Fontsize',16), hold on, yline(329.63,'w','E4','Fontsize',16),hold on
% yline(369.99,'w','F#4','Fontsize',16), hold on, yline(392,'w','G4','Fontsize',16), hold on
% yline(415.3,'w','G#4','Fontsize',16),hold on 
% yline(440,'w','A4','Fontsize',16),hold on, yline(493.88,'w','B4','Fontsize',16), hold on
% yline(554.37,'w','C#5','Fontsize',16),hold on, yline(587.33,'w','D5','Fontsize',16); 
% 

