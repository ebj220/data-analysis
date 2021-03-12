%Elliot Jennis
%AMATH 582

clc;clear all; close all;
load subdata.mat; % Imports the data as a 262144x49 (space by time) matrix called subdata 5
%% Initialzing data resolution and relevent variables

L = 10; % spatial domain
n = 64; % Fourier modes
x2 = linspace(-L,L,n+1); x = x2(1:n); y =x; z = x;
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; ks = fftshift(k);

[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);
Utave = zeros(n,n,n);

%% Averaging Algorithm

for j=1:49
    
Un(:,:,:)=reshape(subdata(:,j),n,n,n); %spacial matrix
Ut = fftn(Un); %transform of data
Uf = fftshift(Ut);  %shifted data
Utave = Utave + Uf;  %averaging to remove noise

end
Mave = max(abs(Utave),[],'all');

isosurface(Kx,Ky,Kz,abs(Utave)/Mave,0.7) 
axis([-10 10 -10 10 -10 10]), grid on;
xlabel('Kx');
ylabel('Ky');
zlabel('Kz');
hold on;

%% Data filtering

Xfreq = 5;
Yfreq = -7;
Zfreq = 2;
band = .6;
%gaussian filter
gx=exp(-band*(Kx-Xfreq).^2);
gy=exp(-band*(Ky-Yfreq).^2);
gz=exp(-band*(Kz-Zfreq).^2);
filter = gx.*gy.*gz;
%plot of filter in signal domain
isosurface(Kx,Ky,Kz,filter,0.7); 
axis([-10 10 -10 10 -10 10]), grid on;
xlabel('Kx');
ylabel('Ky');
zlabel('Kz');

xloc = (1:49);
yloc = (1:49);
zloc = (1:49);

for j=1:49   
Un1(:,:,:)=reshape(subdata(:,j),n,n,n); %spacial matrix
Ut1 = fftshift(fftn(Un1));  %transform of data
Utf = Ut1.*filter;          %filtering data
Unf = ifftn(fftshift(Utf));
[M,Index] = max(abs(Unf(:))); %max value of filtered data and index of location
xloc(j) = X(Index);           %matricies of location data
yloc(j) = Y(Index);
zloc(j) = Z(Index); 
end

figure(2)
plot3(xloc,yloc,zloc,'b','Linewidth',3); %3d visualization of sub path
grid on
axis([-8 8 -8 8 -8 8]);
hold on
plot3(xloc(49),yloc(49),zloc(49),'*'); %final location of sub
xlabel('X'), ylabel('Y'), zlabel('Z');