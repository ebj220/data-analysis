%function DMD_SRP
% Author: Lauren Jones, Elliot Jennis, Hari Annamalai
% 
%function DMD_SRP(File,scale,sizeX,Reconstruct,usePercent,numModes,2ndProjFile)

function DMD_SRP(varargin)

maxi=1;

set(0,'DefaultFigureColormap',jet())

close all

if nargin ~= 0%~isempty(varargin)
    Folder = [pwd '/'];
    File = varargin{1};
    Reconst = varargin{2};
    %scale = varargin{3};
    %sizeX = varargin{4}/scale;
    if Reconst
        UsePercent = varargin{3};                         %Use pecentage of total energy when reconstructing
        num_modes = varargin{4};
        if length(varargin)==7
            SecondProjFile = varargin{5};
        end
    else
        UsePercent = 0;
        num_modes = varargin{3};
    end
else
    [File,Folder]=uigetfile('./*.mat', 'Pick dataset: ');
    Reconst = 1;
    
    %scale = 3;
    %sizeX = 45/scale;                %P5=125, P1=45, P05=50?
end

load(strcat(Folder, File));

%if ~isempty(goodImgs)
%    S = S(:,:,goodImgs);
 %   N = size(S,3);
%end


if ~isempty(varargin)
    Name = [];
else
    Name = input('Save as (default: same as input): ');
end

if ~isempty(varargin)
    if Reconst && numel(varargin) == 6
        saveFolder = varargin{6};
        rnames = {Name};
    elseif ~Reconst && numel(varargin) == 4
        saveFolder = varargin{4};
        rnames = {Name};
    else
        
        if isempty(Name)
            saveFolder = [Folder strtok(File,'.') ' Results/POD'];
            rnames = {['./' strtok(File,'.') ' Results/POD']};
        else
            saveFolder = [Folder Name];
            rnames = {[Folder Name]};
        end
    end
else
    if isempty(Name)
        saveFolder = [Folder strtok(File,'.') ' Results/POD'];
        rnames = {['./' strtok(File,'.') ' Results/POD']};
    else
        saveFolder = [Folder Name];
        rnames = {[Folder Name]};
    end
end

if exist(saveFolder, 'dir')~=7
    mkdir(saveFolder)
end

if exist([saveFolder '/MatFigs/'], 'dir')~=7
    mkdir([saveFolder '/MatFigs/'])
end

xmax = 2.15;
minx = round(Stab.TipXfit-xmax*NozWidth);
if minx<1
   indx = 1;
else
   indx = minx;
end
U = double(S(round(Stab.TipYfit-1.9*NozWidth):round(Stab.TipYfit+1.9*NozWidth),indx:round(Stab.TipXfit+NozWidth),:));
Umean = nanmean(U,3);

%implay(framesmt)

%% Find bounds of PCA
figure
clf
pcolor(Umean)
shading interp
axis('equal', 'tight'); title('Mean Image of Data')

%PrepData

[X, Y] = meshgrid(1:size(U,2),1:size(U,1));
Uf = U-repmat(nanmean(U,3),[1,1,size(U,3)]);
Umean = nanmean(U,3);

%% create matrix will all fluctuating velocity components for each snapshot in a column
[uSize] = size(Uf);
Uall=reshape(Uf,uSize(1)*uSize(2),uSize(3));

%Do SVD analysis
[u,s,v]=svd(Uall.','econ');
figure
plot(diag(s)/(sum(diag(s))),'ro')
xlabel('Modes')
A = diag(s)/sum(diag(s));
energy_total = zeros(1,length(A));
energy_total(1) = A(1);
for i=2:length(A)
    energy_total(i) = A(i)+energy_total(i-1);
end
figure
subplot(3,1,1),yyaxis left, plot(diag(s)/sum(diag(s)),'ko','Linewidth',[2]), ylabel('% of Total Energy'),...
    yyaxis right, plot(energy_total,'--'), xlabel('Modes'), ylabel('Integrated Energy'),...
    ylim([0 1])
subplot(3,1,2), plot(u(1:uSize(3),1:3),'Linewidth',[2]);title('Dominant Three Modes');xlabel('Frames');  %modes
subplot(3,1,3), plot(v(1:uSize(3),1:3),'Linewidth',[2]); title('Time Dynamics');xlabel('Frames')  %time dynamics
%subplot(4,1,2), plot(t,v(:,1)/max(v:,1)),t,v(:,2)/max(v(:,2)),'Linewidth',[2])
%subplot(4,1,3), plot(x,u(:,1)/max(u(:,1)),'Linewidth',[2])
%subplot(4,1,4), plot(x,u(:,2)/max(u(:,2)),'Linewidth',[2])

%DMD
t = linspace(0,length(goodImgs)/76000,length(goodImgs));
dt = t(2)-t(1);
X = Uall; 
X1 = X(:,1:end-1); 
X2 = X(:,2:end);
r=6;
[U,S,V] = svd(X1,'econ'); 
Ur=U(:,1:r); Sr=S(1:r,1:r); Vr=V(:,1:r);
    
Atilde = Ur'*X2*Vr/Sr;    
[W,D] = eig(Atilde);    %eigen decomp of that matrix    
Phi = X2*Vr/Sr*W;   %DMD modes
    
lambda = diag(D);   %DMD eigenvalues
omega = log(lambda)/dt; %DMD frequencies
%DMD modes are the purple and gold projected on old plots

x1=X(:,1);
b = Phi\x1;  % pseudo-inverse initial conditions
u_modes = zeros(r,length(t));   %time dynamics-r rows, t columns

for iter = 1:length(t)
    %outer product
     u_modes(:,iter) =(b.*exp(omega*t(iter)));  %for each moment in time, compute this(b*e^iwt)
end
figure
plot(real(u_modes')),xlabel('Frames'),legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6')
%xlim([0 50])
X_dmd = Phi*u_modes; %DMD solution
for i = 1:r
figure
clf
%size(Umean)
%size(Phi)
size(X)
pcolor(reshape(real(Phi(:,i)),size(RMSCrop)))
shading interp
colorbar('SouthOutside')
axis('equal', 'tight'); title(['DMD Mode # ', num2str(i)])
end
DMDfreq = omega/(2*pi);
alpha = X'*Phi; %projection of u onto phi (from Lusseryan)
DMDpower = vecnorm(alpha);
figure
plot(abs(DMDfreq),-10*log10(DMDpower/max(DMDpower)),'*','LineWidth',2);
xlabel('Frequency (Hz)')
ylabel('- Power (dB)')
grid on;

figure
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--');
hold on, grid on
scatter(real(lambda),imag(lambda),'ok');
xlabel('Real')
ylabel('Imaginary')
