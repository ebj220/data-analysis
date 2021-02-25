%582 PCA Part 2 
clear all;
clc;

%% Camera 1
camera_1_data = load('cam1_1.mat');

camera_1_data.vidFrames1_2 = double(camera_1_data.vidFrames1_1);
[l w d numFrames1_1] = size(camera_1_data.vidFrames1_1);

for k = 1:numFrames1_1
    mov1_1(k).cdata = uint8(camera_1_data.vidFrames1_1(:,:,:,k));
    mov1_1(k).colormap = [];
end

for j = 1:numFrames1_1
   X_1_1 = rgb2gray(frame2im(mov1_1(j)));
   X_1bw = bwareaopen(X_1_1,50);
   X1(:,:,j) = double(X_1_1);
   %  imshow(X_1_2); drawnow;
end



%% Camera 2
camera_2_data = load('cam2_1.mat');
camera_2_data.vidFrames2_1 = double(camera_2_data.vidFrames2_1);
[l w d numFrames2_1] = size(camera_2_data.vidFrames2_1);

for k = 1:numFrames2_1
    mov2_1(k).cdata = uint8(camera_2_data.vidFrames2_1(:,:,:,k));
    mov2_1(k).colormap = [];
end
for j = 1:numFrames1_1
     X_2_1 = rgb2gray(frame2im(mov2_1(j)));
   X2(:,:,j) = double(X_2_1);
    % imshow(X_2_1); drawnow;
end


%% Camera 3
camera_3_data = load('cam3_1.mat');
camera_3_data.vidFrames3_1 = double(camera_3_data.vidFrames3_1);
[l w d numFrames3_1] = size(camera_3_data.vidFrames3_1);

for k = 1:numFrames3_1
    mov3_1(k).cdata = uint8(camera_3_data.vidFrames3_1(:,:,:,k));
    mov3_1(k).colormap = [];
end
for j = 1:numFrames1_1
    X_3_1  = rgb2gray(frame2im(mov3_1(j)));
   X3(:,:,j) = double(X_3_1);
%    imshow(X_3_2); drawnow;
end


%% Data Extraction
%Camera 1
XCrop = X1(200:480,:,:);
for n = 1:numFrames1_1
    if n ==1
    [row col] = find(XCrop(:,:,n)>=250); %find white spots
    for i = 1:length(row)
        row(i) = row(i) + 200;
    end
    Arow_1(n) = sum(row)/length(row);
    Acol_1(n) = sum(col)/length(col);
    else
        [row col]=find(X1(:,:,n)>=240);
        i = 1;
        while i <= length(row)
            if row(i)<(Arow_1(n-1)-30)|row(i)>(Arow_1(n-1)+30)
                row(i)=[];
                col(i)=[];
            else
                i = i+1;
            end
        end
        i = 1;
         while i <= length(row)
            if col(i)<(Acol_1(n-1)-30)|col(i)>(Acol_1(n-1)+30)
                row(i)=[];
                col(i)=[];
            else
                i = i+1;
            end
        end
        Arow_1(n) = sum(row)/length(row);
        Acol_1(n) = sum(col)/length(col);
    end
%         
%      %   scatter(Acol1(n),Arow1(n));
%         scatter(col,row);
%           axis([0 640 0 480]);
%         pause(.05)
% %       hold on;
end
%Camera 2
for n = 1:numFrames1_1
    if n ==1
    [row col] = find(X2(:,:,n)>=255); %find white spots
    Arow_2(n) = sum(row)/length(row);
    Acol_2(n) = sum(col)/length(col);
    else
        [row col]=find(X2(:,:,n)>=247);
        i = 1;
        while i <= length(row)
            if row(i)<(Arow_2(n-1)-30)|row(i)>(Arow_2(n-1)+30)
                row(i)=[];
                col(i)=[];
            else
                i = i+1;
            end
        end
        i = 1;
         while i <= length(row)
            if col(i)<(Acol_2(n-1)-20)|col(i)>(Acol_2(n-1)+20)
                row(i)=[];
                col(i)=[];
            else
                i = i+1;
            end
        end
        Arow_2(n) = sum(row)/length(row);
        Acol_2(n) = sum(col)/length(col);
    end
        
%      %   scatter(Acol1(n),Arow1(n));
%        scatter(col,row);
%        axis([0 640 0 480]);
%        pause(.05)
%        grid on;
%       hold on;
end
%Camera 3
for n = 1:numFrames1_1
    if n ==1 %First frame search
    [row col] = find(X3(:,:,n)>=255); %find white spots
    Arow_3(n) = sum(row)/length(row);
    Acol_3(n) = sum(col)/length(col);
    else
        [row col]=find(X3(:,:,n)>=235);
        i = 1;
        while i <= length(row) %Search location dictated by averaged value from previous search
            if row(i)<(Arow_3(n-1)-45)|row(i)>(Arow_3(n-1)+45)
                row(i)=[];
                col(i)=[];
            else
                i = i+1;
            end
        end
        i = 1;
         while i <= length(row)
            if col(i)<(Acol_3(n-1)-45)|col(i)>(Acol_3(n-1)+45)
                row(i)=[];
                col(i)=[];
            else
                i = i+1;
            end
        end
        Arow_3(n) = sum(row)/length(row);
        Acol_3(n) = sum(col)/length(col);
    end
        
%      %   scatter(Acol1(n),Arow1(n));
%        scatter(col,row);
%        axis([0 640 0 480]);
%        pause(.05)
%        grid on;
%       hold on;
end
% figure(1)
% subplot(3,1,1)
% scatter(Acol1,Arow1)
%  axis([0 640 0 480]);
% subplot(3,1,2)
% scatter(Acol2,Arow2)
%  axis([0 640 0 480]);
% subplot(3,1,3)
% scatter(Acol3,Arow3)
%  axis([0 640 0 480]);
 %% PCA Analysis
X = [Acol_1;Arow_1;Acol_2;Arow_2;Acol_3;Arow_3]; %data matrix
[m,n]=size(X); %size of X
Xave = mean(X,2);
Xm=X-repmat(Xave,1,n); %subtract mean of data
% scatter3(Xm(1,:),Xm(2,:),Xm(3,:));

[u,s,v]=svd(Xm'/sqrt(n-1),'econ'); % perform the SVD
Y=u.*X'; % produce the principal components projection
figure(2)
plot(Y(:,1));
hold on
plot(Y(:,2));
xlabel('Frames','Fontsize',14);
ylabel('Travel Distance','Fontsize',14);
title('Dominant Mode Projection','Fontsize',20);
axis([0 numFrames1_1 -50 50]);
legend('DOF 1', 'DOF 2');
%Modal Energies
figure(3)
semilogy(s/sum(s),'o','MarkerFaceColor',[0 0.447 0.741]);
grid on
xticks([1 2 3 4 5 6]);
xlabel('Modes','Fontsize',12);
ylabel('Modal Energy Percentage','Fontsize',14);
title('Modal Energy','Fontsize',20);
