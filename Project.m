clc
clear
close all 
%cd 'C:\Users\scrpa\OneDrive - Politecnico di Milano\Desktop\Poli\Magistrale\Primo anno\BSPMI\MI\project repo\MI-Project'
%% Animation on: a=1 Animation off: a=0;
a=1;
%% dati 
load MRIdata.mat

%% AXIAL PLANE: Pre-processing of MRI, segmentation of lesion and total volume

if (a==1)
    figure(1)
        montage(vol)
        title("MRI iniziale.")
    figure(2)
        montage(vol(:,:,64:90))
        title('Tumor From slice 64 to 90.')
end 

% the tumor is visible from 65 to 89 slice
%% 
% isolate the region of the lesion
figure()
subplot(2,1,1)
imshow(vol(:,:,75))
impixelinfo
title("Isolate the tumor")
subplot(2,1,2)
[Cropped_vol d]= imcrop(vol(:,:,75), [130 102 51 45]); %%% qualcosa di strano

% Dimensioni del taglio 
v1=round(d(2)):(round(d(2))+length(Cropped_vol(:,1)));
v2=round(d(1)):(round(d(1))+length(Cropped_vol(1,:)));
v3=64:90;
VOI=vol(v1,v2,v3);

figure('Name',"Confronto funzione matlab e dimensioni ottenute")
subplot(2,1,1)
imshow(VOI(:,:,(12)));
title('Immagine ricavata dalle dimensioni')
subplot(2,1,2)
imshow(Cropped_vol)
title('Immagine ricavata dalla funzione')

%% 
% studio histogrammi 
figure('Name', "Istogrammi")
for i=1:27
    subplot(2,1,1)
    imshow(VOI(:,:,i))
    colorbar
    subplot(2,1,2)
    histogram(VOI(:,:,i),255); 
    xlim([0,255])
    grid on 
    pause(1)
end 

%%
% manual correction of white region
MASK=VOI(:,:,25:27);
MASK(MASK>240)=0;
VOI(:,:,25:27)=MASK;
figure()
montage(VOI)


%Aumento del contrasto
gamma=[0.5,1,1.5,2];
gammas={'0.5','1.0','1.5','2.0'};
figure()

LOWin = 150/255;
HIGHin = max(max(VOI(:,:,i)))/255;

for z=1:length(gamma)
   
    for i=1:size(VOI,3)
        vol_imadjusted(:,:,i) = imadjust(VOI(:,:,i),[0 0.5882],[0 1],gamma(z));
    end 
    subplot(1,length(gamma),z)
    montage(vol_imadjusted)
    title(['Gamma =' gammas(z)])
end 
% scegliamo gamma --> 2 


figure('Name', "Istogrammi")
for i=1:27
    subplot(2,2,1)
        imshow(vol_imadjusted(:,:,i))
        colorbar
        title("Modified image")
    subplot(2,2,2)
        histogram(vol_imadjusted(:,:,i),255); 
        xlim([0,255])
        grid on 
    subplot(2,2,3)
        imshow(VOI(:,:,i))
        colorbar
    subplot(2,2,4)
        histogram(VOI(:,:,i),255); 
        xlim([0,255])
        xline(150,'r')
        grid on 
        title("Original image")
    pause(1)
end 


%%
%salt & pepper filtering
for i=1:length(v3)
    vol_pn(:,:,i)=medfilt2(vol_imadjusted(:,:,i), [6 6]);
end

figure()
montage(vol_pn)

%%
%Binarizzazione 
bin_vol=imbinarize(vol_pn,0.8);

figure()
subplot(1,2,1)
montage(vol_pn)
title('BEFORE BIN')
subplot(1,2,2)
montage(bin_vol)
title("AFTER BIN")

%%
% Prendo i contorni 

%vidfile = VideoWriter('testmovie.mp4','MPEG-4');
%open(vidfile);
figure()
j=0;
for i=2:25
    j=j+1;
    imshow(VOI(:,:,i))
    title("Contour of the tumor")
    hold on
    imcontour(bin_vol(:,:,i),5,'m');
    %F(j) = getframe(gcf); 
    %writeVideo(vidfile,F(j));
    pause (1)
end 
%close(vidfile)

%%
% area of the binarized image for every slice
Axial_num_pixel=[];
for i=2:25
    Axial_num_pixel=[Axial_num_pixel sum(sum(bin_vol(:,:,i)==1))]; %conta i pixel bianchi 
end 

% total volume of the lesion in mm^3
volume_l=0;
for i=1:24
     volume_l=volume_l+Axial_num_pixel(1,i).*pixdim(1,3).*pixdim(1,1).*pixdim(1,2);
end

%%
%3D Visualization 
volumeViewer(vol(v1,v2,v3))

%% SAGITTAL PLANE: Pre-processing of MRI, segmentation of lesion and total volume
%from axial to sagittal plane:
for i=1:dim(1)
    vol_s(:,:,i)=vol(i,:,:); %% perchè vol_ax se siamo sul piano sagittale?
end

%display images
if (a==1)
    figure(1)
        montage(vol_s)
        title("MRI Sagittale.")
    figure(2)
        montage(vol_s(:,:,107:144))
        title('Tumor From slice 107 to 144.')
end 

% the lesion is visible from slice 108 to 143
%%
[Cropped_vol_s d_s]= imcrop(vol_s(:,:,126), [60 138 30 38]);
%55 135 35 60
% Dimensioni del taglio 
v1=round(d_s(2)):(round(d_s(2))+length(Cropped_vol_s(:,1)));
v2=round(d_s(1)):(round(d_s(1))+length(Cropped_vol_s(1,:)));
v3=107:144;

VOI_s=vol_s(v1,v2,v3);
figure()
    montage(VOI_s)

%% Istogramma 
figure('Name', "Istogrammi")
for i=1:size(VOI_s,3)
    subplot(2,1,1)
    imshow(VOI_s(:,:,i))
    colorbar
    subplot(2,1,2)
    histogram(VOI_s(:,:,i),255); 
    xlim([0,255])
    grid on 
    pause(1)
end
%%
MASK=VOI_s(:,:,22:25);
MASK(MASK>230)=0;
VOI_s(:,:,22:25)=MASK;
figure()
montage(VOI_s)

%% Aumento del contrasto

 for i=1:size(VOI_s,3)
     VOI_adj(:,:,i) = imadjust(VOI_s(:,:,i),[0 0.5882],[0 1],2);
 end 

 figure()
 montage(VOI_adj)

%% salt & pepper filtering
for i=1:size(VOI_s,3)
    VOI_pn(:,:,i)=medfilt2(VOI_adj(:,:,i), [6 6]);
end

figure()
subplot(2,1,1)
    montage(VOI_s)
 subplot(2,1,2)
 montage(VOI_pn)


%% histogramm 

figure('Name', "Istogrammi")
for i=1:size(VOI_s,3)
    subplot(2,1,1)
    imshow(VOI_pn(:,:,i))
    colorbar
    subplot(2,1,2)
    histogram(VOI_pn(:,:,i),255); 
    xlim([0,255])
    grid on 
    pause(1)
end 

%% Binarizzazione 
bin_vol=imbinarize(VOI_pn,0.8);

figure()
subplot(1,2,1)
montage(VOI_pn)
title('Enhanced contrast')
subplot(1,2,2)
montage(bin_vol)
%%
%Prendo i contorni 
figure()
for i=1:size(bin_vol,3)
    imshow(VOI_s(:,:,i))
    hold on
    imcontour(bin_vol(:,:,i),4,'m')
    pause (1)
end 
title("Contours of the tumor")

%%
num_pixel_s=[];
for i=2:25
    num_pixel_s=[num_pixel_s sum(sum(bin_vol(:,:,i)==1))]; %conta i pixel bianchi 
end 

% total volume of the lesion in mm^3
volume_l_s=0;
for i=1:24
     volume_l_s=volume_l_s+num_pixel_s(1,i).*pixdim(1,3).*pixdim(1,1).*pixdim(1,2);
end

perc=(volume_l_s/volume_l)*100


%% 5. Add noise to the original dataset and check the performances of your implemented workflow with respect to different levels of noise.
rand_IM = rand(256,256);
figure, 
imshow(rand_IM)
title('additional noise')

rum = double(vol(:,:,:))./255+rand_IM;
rum_sc = rescale(rum(:,:,:),0,1);

figure, 
subplot(1,2,1)
imshow(vol(:,:,75)) % uint8 da 0 a 255
title('original image')
subplot(1,2,2)
imshow(rum_sc(:,:,75)) % double da 0 a 1
title('image with additional noise')

%%
[Cropped_vol d]= imcrop(rum_sc(:,:,75), [130 102 51 45]);
v1=round(d(2)):(round(d(2))+length(Cropped_vol(:,1)));
v2=round(d(1)):(round(d(1))+length(Cropped_vol(1,:)));
v3=64:90;
VOI=rum_sc(v1,v2,v3);

%%
% MASK=VOI(:,:,25:27);
% MASK(MASK>(240/255))=0;
% VOI(:,:,25:27)=MASK;

clear vol_imadjusted

for i=1:size(VOI,3)
        vol_imadjusted(:,:,i) = imadjust(VOI(:,:,i),[0 0.5882],[0 1],2);
end 
figure()
montage(vol_imadjusted)

clear vol_pn
for i=1:length(v3)
    vol_pn(:,:,i)=medfilt2(vol_imadjusted(:,:,i), [6 6]);
end

figure()
montage(vol_pn)

%%
bin_vol=imbinarize(vol_pn,0.8);

figure()

for i=2:25
    
    imshow(VOI(:,:,i))
    title("Contour of the tumor")
    hold on
    imcontour(bin_vol(:,:,i),5,'m');
    pause (1)
end

%%
% contorni troppo irregolari --> aggiungiamo un low-pass filter e
% riapplichiamo la procedura
img_filt=zeros(size(rum_sc,1),size(rum_sc,2),size(rum_sc,3));
k=[1/9 1/9 1/9;
   1/9 1/9 1/9;
   1/9 1/9 1/9];
for i=1:size(rum_sc,3)
    img_filt(:,:,i) = imfilter(rum_sc(:,:,i),k,'conv','symmetric');
end

figure, 
subplot(1,2,1)
imshow(rum_sc(:,:,75)) % uint8 da 0 a 255
title('original image')
subplot(1,2,2)
imshow(img_filt(:,:,75)) % double da 0 a 1
title('image with additional noise')

VOI_filt=img_filt(v1,v2,v3);

MASK=VOI_filt(:,:,25:27);
MASK(MASK>(240/255))=0;
VOI_filt(:,:,25:27)=MASK;


clear vol_imadjusted
for i=1:size(VOI_filt,3)
        vol_imadjusted(:,:,i) = imadjust(VOI_filt(:,:,i),[0 0.5882],[0 1],2);
end 

clear vol_pn
for i=1:length(v3)
    vol_pn(:,:,i)=medfilt2(vol_imadjusted(:,:,i), [6 6]);
end

figure()
montage(vol_pn)

bin_vol=imbinarize(vol_pn,0.8);
figure()
montage(bin_vol)

figure()
for i=2:25
    
    imshow(VOI(:,:,i))
    title("Contour of the tumor")
    hold on
    imcontour(bin_vol(:,:,i),5,'m');
    pause (1)
end

num_pixel_noise=[];
for i=2:25
    num_pixel_noise=[num_pixel_noise sum(sum(bin_vol(:,:,i)==1))]; %conta i pixel bianchi 
end 

% total volume of the lesion in mm^3
volume_l_noise=0;
for i=1:24
     volume_l_noise=volume_l_noise+num_pixel_noise(1,i).*pixdim(1,3).*pixdim(1,1).*pixdim(1,2);
end

perc=(volume_l_noise/volume_l)*100

