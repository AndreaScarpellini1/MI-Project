clc
clear
close all 
%cd 'C:\Users\scrpa\OneDrive - Politecnico di Milano\Desktop\Poli\Magistrale\Primo anno\BSPMI\MI\project repo\MI-Project'
%% Animation on: a=1 Animation off: a=0;
a=1;
%% Visualizzazione dei dati 
load MRIdata.mat

if (a==1)
    figure()
    for i=1:length(vol(1,1,:))
        subplot(2,1,1)
           imshow(vol(:,:,i))
           colorbar
           title("MRI iniziale.")
        subplot(2,1,2)
           imshow(vol(:,:,i),'Colormap', cmap)
           colorbar
           title("MRI colored")
           pause(0.01)
    end 
end

%% crop
figure()
imshow(vol(:,:,90))

clear img
for i = 1:(90-64)
    img(:,:,i) = imcrop(vol(:,:,64+i),[140 100 40 40]);
end
montage(img)

% da vedere le misure

%% Binarize image
for i=1:length(vol(1,1,:))
    bin_vol(:,:,i) = imbinarize(vol(:,:,i));
end
%%
figure()
for i=1:length(vol(1,1,:))
    imshow(bin_vol(:,:,i))
    pause(0.01)
end
%%
figure()
imshow(fant_pet_bin)
% Binarize the image using as threshold half of the maximum
%%
imcontour(fant_pet_bin, 1, 'm');

%%
figure()
for i=1:length(vol(1,1,:))
       imshow(vol_cont(:,:,i))
       title("MRI binarized.")
       pause(0.01)
end

%% 1. Briefly review the topic of tissue and lesion segmentation over MRI images

if (a==1)
    figure(1)
        montage(vol)
        title("MRI iniziale.")
    figure(2)
        montage(vol(:,:,64:90))
        title('Tumor From slice 64 to 90.')
end 



%% SEGMENTATION:
%The goal of image segmentation is to divide an image into a set of 
% semantically meaningful,homogeneous,nonoverlapping 
% regions of similar attributes such as intensity, depth, color, or
% texture. The segmentation result is either an image of labels identifying each 
% homogeneous region or a set ofcontours which describe the region boundaries.

figure()
subplot(2,1,1)
imshow(vol(:,:,75))
impixelinfo
title("Isolate the tumor")
subplot(2,1,2)
[Cropped_vol d]= imcrop(vol(:,:,75));

% Dimensioni del taglio 
v1=round(d(2)):(round(d(2))+length(Cropped_vol(:,1)));
v2=round(d(1)):(round(d(1))+length(Cropped_vol(1,:)));
v3=64:90;

figure()
subplot(2,1,1)
imshow(vol(v1,v2,75));
title('Immagine ricavata dalle dimensioni')
subplot(2,1,2)
title('Immagine ricavata dalla funzione')
imshow(Cropped_vol)

%Aumento del contrasto
j=0;
for i=v3
    j=j+1;
    vol_imadjust(:,:,j) = imadjust(vol(v1,v2,i));
end 

%salt & pepper filtering
for i=1:length(v3)
    vol_imadjust(:,:,i)=medfilt2(vol_imadjust(:,:,i), [5 5]);
end

%Binarizzazione 
bin_vol=imbinarize(vol_imadjust,0.4);

figure()
subplot(1,2,1)
montage(vol_imadjust)
title('Enhanced contrast')
subplot(1,2,2)
montage(bin_vol)

%Prendo i contorni 
figure()
for i=2:26
    imshow(bin_vol(:,:,i))
    hold on
    imcontour(bin_vol(:,:,i),3,'m')
    pause (1)
end 
title("Contours of the tumor")
%%
%3D Visualization 
volumeViewer(vol(v1,v2,v3))
%% 3. Segment the lesion and calculate the respective cross-sectional area over sagittal slice number 135 4. Identify sagittal slices that contain the lesion and extend the quantification of its cross-sectional area to the whole volume. Try to repeat this process across axial slices. What are the main challenges of segmenting this lesion with respect to other cerebral tissues and orthogonal views?

%% 5. Add noise to the original dataset and check the performances of your implemented workflow with respect to different levels of noise.

%% 6. [Optional] manually segment the lesion starting from sagittal slice number 135, hence quantify segmentation performances in terms of sensitivity, specificity and Dice coefficient.
