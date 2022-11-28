clc
clear
close all 
%% Animation on: a=1 Animation off: a=0;
a=0;
%% Visualizzazione dei dati 
load MRIdata.mat
cmap = colormap('jet');

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


for i = 1:(90-64)
    img(:,:,i) = imcrop(vol(:,:,64+i),[140 85 55 55]);
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


%% 2. Implement a workflow, motivating its main sub-steps, to pre-process an MRI image and segment a lesion based on methods we have seen during lessons.
%SEGMENTATION:
%The goal of image segmentation is to divide an image into a set of 
% semantically meaningful,homogeneous,nonoverlapping 
% regions of similar attributes such as intensity, depth, color, or
% texture. The segmentation result is either an image of labels identifying each 
% homogeneous region or a set ofcontours which describe the region boundaries.


%% 3. Segment the lesion and calculate the respective cross-sectional area over sagittal slice number 135 4. Identify sagittal slices that contain the lesion and extend the quantification of its cross-sectional area to the whole volume. Try to repeat this process across axial slices. What are the main challenges of segmenting this lesion with respect to other cerebral tissues and orthogonal views?

%% 5. Add noise to the original dataset and check the performances of your implemented workflow with respect to different levels of noise.

%% 6. [Optional] manually segment the lesion starting from sagittal slice number 135, hence quantify segmentation performances in terms of sensitivity, specificity and Dice coefficient.
