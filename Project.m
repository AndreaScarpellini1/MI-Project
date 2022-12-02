clc
clear
close all 
cd 'C:\Users\scrpa\OneDrive - Politecnico di Milano\Desktop\Poli\Magistrale\Primo anno\BSPMI\MI\project repo\MI-Project'
%% Animation on: a=1 Animation off: a=0;
a=1;
%% dati 
load MRIdata.mat
%% Segmentation and AMOUNT AXIAL PIXEL
if (a==0)
    figure(1)
        montage(vol)
        title("MRI iniziale.")
    figure(2)
        montage(vol(:,:,64:90))
        title('Tumor From slice 64 to 90.')
end 

figure()
subplot(2,1,1)
imshow(vol(:,:,75))
impixelinfo
title("Isolate the tumor")
subplot(2,1,2)
[Cropped_vol d]= imcrop(vol(:,:,75), [130 102 51 45]);

% Dimensioni del taglio 
v1=round(d(2)):(round(d(2))+length(Cropped_vol(:,1)));
v2=round(d(1)):(round(d(1))+length(Cropped_vol(1,:)));
v3=64:90;

figure('Name',"Confronto funzione matlab e dimensioni ottenute")
subplot(2,1,1)
imshow(vol(v1,v2,75));
title('Immagine ricavata dalle dimensioni')
subplot(2,1,2)
imshow(Cropped_vol)
title('Immagine ricavata dalla funzione')

%% studio histogrammi 
VOI=vol(v1,v2,v3);
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
% gamma --> 2 


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

%% Prendo i contorni 

vidfile = VideoWriter('testmovie.mp4','MPEG-4');
open(vidfile);
figure()
j=0;
for i=2:26
    j=j+1;
    imshow(VOI(:,:,i))
    title("Contour of the tumor")
    hold on
    imcontour(bin_vol(:,:,i),5,'m');
    F(j) = getframe(gcf); 
    writeVideo(vidfile,F(j));
    pause (1)
end 
close(vidfile)


%%
% area of the binarized image 
Axial_num_pixel=0;
for i=2:26
    Axial_num_pixel=Axial_num_pixel+sum(sum(bin_vol(:,:,i)==1)); %conta i pixel bianchi 
end 

%%
%3D Visualization 
volumeViewer(vol(v1,v2,v3))
%% 3. Segment the lesion and calculate the respective cross-sectional area over sagittal slice number 135 4. Identify sagittal slices that contain the lesion and extend the quantification of its cross-sectional area to the whole volume. Try to repeat this process across axial slices. What are the main challenges of segmenting this lesion with respect to other cerebral tissues and orthogonal views?
%from axial to sagittal plane:
for i=1:dim(1)
    vol_ax(:,:,i)=vol(i,:,:);
end

clear vol_imadjust
%display images
if (a==0)
    figure(1)
        montage(vol_ax)
        title("MRI Sagittale.")
    figure(2)
        montage(vol_ax(:,:,107:144))
        title('Tumor From slice 107 to 144.')
end 
[Cropped_vol_ax d_ax]= imcrop(vol_ax(:,:,126), [55 135 35 60]);

% Dimensioni del taglio 
v1=round(d_ax(2)):(round(d_ax(2))+length(Cropped_vol_ax(:,1)));
v2=round(d_ax(1)):(round(d_ax(1))+length(Cropped_vol_ax(1,:)));
v3=107:144;

figure()
subplot(2,1,1)
imshow(vol_ax(v1,v2,126));
title('Immagine ricavata dalle dimensioni')
subplot(2,1,2)
title('Immagine ricavata dalla funzione')
imshow(Cropped_vol_ax)

%Aumento del contrasto
j=0;
for i=v3
    j=j+1;
    vol_imadjust(:,:,j) = imadjust(vol_ax(v1,v2,i));
end 
%salt & pepper filtering
for i=1:length(v3)
    vol_imadjust(:,:,i)=medfilt2(vol_imadjust(:,:,i), [5 5]);
end

%Binarizzazione 
bin_vol=imbinarize(vol_imadjust,0.6);

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
    imcontour(bin_vol(:,:,i),4,'m')
    pause (1)
end 
title("Contours of the tumor")
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

%% 6. [Optional] manually segment the lesion starting from sagittal slice number 135, hence quantify segmentation performances in terms of sensitivity, specificity and Dice coefficient.
