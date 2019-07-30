clear all;
close all;
clc;
% Read an original spec and its version for segmentation.
spec_file = 'wind_sea_NODE008919_201001.png';
spec_file2 = 'wind_sea2_NODE008919_201001.png';

directory = sprintf('/homes/v18nguye/Documents/intern2019/sparse-waves/sparse-waves/data/WW3/%s',spec_file);
directory2 = sprintf('/homes/v18nguye/Documents/intern2019/sparse-waves/sparse-waves/data/WW3/%s',spec_file2);
%directory2 = sprintf('homes\\v18nguye\\Documents\\intern2019\\sparse-waves\\sparse-waves\\data\\WW3\\%s',spec_file2);
%directory2 = sprintf('E:\\IMT\\2019 Summer Internship at IMT Atlantique\\sparse-waves\\data\\WW3\\%s',spec_file2);

spec = imread(directory);
spec2 = imread(directory2);
% Select a zone of interest of the image for segmentation
figure;
imshow(spec2);
h_rect = imrect();
%%% get rectangle position
pos_rect = h_rect.getPosition();
%%% round off the coordinates can be used as indices
pos_rect = round(pos_rect);
%%% selected region of the image
spec_cropped = spec(pos_rect(2)+(0:pos_rect(4)),pos_rect(1)+(0:pos_rect(3)),:);
spec_cropped2 = spec2(pos_rect(2)+(0:pos_rect(4)),pos_rect(1)+(0:pos_rect(3)),:);


% Covert the cropped image into grayscale
I = rgb2gray(spec_cropped2);

% Compute the gradient magnitude of this type of the image
%gmag = imgradient(I,'intermediate');
gmag = imgradient(I,'central');

% Mark the foreground objects (morphological technique called "opening-by-reconstruction
% and closing-by-reconstruction")
%%% structuring element
se = strel('disk',10);
%%% compute the opening-by-reconstruction
Ie = imerode(I,se); % marker
Iobr = imreconstruct(Ie,I); % marker, mask
figure
imshow(Iobr)
%%% compute closing-by-reconstruction
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
figure
imshow(Iobrcbr)

% Obtain the good foreground markers
fgm = imregionalmax(Iobrcbr);
se2 = strel(ones(5,5));
fgm2 = imclose(fgm,se2);
fgm3 = imerode(fgm2,se2);

% Exclude stray isolated pixels
fgm4 = bwareaopen(fgm3,3);

%Compute background markers
%%% start thresholding the image previously obtained
bw = imbinarize(Iobrcbr);
figure
imshow(Iobrcbr)
figure
imshow(bw)
%%% Compute the watershed transform of the distance transform to
%%% obtain the watershed ridge lines from the background markers.
D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;
figure
imshow(bgm)
% Compute the watershed Transform of the segmentatioin function
gmag2 = imimposemin(gmag, bgm | fgm4);
figure
imshowpair(gmag, gmag2, 'montage')
L = watershed(gmag2);
figure 
imshow(L)

% Visalilize the result
labels = imdilate(L==0,ones(1,1))+ 3*fgm4;

I2 = labeloverlay(spec_cropped,labels);
spec(pos_rect(2)+(0:pos_rect(4)),pos_rect(1)+(0:pos_rect(3)),:) = I2;
figure
imshow(spec)




