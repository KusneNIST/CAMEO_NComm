function gc_example()
% An example of how to segment a color image according to pixel colors.
% Fisrt stage identifies k distinct clusters in the color space of the
% image. Then the image is segmented according to these regions; each pixel
% is assigned to its cluster and the GraphCut poses smoothness constraint
% on this labeling.
close all

% read an image
%im = im2double(imread('outdoor_small.jpg'));

 im = imread('C:\Users\akusne\Dropbox\Combi_work\Materials Informatics\Written Sections\3-phase_material.png');
im = double(im) + .01*rand(size(im));
im = im / max(im(:));
im = im(100:130,30:100,:);
sz = size(im);
%im(1:5,1:10,:) = 40;
figure(1); clf; image(im); colorbar;
% figure(1); imagesc(im);
sz = size(im);

% try to segment the image into k different regions
k = 6;

% color space distance
distance = 'sqEuclidean';

try_number = 0;
kmeans_success = 0;

% cluster the image colors into k regions
data = double(ToVector(im));

while(kmeans_success == 0 && try_number < 100)
    try
        try_number = try_number + 1;
        [idx c] = kmeans(data, k); %, 'distance', distance,'maxiter',500);
        kmeans_success = 1;
    catch err
        err.message
    end
end
% calculate the data cost per cluster center
Dc = zeros([sz(1:2) k],'single');
for ci=1:k
    % use covariance matrix per cluster
    icv = inv(cov(data(idx==ci,:)));    
    dif = data - repmat(c(ci,:), [size(data,1) 1]);
    % data cost is minus log likelihood of the pixel to belong to each
    % cluster according to its RGB value
    Dc(:,:,ci) = reshape(sum((dif*icv).*dif./2,2),sz(1:2));
end

% cut the graph

% smoothness term: 
% constant part
Sc = ones(k) - eye(k);
% spatialy varying part
% [Hc Vc] = gradient(imfilter(rgb2gray(im),fspecial('gauss',[3 3]),'symmetric'));
[Hc Vc] = SpatialCues(im);

gch = GraphCut('open', Dc, 10*Sc, exp(-Vc*5), exp(-Hc*5));
[gch L] = GraphCut('expand',gch);
gch = GraphCut('close', gch);

% show results
imshow(im);
hold on;
PlotLabels(L);
hold off;
figure(2); imagesc(reshape(L,sz(1),sz(2)));


%---------------- Aux Functions ----------------%
function v = ToVector(im)
% takes MxNx3 picture and returns (MN)x3 vector
sz = size(im);
v = reshape(im, [prod(sz(1:2)) 3]);

%-----------------------------------------------%
function ih = PlotLabels(L)

L = single(L);

bL = imdilate( abs( imfilter(L, fspecial('log'), 'symmetric') ) > 0.1, strel('disk', 1));
LL = zeros(size(L),class(L));
LL(bL) = L(bL);
Am = zeros(size(L));
Am(bL) = .5;
ih = imagesc(LL); 
set(ih, 'AlphaData', Am);
colorbar;
colormap 'jet';

%-----------------------------------------------%
function [hC vC] = SpatialCues(im)
g = fspecial('gauss', [13 13], sqrt(13));
dy = fspecial('sobel');
vf = conv2(g, dy, 'valid');
sz = size(im);

vC = zeros(sz(1:2));
hC = vC;

for b=1:size(im,3)
    vC = max(vC, abs(imfilter(im(:,:,b), vf, 'symmetric')));
    hC = max(hC, abs(imfilter(im(:,:,b), vf', 'symmetric')));
end