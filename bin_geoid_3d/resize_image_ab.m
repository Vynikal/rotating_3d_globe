function img_new=resize_image_ab(img,zoom_factor)
%IMG_NEW=RESIZE_IMAGE_AB(IMG,ZOOM_FACTOR)
% Resizing of images using Gaussian blurring given by 'zoom_factor'
% 'img', 'img_new' are input and output images

% Uses: conv2padded by Pascal Getreuer, http://www.getreuer.info/tutorials/matlabimaging
% Ales Bezdek, bezdek@asu.cas.cz, 1/2015

%% Gaussian filtering
% I found the value of sigma empirically
sigma=1/zoom_factor/3.5;
radius=ceil(3*sigma);    % limit the Gaussian kernel

h1=exp(-(-radius:radius).^2/(2*sigma^2));
h1=h1/sum(h1);                 % normalization
%       img_new=conv2padded(h1,h1,img);

oldClass=class(img);  % original image type
img=double(img);      % convert image to double precision for processing
img_size=size(img);

% if numel(img_size)==2 %greyscale
%    img=conv2(h1,h1,img,'same');
% else
%    img(:,:,1)=conv2(h1,h1,img(:,:,1),'same');
%    img(:,:,2)=conv2(h1,h1,img(:,:,2),'same');
%    img(:,:,3)=conv2(h1,h1,img(:,:,3),'same');
% end
img = conv2padded(h1,h1,img);
%% resizing of slightly blurred image using the bicubic interpolation
new_size=round(img_size(1:2)*zoom_factor);
x1=((1:new_size(2))-0.5)./zoom_factor+0.5;  % new image pixel X coordinates
y1=((1:new_size(1))-0.5)./zoom_factor+0.5;  % new image pixel Y coordinates
if numel(img_size)==2 %greyscale
   img_new=interp2(img,x1,y1(:),'cubic');
else
   img_new=zeros([new_size 3]);
   for i=1:3
      img_new(:,:,i)=interp2(img(:,:,i),x1,y1(:),'cubic');
   end
end

img_new=cast(img_new,oldClass);  % Convert back to original image type
end
