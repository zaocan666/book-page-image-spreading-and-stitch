function img=pad_resize(raw_img, pad_radio)
small_raw_image = imresize(raw_img, 0.4);
img = rgb2gray(small_raw_image);
img = padarray(img, [round(size(img, 1)*pad_radio) round(size(img, 2)*pad_radio)], 'both');
end