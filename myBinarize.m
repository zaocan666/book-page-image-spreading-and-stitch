function img_bw = myBinarize(bound_img)
balanced_cross_img = bound_img;

%横向均衡
part_num = 8;
img_1 = bound_img(:, 1:round(size(bound_img,2)/part_num));
average_1 = mean(mean(img_1));

for i =1:part_num-1
    left = size(bound_img,2)/part_num*i;
    right = size(bound_img,2)/part_num*(i+1);
    img_i = bound_img(:, round(left):round(right));
    average_i = mean(mean(img_i));
    img_i = double(img_i);
    img_i = img_i * average_1/average_i;
    balanced_cross_img(:, round(left):round(right)) = uint8(img_i);
end

balanced_cross_img = uint8(balanced_cross_img);
% imshow(balanced_cross_img);
% img_bw = imbinarize(balanced_cross_img,160/255);
% imshow(img_bw)

balanced_length_img = balanced_cross_img;
%纵向均衡
part_num = 6;
img_1 = balanced_cross_img(1:round(size(bound_img,1)/part_num), :);
average_1 = mean(mean(img_1));

for i =1:part_num-1
    up = size(bound_img,1)/part_num*i;
    down = size(bound_img,1)/part_num*(i+1);
    img_i = balanced_cross_img(round(up):round(down), :);
    average_i = mean(mean(img_i));
    img_i = double(img_i);
    img_i = img_i * average_1/average_i;
    balanced_length_img(round(up):round(down), :) = uint8(img_i);
end

balanced_length_img = uint8(balanced_length_img);
% imshow(balanced_length_img);
img_bw = imbinarize(balanced_length_img,160/255);
end