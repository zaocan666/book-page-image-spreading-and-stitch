function [f, d]=sift(img)
if length(size(img))==3
    img = rgb2gray(img);
end
img = single(img); %灰度值归一化
[f,d] = vl_sift(img) ;
max_num = 22000;
sub = vl_colsubset(1 : size(f,2), max_num) ;
f = f(:, sub);
d = d(:, sub);
end