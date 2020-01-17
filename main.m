clear;
auto_flag = false; %是否自动取点

raw_image_1=imread('img1.jpg');
raw_image_2=imread('img2.jpg');
pad_radio = 0.1;

%将书页图片展平
p_dis=0;
if auto_flag==false
    p_dis = 0;
    new_img = pad_resize(raw_image_1, pad_radio);
    load('select_ps_1.mat');
    undistort_img_1 = undistort(select_ps_1, new_img, size(select_ps_1,1)/10, 10, 1, pad_radio);
    close all;
    %imshow(uint8(undistort_img_1));
    %imwrite(uint8(undistort_img_1), 'undistort_1.jpg');

    new_img = pad_resize(raw_image_2, pad_radio);
    load('select_ps_2.mat');
    undistort_img_2 = undistort(select_ps_2, new_img, size(select_ps_2,1)/8, 8, 2, pad_radio);
    close all;
    %imshow(uint8(undistort_img_2));
    %imwrite(uint8(undistort_img_2), 'undistort_2.jpg');
else
    p_dis = 20;
    undistort_img_1 = auto_point(raw_image_1, pad_radio, 7);
    undistort_img_2 = auto_point(raw_image_2, pad_radio, 9);
    %imwrite(uint8(undistort_img_1), 'undistort_1.jpg');
    %imwrite(uint8(undistort_img_2), 'undistort_2.jpg');
end

%两张图片进行拼接
% undistort_img_1 = imread('undistort_1.jpg');
% undistort_img_2 = imread('undistort_2.jpg');
stitch_img = stitch(undistort_img_1, undistort_img_2, p_dis);
[x_min, x_max, y_min, y_max] = get_bound(stitch_img);
bound_img = stitch_img(x_min:x_max, y_min:y_max);
bound_img = uint8(round(bound_img));
bound_img(bound_img<30)=200;

subplot(1,2,1)
imshow(bound_img);
imwrite(bound_img, 'result.jpg');

%二值化
% bound_img = imread('result.jpg');
img_bw = myBinarize(bound_img);
subplot(1,2,2)
imshow(img_bw);
imwrite(img_bw, 'binarize_result.jpg');

function [x_min, x_max, y_min, y_max] = get_bound(img)
    logic_img = (img ~= 0); %有图像的地方为1
    
    x_min=1;
    x_max=size(logic_img, 1);
    for j=1:size(logic_img,2)
        for i=1:round(size(logic_img,1)/4)
            if logic_img(i,j)==1
                if i>x_min
                    x_min=i;
                end
                break;
            end
        end
        
        for i=size(logic_img,1):-1:round(size(logic_img,1)*3/4)
            if logic_img(i,j)==1
                if i<x_max
                    x_max=i;
                end
                break;
            end
        end
    end
    
    y_min=1;
    y_max=size(logic_img, 2);
    for i=x_min:x_max
        for j=1:round(size(logic_img,2)/4)
            if logic_img(i,j)==1
                if j>y_min
                    y_min=j;
                end
                break;
            end
        end
        
        for j=size(logic_img,2):-1:round(size(logic_img,2)*3/4)
            if logic_img(i,j)==1
                if j<y_max
                    y_max=j;
                end
                break;
            end
        end
    end
    
end