close all;
clear all;
pic_name = 'img1.jpg';
point_name = 'select_ps_1.mat';

if strcmp(point_name,'')==0
    load(point_name);
end

if exist('select_ps_1')
    select_ps = select_ps_1;
else
    select_ps = select_ps_2;
end

raw_image=imread(pic_name);
new_img=pad_resize(raw_image,0.1);

imshow(new_img,'border','tight','initialmagnification','fit');
%set (gcf,'Position',[0,0,size(new_img,2),size(new_img,1)]);
axis normal;
hold on;
plot( select_ps(:,1),select_ps(:,2),'r.' );
frame = getframe(gca);
plot_img = frame2im(frame);
plot_img=imresize(plot_img, [size(new_img,1) size(new_img,2)]);
cpselect(plot_img,plot_img);

new_select_ps=[select_ps; new_ps];
close all;
imshow(new_img,'border','tight','initialmagnification','fit');
hold on;
plot( new_select_ps(:,1),new_select_ps(:,2),'r.' );

if exist('select_ps_1')
    select_ps_1 = new_select_ps;
else
    select_ps_2 = new_select_ps;
end

