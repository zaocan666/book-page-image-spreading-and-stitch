function mix_img = multiFre_mix(img1, img2, k_mask1)
%来源：https://blog.csdn.net/ccblogger/article/details/70665552
%经修改
%高斯滤波
kernel=fspecial('gaussian',[5 5],1);

img1(isnan(img1))=0;
img2(isnan(img2))=0;
if max(max(img2))<=1
    img2=img2*255;
end

%高斯金字塔
G_A0 = img1;
G_A1 = conv2(G_A0,kernel,'same');
G_A1 = G_A1(2:2:size(G_A1,1),2:2:size(G_A1,2));
G_A2 = conv2(G_A1,kernel,'same');
G_A2 = G_A2(2:2:size(G_A2,1),2:2:size(G_A2,2));
G_A3 = conv2(G_A2,kernel,'same');
G_A3 = G_A3(2:2:size(G_A3,1),2:2:size(G_A3,2));
G_A4 = conv2(G_A3,kernel,'same');
G_A4 = G_A4(2:2:size(G_A4,1),2:2:size(G_A4,2));
G_A5 = conv2(G_A4,kernel,'same');
G_A5 = G_A5(2:2:size(G_A5,1),2:2:size(G_A5,2));
 
G_B0 = img2;
G_B1 = conv2(G_B0,kernel,'same');
G_B1 = G_B1(2:2:size(G_B1,1),2:2:size(G_B1,2));
G_B2 = conv2(G_B1,kernel,'same');
G_B2 = G_B2(2:2:size(G_B2,1),2:2:size(G_B2,2));
G_B3 = conv2(G_B2,kernel,'same');
G_B3 = G_B3(2:2:size(G_B3,1),2:2:size(G_B3,2));
G_B4 = conv2(G_B3,kernel,'same');
G_B4 = G_B4(2:2:size(G_B4,1),2:2:size(G_B4,2));
G_B5 = conv2(G_B4,kernel,'same');
G_B5 = G_B5(2:2:size(G_B5,1),2:2:size(G_B5,2));
 
%拉普拉斯金字塔
L_A0 = double(G_A0)-imresize(G_A1,size(G_A0));
L_A1 = double(G_A1)-imresize(G_A2,size(G_A1));
L_A2 = double(G_A2)-imresize(G_A3,size(G_A2));
L_A3 = double(G_A3)-imresize(G_A4,size(G_A3));
L_A4 = double(G_A4)-imresize(G_A5,size(G_A4));
L_A5 = double(G_A5);
 
L_B0 = double(G_B0)-imresize(G_B1,size(G_B0));
L_B1 = double(G_B1)-imresize(G_B2,size(G_B1));
L_B2 = double(G_B2)-imresize(G_B3,size(G_B2));
L_B3 = double(G_B3)-imresize(G_B4,size(G_B3));
L_B4 = double(G_B4)-imresize(G_B5,size(G_B4));
L_B5 = double(G_B5);
 
mask0 = k_mask1;
mask1 = imresize(mask0,size(L_A1));
mask2 = imresize(mask0,size(L_A2));
mask3 = imresize(mask0,size(L_A3));
mask4 = imresize(mask0,size(L_A4));
mask5 = imresize(mask0,size(L_A5));
 
%获得结果
L_C0 = L_A0 .* mask0 + L_B0 .* (1-mask0);
L_C1 = L_A1 .* mask1 + L_B1 .* (1-mask1);
L_C2 = L_A2 .* mask2 + L_B2 .* (1-mask2);
L_C3 = L_A3 .* mask3 + L_B3 .* (1-mask3);
L_C4 = L_A4 .* mask4 + L_B4 .* (1-mask4);
L_C5 = L_A5 .* mask5 + L_B5 .* (1-mask5);
size0=size(L_C0);
mix_img = L_C0+imresize(L_C1,size0)+imresize(L_C2,size0)+imresize(L_C3,size0)+imresize(L_C4,size0)+imresize(L_C5,size0);
 
% figure(1);
% imshow(img1);
% figure(2);
% imshow(img2);
% figure(3);
% imshow(uint8(mix_img));
 
end