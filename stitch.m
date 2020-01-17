function stitch_img=stitch(img1, img2, p_dis)
% 来自MATLAB官网vlfeat包
% 经修改
run('vlfeat-0.9.21-bin\vlfeat-0.9.21\toolbox\vl_setup.m');
img1 = imresize(img1, 0.6);
img2 = imresize(img2, 0.6);

if size(img1,3) > 1
    img1_g = rgb2gray(img1) ; 
else
    img1_g = img1 ;
end
if size(img2,3) > 1
    img2_g = rgb2gray(img2); 
else
    img2_g = img2 ; 
end

img1_g = single(img1_g);
img2_g = single(img2_g);
[f1,d1] = sift(img1_g) ;%图像特征点提取 f是坐标 d是对象描述
[f2,d2] = sift(img2_g) ;

[H, ok, P1, P2, matches] = match_features(d1, d2, f1, f2);

matches_p1 = P1(:, ok);
matches_p2 = P2(:, ok);
[max1y, index_max1y] = max(matches_p1(2,:));
[min1y, index_min1y] = min(matches_p1(2,:));

radio=(max1y-min1y)/(matches_p2(2, index_max1y)-matches_p2(2, index_min1y));
img2=imresize(img2, radio); %使im1与im2高度比例合适，然后重新匹配
if size(img2,3) > 1
    img2_g = rgb2gray(img2); 
else
    img2_g = img2 ; 
end
img2_g = single(img2_g);
[f2,d2] = sift(img2_g) ;
[H, ok, P1, P2, matches] = match_features(d1, d2, f1, f2);

P2_ok = f2(1:2,matches(2,ok)); 
P2_ok(3,:)=1;
P1_2ok = inv(H)*P2_ok;

%再次TPS使得两张图公共部分重合
%p_subset = vl_colsubset(1:size(P1_2ok, 2), 35);
num_ps = 60;%15;
source_ps = f1(1:2,matches(1,ok));
aim_ps = P1_2ok(1:2,:);

% imshow(img1_g);hold on;
% plot( source_ps(1,:),source_ps(2,:),'r.' );
% plot( aim_ps(1,:),aim_ps(2,:),'g.' );
[source_ps, aim_ps] = eliminate_tooClose(source_ps, aim_ps, p_dis);
[~, s_index]= sort(source_ps(2, :));
source_ps = source_ps(:, s_index);
aim_ps = aim_ps(:, s_index);
index_delta = round(size(source_ps, 2)/num_ps);
if index_delta<1
    index_delta=1;
end
num_index = 1:index_delta:size(source_ps, 2);

source_ps = source_ps(:, num_index);
aim_ps = aim_ps(:, num_index);
% source_ps = source_ps(:, p_subset);
% aim_ps = P1_2ok(1:2,:);
% aim_ps = aim_ps(:, p_subset);
[img1,mask1] = rbfwarp2d( img1, (source_ps)', (aim_ps)','thin');
%imagesc(img1) ; axis image off ;
if size(img1,3) > 1
    img1_g = rgb2gray(img1); 
else
    img1_g = img1 ; 
end
img1_g = single(img1_g);
[f1,d1] = sift(img1_g) ;
[H, ok, P1, P2, matches] = match_features(d1, d2, f1, f2);

% dx1 = max(size(img2,1)-size(img1,1),0) ;
% dx2 = max(size(img1,1)-size(img2,1),0) ;
%匹配的所有点
% figure(1); 
% clf;
% imagesc([padarray(img1,dx1,'post') padarray(img2,dx2,'post')]) ;
% pad = size(img1,2) ;
% line([f1(1,matches(1,:));f2(1,matches(2,:))+pad], ...
%      [f1(2,matches(1,:));f2(2,matches(2,:))]) ;
% axis image off ;

%在变换H下仍然匹配的点
% figure(2); 
% clf;
% imagesc([padarray(img1,dx1,'post') padarray(img2,dx2,'post')]) ;
% pad = size(img1,2) ;
% line([f1(1,matches(1,ok));f2(1,matches(2,ok))+pad], ...
%      [f1(2,matches(1,ok));f2(2,matches(2,ok))]) ;
% axis image off ;
%  
% print(gcf,'-djpeg','matchLine.jpeg');
% drawnow ;
 
% --------------------------------------------------------------------
%                                                               Mosaic
% --------------------------------------------------------------------
 
bound2 = [1  size(img2,2) size(img2,2)  1 ;
        1  1           size(img2,1)  size(img2,1) ;
        1  1           1            1 ] ;%im2四个边界点
bound2_ = inv(H) * bound2 ;
bound2_(1,:) = bound2_(1,:) ./ bound2_(3,:) ;
bound2_(2,:) = bound2_(2,:) ./ bound2_(3,:) ;
ur = min([1 bound2_(1,:)]):max([size(img1,2) bound2_(1,:)]) ;
vr = min([1 bound2_(2,:)]):max([size(img1,1) bound2_(2,:)]) ;
 
[u,v] = meshgrid(ur,vr) ;
img1_im = vl_imwbackward(im2double(img1),u,v) ;
 
z_ = H(3,1) * u + H(3,2) * v + H(3,3) ;
u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;%u_, v_为把p1投影到p2

img2_im = vl_imwbackward(im2double(img2),u_,v_) ; %在特征点p1取特征点p2的灰度值

img1_im(img1_im==0)=NaN;
img2_im(img2_im==0)=NaN;
img_12_mask=~isnan(img1_im);
img_12_mask(~isnan(img2_im))=1; %有图像的部分
[k_mask1]=get_stitch_mask(img1_im, img2_im);
mix_img = multiFre_mix(img1_im, img2_im, k_mask1); %多频段融合
mix_img(~img_12_mask)=0; %无图像部分为0
% imagesc(mix_img) ; axis image off ;
stitch_img=mix_img;

% common = ((~isnan(img1_im) + ~isnan(img2_im))==2) ; %二者共同部分为1
% k=~isnan(img1_im) + ~isnan(img2_im);%二者共同部分为2，私有部分为1，其他部分为0
% 
% img1_im_ave = img1_im;
% img1_im_ave(~common)=0;
% img2_im_ave = img2_im;
% img2_im_ave(~common)=0;
% img1_com_average = sum(sum(img1_im_ave));
% img2_com_average = sum(sum(img2_im_ave));
% img2_im = img2_im*img1_com_average/img2_com_average; %根据共同部分的亮度平均值，使两张图片亮度均衡
% 
% img1_im(isnan(img1_im)) = 0 ;
% img2_im(isnan(img2_im)) = 0 ;
% img1_com = img1_im;
% img1_com(common)=0;
% stitch_img = (img1_im + img2_im)./k;
% stitch_img(isnan(stitch_img))=0;
%  
% figure(3) ; clf ;
% imagesc(stitch_img) ; axis image off ;

end

function [H, ok, P1, P2, matches] = match_features(d1, d2, f1, f2)
    [matches, scores] = vl_ubcmatch(d1,d2) ; %特征点匹配

    num_match = size(matches,2) ;

    P1 = f1(1:2,matches(1,:)) ; %相互匹配的点
    P1(3,:) = 1 ;
    P2 = f2(1:2,matches(2,:)) ; 
    P2(3,:) = 1 ;
    
    clear H score ok ;
%     for t = 1:100   % 迭代100次，选择最佳的H
%       subset = vl_colsubset(1:num_match, 4) ;%随机取四个数
%       A = [] ;
%       for i = subset
%         A = cat(1, A, kron(P1(:,i)', vl_hat(P2(:,i)))) ;%cat(1,A,B)在列方向上拼接矩阵，kron(A,B)表示[A(1,1)*B,A(1,2)*B,...]每个矩阵在行方向上拼接，vl_hat(A)把A变成一个3*3的反对称矩阵
%       end
%       [U,S,V] = svd(A) ;%U*S*V'=A
%       H{t} = reshape(V(:,9),3,3) ;%按列重新组合成3*3的矩阵
%      
%       
%       P2_ = H{t} * P1 ;
%       du = P2_(1,:)./P2_(3,:) - P2(1,:)./P2(3,:) ;
%       dv = P2_(2,:)./P2_(3,:) - P2(2,:)./P2(3,:) ;
%       ok{t} = (du.*du + dv.*dv) < 6*6 ;
%       score(t) = sum(ok{t}) ;
%     end
    
    for t = 1:100   % 迭代100次，选择最佳的H
      subset = vl_colsubset(1:num_match, 1) ;

      p1=P1(:,subset(1));
      p2=P2(:,subset(1));
      x0=p2(1)-p1(1);
      y0=p2(2)-p1(2);
      H{t} = [1 0 x0; 0 1 y0; 0 0 1];

      P2_ = H{t} * P1 ;
      du = P2_(1,:)./P2_(3,:) - P2(1,:)./P2(3,:) ;
      dv = P2_(2,:)./P2_(3,:) - P2(2,:)./P2(3,:) ;
      ok{t} = (du.*du + dv.*dv) < 20*20 ;
      score(t) = sum(ok{t}) ;
    end

    [score1, best] = max(score) ;
    H = H{best} ;
    ok = ok{best} ;
end

function [new_source_ps, new_aim_ps] = eliminate_tooClose(source_ps, aim_ps, p_dis)
    new_source_ps = source_ps;
    new_aim_ps = aim_ps;
    num = 1;
    for i=1:size(source_ps,2)
        p_i=source_ps(:, i);
        exist_flag=false;
        for j=i+1:size(source_ps,2)
            p_j=source_ps(:,j);
            if ((p_i(1)-p_j(1))^2+(p_i(2)-p_j(2))^2)< p_dis*p_dis
                exist_flag=true;
                break;
            end
        end
        if exist_flag == false
            new_source_ps(:, num) = p_i(1:2);
            new_aim_ps(:, num) = aim_ps(:, i);
            num=num+1;
        end
    end
    num=num-1;
    new_source_ps=new_source_ps(:, 1:num);
    new_aim_ps=new_aim_ps(:, 1:num);
end