function [result_img] = auto_point(raw_image, pad_radio, l_num)
run('vlfeat-0.9.21-bin\vlfeat-0.9.21\toolbox\vl_setup.m');
%raw_image=imread('img1.jpg');
%new_img = pad_resize(raw_image, 0.1);
new_image = imresize(raw_image, 0.4);
new_image = rgb2gray(new_image);

G = fspecial('gaussian', [25 25], 15);
Ig = imfilter(new_image,G,'same');
img_single=single(Ig);
[f,d] = sift(img_single) ;
points = f(1:2,:);

p_grey = Ig(round(points(2,:)),round(points(1,:)));
p_grey = diag(p_grey);
p_grey=double(p_grey);
points_grey=[points;p_grey'];

new_points = [];
for i=1:size(points_grey,2)  %去除白点
   p = points_grey(:,i);
   if points_grey(3,i)<190
       new_points(:,size(new_points,2)+1)=p;
   end
end

[~,sort_index]=sort(new_points(3,:));
new_points=new_points(:,sort_index); %灰度值较小的排在前面，优先被选到

dis_points = []; %去除距离过近的点
for i=1:size(new_points,2)
   p = new_points(:,i);
   exist_flag=false;
   for j = 1:size(dis_points,2)
       p_j = dis_points(:,j);
       if get_dis(p(1:2),p_j(1:2))<30
           exist_flag=true;
           break;
       end
   end
   if exist_flag==false
       dis_points(:,size(dis_points,2)+1)=p;
   end
end

% imshow(Ig);
% hold on;
% plot( dis_points(1,:),dis_points(2,:),'g.' );

p_struct_1=struct('p',0, 'right_neibor',0, 'left_neibor',0, 'right_arc', 0);
points_rightN=repmat(p_struct_1,[size(dis_points,2) 1]);
for i=1:size(dis_points,2)
    p = dis_points(1:2,i);
    right_neibor = [];
    for j=1:size(dis_points,2)
        p_j = dis_points(1:2,j);
        if ~(p_j(1)>p(1) && p_j(1)<p(1)+140) %找到在p右边一个方框内的点
            continue
        end
        if ~(p_j(2)>p(2)-70 && p_j(2)<p(2)+70)
            continue
        end
        right_neibor(1,size(right_neibor,2)+1)=j;
        arc=atan((p_j(2)-p(2))/(p_j(1)-p(1)))*180/pi;
        right_neibor(2,size(right_neibor,2))=arc;
    end
    
    points_rightN(i).p = i;
    
    if size(right_neibor,1)==0
        continue;
    end
    
    [~,index_pn]=min(abs(right_neibor(2,:)));
    right_neibor_p = right_neibor(1, index_pn);
    arc = right_neibor(2, index_pn);
    
    if size(right_neibor_p, 1) ~= 0 && abs(arc)<15 %角度不能太大
        points_rightN(i).right_neibor = right_neibor_p(1);
        points_rightN(i).right_arc = arc;
        points_rightN(right_neibor_p(1)).left_neibor = i;
    else
        points_rightN(i).right_neibor = 0;
    end
    
end

% figure(2);
% imshow(Ig);
% hold on;
% 
% for i=1:size(points_rightN, 1)
%     p_struct = points_rightN(i);
%     p = dis_points(1:2, p_struct.p);
%     if p_struct.right_neibor ~= 0
%         p_neibor = dis_points(1:2, p_struct.right_neibor);
%         line([p(1),p_neibor(1)], [p(2),p_neibor(2)]) ;
%     end
%     plot( p(1), p(2),'g.' );
% end

line_struct_1=struct('ps',[], 'length',0);
line_structN=repmat(line_struct_1,[round(size(dis_points,2)/8) 1]);
line_len = 1;
for i=1:size(dis_points,2)
    p = dis_points(1:2,i);
    if p(1)>size(Ig, 2)/4  %寻找起始点
        continue;
    end
    if points_rightN(i).left_neibor ~=0
        continue;
    end
    
    lines_p = [];
    lines_p(1)=i;
    last_arc = 0;
    delta_arc_flag = false;
    while true
        next = points_rightN(lines_p(end)).right_neibor;
        arc_delta = points_rightN(lines_p(end)).right_arc-last_arc;
        last_arc=points_rightN(lines_p(end)).right_arc;
        if abs(arc_delta)>12
            delta_arc_flag = true;
            break;
        end
        if next==0
            break
        end
        lines_p(length(lines_p)+1)=next;
    end
    
    if length(lines_p)<8  %曲线长度不够
        continue;
    end
    line_structN(line_len).ps = lines_p;
    line_structN(line_len).length = length(lines_p);
    line_len = line_len+1;
end
line_len = line_len-1;

[~, line_index]=sort([line_structN.length]);
line_structN = line_structN(line_index);

line_struct_1=struct('ps',[], 'length',0);
new_line_structN=repmat(line_struct_1,[line_len 1]);
new_line_len = 1;
for i=1:size(line_structN,1)
    ps_i = line_structN(i).ps;
    common = false;
    for j=i+1:size(line_structN,1)
        ps_j = line_structN(j).ps;
        if ~isempty(intersect(ps_i, ps_j))
            common = true;
            break;
        end
    end
    if common==true
        continue;
    end
    new_line_structN(new_line_len) = line_structN(i);
    new_line_len = new_line_len+1;
end
new_line_len = new_line_len-1;

[~, new_line_index]=sort([new_line_structN.length], 'descend');
new_line_structN = new_line_structN(new_line_index);
%l_num = 7;%11;
new_line_structN = new_line_structN(1:l_num);
new_line_len=l_num;

% figure(3);
% imshow(Ig);
% hold on;
% 
% for i=1:new_line_len
%     ps = new_line_structN(i).ps;
%     for p_n=1:length(ps)-1
%         p1 = dis_points(1:2, ps(p_n));
%         p2 = dis_points(1:2, ps(p_n+1));
%         line([p1(1),p2(1)], [p1(2),p2(2)]) ;
%         plot( p1(1), p1(2),'g.' );
%         plot( p2(1), p2(2),'g.' );
%     end
% end

select_ps = []; %得到原来的点和矫正后的点
aim_ps = [];
for i=1:new_line_len
    ps = new_line_structN(i).ps;
    p_i = size(aim_ps, 1)+1;
    aim_ps(p_i,:) = dis_points(1:2, ps(1));
    select_ps(p_i,:) = dis_points(1:2, ps(1));
    for p_n=2:length(ps)
        p_i = p_i+1;
        p_select = dis_points(1:2, ps(p_n));
        select_ps(p_i,:) =p_select;
        dis = get_dis(p_select, select_ps(p_i-1,:));
        aim_ps(p_i,1) = aim_ps(p_i-1, 1)+dis;
        aim_ps(p_i,2) = aim_ps(p_i-1, 2);
    end
end

new_image_pad = padarray(new_image, [round(size(Ig, 1)*pad_radio) round(size(Ig, 2)*pad_radio)], 'both');
select_ps(:,1) = select_ps(:,1) +  round(size(Ig, 2)*pad_radio);
select_ps(:,2) = select_ps(:,2) +  round(size(Ig, 1)*pad_radio);
aim_ps(:,1) = aim_ps(:,1) +  round(size(Ig, 2)*pad_radio);
aim_ps(:,2) = aim_ps(:,2) +  round(size(Ig, 1)*pad_radio);

% figure(4);
% imshow(new_image_pad);
% hold on;
% plot( select_ps(:,1),select_ps(:,2),'r.' );
% plot( aim_ps(:,1),aim_ps(:,2),'g.' );

[imo1,mask1] = rbfwarp2d( new_image_pad, select_ps, aim_ps,'thin');
% imshow(uint8(imo1));

delta_x = round(size(imo1, 1)*pad_radio/2);
delta_y = round(size(imo1, 2)*pad_radio/2);
cut_imo1(:,:)=imo1(delta_x: size(imo1, 1)-delta_x, delta_y: size(imo1, 2)-delta_y);
result_img=uint8(cut_imo1);

end

function dis=get_dis(p1,p2)
    dis=sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2);
end