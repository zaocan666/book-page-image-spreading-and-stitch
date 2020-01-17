function cut_imo1 = undistort(select_ps, raw_img, row, col, index, pad_radio)

[select_ps, aim_p] = get_aim_p(select_ps, raw_img, row, col, index);

% subplot(1,2,1);
imshow(raw_img);hold on;
scatter( select_ps(:,1),select_ps(:,2), 80, 'r.' );
scatter( aim_p(:,1),aim_p(:,2), 80, 'g.' );

%close all
[imo1,mask1] = rbfwarp2d( raw_img, select_ps, aim_p,'thin');
delta_x = round(size(imo1, 1)*pad_radio/2);
delta_y = round(size(imo1, 2)*pad_radio/2);
cut_imo1(:,:)=imo1(delta_x: size(imo1, 1)-delta_x, delta_y: size(imo1, 2)-delta_y);

% subplot(1,2,2);
% imshow(uint8(cut_imo1));

end

function [select_ps, aim_p] = get_aim_p(select_ps, raw_img, row, col, index)
    close all;
    %imshow(raw_img);
    %select_ps=ginput(row*col);
    
    if index==1
        center = round(8*col/10);
    else
        center = 1;
    end
    
    p=[0 0];
    aim_p = repmat(p, [row*col 1]);
   
    for i=1:row
        aim_p((i-1)*col+center,:)=select_ps((i-1)*col+center,:);
        for j=center:-1:2
            p1=select_ps((i-1)*col+j,:);
            p2=select_ps((i-1)*col+j-1,:);
            dis=get_dis(p1,p2);
            extra = 0;
            if index==2 && j==2
                extra=-40;
            end
            aim_p((i-1)*col+j-1,1)=aim_p((i-1)*col+j,1)-dis-extra;
            if aim_p((i-1)*col+j-1,1)<1
                aim_p((i-1)*col+j-1,1)=1;
            end
            aim_p((i-1)*col+j-1,2)=aim_p((i-1)*col+j,2);
        end
        
        for j=center:col-1
            p1=select_ps((i-1)*col+j,:);
            p2=select_ps((i-1)*col+j+1,:);
            dis=get_dis(p1,p2);
            extra = 0;
            if index==2 && j==1
                extra=40;
            end
            aim_p((i-1)*col+j+1,1)=aim_p((i-1)*col+j,1)+dis + extra;
            if aim_p((i-1)*col+j+1,1)>size(raw_img, 1)
                aim_p((i-1)*col+j+1,1)=size(raw_img, 1);
            end
            aim_p((i-1)*col+j+1,2)=aim_p((i-1)*col+j,2);
        end
    end
end

function dis=get_dis(p1,p2)
dis = sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2);
end