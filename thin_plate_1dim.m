function [result_img] = thin_plate_1dim( img, p_source, p_y)
% 薄板样条插值
% img:输入图像
% p_source：输入点
% p_y：输入点对应的值

nump = size(p_source,1);
K=zeros(nump,nump);

for i=1:nump
    %Inverse warping from destination!
    %dx = ones(nump,1)*ps(i,:)-pd; 
    dx = ones(nump,1)*p_source(i,:)-p_source; 
    K(:,i) = sum(dx.^2,2);
end

K = ThinPlate_U(K);

% P = [1,xp,yp] where (xp,yp) are n landmark points (nx2)
P = [ones(nump,1),p_source];
% L = [ K  P;
%       P' 0 ]
L = [K,P;P',zeros(3,3)];
% Y = [x,y;
%      0,0]; (n+3)x2
Y = [p_y;zeros(3,1)];
%w = inv(L)*Y;
w = L\Y;

[imh,imw,~] = size(img);
%% Using w
[x,y] = meshgrid(1:imw,1:imh);
pt = [x(:), y(:)];

all_nump = size(pt,1);
Kp = zeros(all_nump,nump);
for i=1:nump
    %dx = ones(nump,1)*ps(i,:)-pt;
    dx = ones(all_nump,1)*p_source(i,:)-pt;
    Kp(:,i) = sum(dx.^2,2);
end

Kp = ThinPlate_U(Kp);    

L = [Kp,ones(all_nump,1),pt];
ptall = L*w;

%reshape to 2d image
result_img = reshape( ptall,imh,imw );

end

function ko = ThinPlate_U(ri)
% k=(r^2) * log(r^2)
    r1i = ri;
    r1i((ri==0))=realmin; % Avoid log(0)=inf
    ko = (ri).*log(r1i);
end