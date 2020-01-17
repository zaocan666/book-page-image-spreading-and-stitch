function [imo,mask] = rbfwarp2d( img, p_source, p_dest, varargin )
% 来自：MATLAB官网 https://ww2.mathworks.cn/matlabcentral/fileexchange/53828-rbf-or-thin-plate-splines-image-warping
% 经修改
% Radial base function/Thin-plate spline 2D image warping.
% [imo,mask] = rbfwarp2d( im, ps, pd, method)
%   input:
%       im: image 2d matrix
%       ps: 2d source landmark [n*2]
%       pd: 2d destin landmark [n*2]
%       method:
%         'gau',r  - for Gaussian function   ko = exp(-|pi-pj|/r.^2);
%         'thin'   - for Thin plate function ko = (|pi-pj|^2) * log(|pi-pj|^2)
%   output:
%       imo  : output matrix
%       mask : mask for output matrix, 0/1 means out/in the border
%
%   Bookstein, F. L. 
%   "Principal Warps: Thin Plate Splines and the Decomposition of Deformations."
%   IEEE Trans. Pattern Anal. Mach. Intell. 11, 567-585, 1989. 
%
%   Code by WangLin
%   2015-11-5
%   wanglin193@hotmail.com

num_required_parameters = 3;
if nargin < num_required_parameters
    help rbfwarp2d.m
    return;
end

% initialize default parameters
[imh,imw,imc] = size(img);
r = 0.1*imw;
imo = zeros(imh,imw,imc);
mask = zeros(imh,imw);

% parse parameters
if nargin > num_required_parameters
    iVarargin = 1;
    while iVarargin <= nargin - num_required_parameters
        switch lower(varargin{iVarargin})
            case 'thin'
                method = 't';
            case 'gau'
                method = 'g';
                r = varargin{iVarargin+1};
                iVarargin = iVarargin + 1;
        end
        iVarargin = iVarargin + 1;
    end
end

%% Training w with L
nump = size(p_dest,1);
num_center = size(p_source,1);
K=zeros(nump,num_center);

for i=1:num_center
    %Inverse warping from destination!
    %dx = ones(nump,1)*ps(i,:)-pd; 
    dx = ones(nump,1)*p_dest(i,:)-p_dest; 
    K(:,i) = sum(dx.^2,2);
end

if( strcmpi(method,'g') )
    K = rbf(K,r);
elseif( strcmpi(method,'t') )
    K = ThinPlate(K);
end

% P = [1,xp,yp] where (xp,yp) are n landmark points (nx2)
P = [ones(num_center,1),p_dest];
% L = [ K  P;
%       P' 0 ]
L = [K,P;P',zeros(3,3)];
% Y = [x,y;
%      0,0]; (n+3)x2
Y = [p_source;zeros(3,2)];
%w = inv(L)*Y;
w = L\Y;

%% Using w
[x,y] = meshgrid(1:imw,1:imh);
pt = [x(:), y(:)];

nump = size(pt,1);
Kp = zeros(nump,num_center);
for i=1:num_center
    %dx = ones(nump,1)*ps(i,:)-pt;
    dx = ones(nump,1)*p_dest(i,:)-pt;
    Kp(:,i) = sum(dx.^2,2);
end
if( strcmpi(method,'g') )
    Kp = rbf(Kp,r);
elseif( strcmpi(method,'t') )
    Kp = ThinPlate(Kp);    
end

L = [Kp,ones(nump,1),pt];
ptall = L*w;

%reshape to 2d image
xd = reshape( ptall(:,1),imh,imw );
yd = reshape( ptall(:,2),imh,imw );

for i = 1:imc
    imt= interp2( single(img(:,:,i)),xd,yd,'linear');
    imo(:,:,i) = uint8(imt);
end


mask = ~isnan(imt);
end

function ko = rbf(d,r) 
    ko = exp(-d/r.^2);
end

function ko = ThinPlate(ri)
% k=(r^2) * log(r^2)
    r1i = ri;
    r1i((ri==0))=realmin; % Avoid log(0)=inf
    ko = (ri).*log(r1i);
end
