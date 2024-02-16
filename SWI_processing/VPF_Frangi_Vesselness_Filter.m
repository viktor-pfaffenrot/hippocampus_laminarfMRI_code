function V = VPF_Frangi_Vesselness_Filter(I,m,s)
sz = size(I);

if ~islogical(m)
    if any(isnan(m))
        m(isnan(m)) = 0;
    end
    m = logical(m);
end
Iw = I.*m;

if any(isnan(Iw))
    Iw(isnan(Iw)) = 0;
end
%%
H = s^2.*VPF_Frangi_Hessian(Iw,sz,s);
E = VPF_Frangi_Eigs(H,m,sz);

alpha = .5;
beta = .5;
c = 0.1*median(Iw(m));
tmp = VPF_Frangi_Vesselness(E,alpha,beta,c);

V = zeros(sz);
V(m) = tmp;
end

function H = VPF_Frangi_Hessian(I,sz,s)

%define Gaussian filter kernel
GIs = imgaussian(I,s);
Hxx = gradient(gradient(GIs));

Hx = gradient(GIs);
[~,Hxy,Hxz] = gradient(Hx);
[~,Hy] = gradient(GIs);
[~,Hyy,Hyz] = gradient(Hy);
[~,~,Hz] = gradient(GIs);
[~,~,Hzz] = gradient(Hz);

Hxx = Hxx(:); Hxy = Hxy(:); Hxz = Hxz(:);
Hyy = Hyy(:); Hyz = Hyz(:);
Hzz = Hzz(:);

H = cat(3,[Hxx Hxy Hxz],[Hxy Hyy Hyz],[Hxz Hyz Hzz]);
H = reshape(H, [sz 3 3]);
end

function I=imgaussian(I,sigma,siz)
% IMGAUSSIAN filters an 1D, 2D color/greyscale or 3D image with an 
% Gaussian filter. This function uses for filtering IMFILTER or if 
% compiled the fast  mex code imgaussian.c . Instead of using a 
% multidimensional gaussian kernel, it uses the fact that a Gaussian 
% filter can be separated in 1D gaussian kernels.
%
% J=IMGAUSSIAN(I,SIGMA,SIZE)
%
% inputs,
%   I: The 1D, 2D greyscale/color, or 3D input image with 
%           data type Single or Double
%   SIGMA: The sigma used for the Gaussian kernel
%   SIZE: Kernel size (single value) (default: sigma*6)
% 
% outputs,
%   J: The gaussian filtered image
%
% note, compile the code with: mex imgaussian.c -v
%
% example,
%   I = im2double(imread('peppers.png'));
%   figure, imshow(imgaussian(I,10));
% 
% Function is written by D.Kroon University of Twente (September 2009)
if(~exist('siz','var')), siz=sigma*6; end
if(sigma>0)
    % Make 1D Gaussian kernel
    x=-ceil(siz/2):ceil(siz/2);
    H = exp(-(x.^2/(2*sigma^2)));
    H = H/sum(H(:));
    % Filter each dimension with the 1D Gaussian kernels\
    if(ndims(I)==1)
        I=imfilter(I,H, 'same' ,'replicate');
    elseif(ndims(I)==2)
        Hx=reshape(H,[length(H) 1]);
        Hy=reshape(H,[1 length(H)]);
        I=imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
    elseif(ndims(I)==3)
        if(size(I,3)<4) % Detect if 3D or color image
            Hx=reshape(H,[length(H) 1]);
            Hy=reshape(H,[1 length(H)]);
            for k=1:size(I,3)
                I(:,:,k)=imfilter(imfilter(I(:,:,k),Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
            end
        else
            Hx=reshape(H,[length(H) 1 1]);
            Hy=reshape(H,[1 length(H) 1]);
            Hz=reshape(H,[1 1 length(H)]);
            I=imfilter(imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate'),Hz, 'same' ,'replicate');
        end
    else
        error('imgaussian:input','unsupported input dimension');
    end
end
end

function E = VPF_Frangi_Eigs(H,m,sz)
% keyboard
H = reshape(permute(H,[4 5 1 2 3]),[3 3 prod(sz)]);
H = H(:,:,m);
E = zeros(3,size(H,3));
for ii = 1:size(H,3)
    E(:,ii) = eig(H(:,:,ii));
end
[~,idx] = sort(abs(E));
for ii = 1:size(H,3)
    E(:,ii) = E(idx(:,ii),ii);
end
end

function V = VPF_Frangi_Vesselness(E,alpha,beta,c)

Ra = abs(E(2,:))./abs(E(3,:)); 
Rb = abs(E(1,:))./sqrt(abs(E(2,:).*E(3,:)));
S  = sqrt(sum(E.^2,1));
V = (1 - exp(-(Ra.^2./(2*alpha^2)))) .* exp(-(Rb.^2./(2*beta^2))) .* (1 - exp(-(S.^2./(2*c^2))));
V(E(2,:) < 0) = 0;
V(E(3,:) < 0) = 0;

end
