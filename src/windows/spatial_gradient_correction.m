
function [gm, gamma] = spatial_gradient_correction(data, r, sigma, varargin)

p = inputParser;

addParameter(p, 'type', 'points',@(x) any(validatestring(x,{'points', 'pairs'})))
parse(p, varargin{:})
type = p.Results.type;

d = r(2)-r(1);

left = min(data(1).spacewin.x);
right = max(data(1).spacewin.x);
top = min(data(1).spacewin.y);
bottom = max(data(1).spacewin.y);

I1 = histcounts2(data(1).x, data(1).y, left:d:right, top:d:bottom)';
if numel(data)==1
    I2 = I1;
elseif numel(data) == 2
    I2 = histcounts2(data(2).x, data(2).y, left:d:right, top:d:bottom)';
else
    error('only supported for 1 or 2 colors')
end

[X, Y] = meshgrid(left+d/2:d:right-d/2, top+d/2:d:bottom-d/2);
mask = inpolygon(X,Y,data(1).spacewin.x, data(1).spacewin.y);

lx = 2*size(I1, 1)+1;
ly = 2*size(I1, 2)+1;

switch type
    case 'points'
        Iblur1 =gausblur(I1.*mask, sigma/d).*mask;
        Iblur2 =gausblur(I2.*mask, sigma/d).*mask;
        Mblur =gausblur(mask, sigma/d).*mask;
        
        Iblur1 = Iblur1./Mblur;
        Iblur1(~Mblur) = 0;
        Iblur2 = Iblur2./Mblur;
        Iblur2(~Mblur) = 0;
        
        I1_fudgefact = mean(Iblur1(logical(mask(:))))/mean(I1(logical(mask(:))));
        I1_densfact = mean(Iblur1(logical(mask(:))));
        I2_fudgefact = mean(Iblur2(logical(mask(:))))/mean(I2(logical(mask(:))));
        I2_densfact = mean(Iblur2(logical(mask(:))));
        
        gamma = fftshift(ifft2(fft2(Iblur1, lx, ly).*conj(fft2(Iblur2, lx, ly))))./ ...
            fftshift(ifft2(abs(fft2(mask, lx, ly)).^2))/...
            I2_densfact/I1_densfact*I1_fudgefact*I2_fudgefact;
    
    case 'pairs'
        xc = ifft2(fft2(I1, lx, ly).*conj(fft2(I2, lx, ly)));
        if numel(data)==1
            xc(1,1) = xc(1,1) - sum(I1(:));
        end
        xc = fftshift(xc);
        
        m2 = fftshift(ifft2(abs(fft2(mask, lx, ly)).^2));
        m2(m2 < 1) = 0; % prevent some problems with rounding errors later.
        mask2 = double(logical(m2));
        edgem2 = gausblur(mask2, sqrt(2)*sigma/d).*mask2;
        
        rat = xc ./ m2;
        rat(~mask2) = 0;
        ratblur = gausblur(rat, sqrt(2)*sigma/d);
        
        I1_densfact = mean(I1(logical(mask(:))));
        I2_densfact = mean(I2(logical(mask(:))));
        
        
        gamma = ratblur./edgem2/I2_densfact/I1_densfact; %.* m2./edgem2/d^2;
        gamma(~mask2) = 0;
end

%%
% angular average
[X, Y] = meshgrid((1:size(gamma, 2))-size(gamma, 2)/2-1, (1:size(gamma, 1))-floor(size(gamma, 1)/2)-1);
R = d*sqrt(X.^2+Y.^2);
gm = zeros(size(r));
for i=1:numel(r)
    inds = R(:)>=r(i)-d/2 &R(:)<r(i)+d/2;
    gm(i) = mean(gamma(inds));
end










% 
% 
% 
% 
% 
% 
% 
% 
% 
% res = r(2)-r(1);
% 
% 
% dims = [256 512]; 
% mask = res_specs.mask;
% 
% % filter for blur
% g_filt = fspecial('gaussian',300, 32);
% 
% movies = 1:numel(alignment_data1.alldata_track);
% data_slice = alignment_data1.alldata_raw(movies); 
% [scrap Ctot_raw1]=generate_STORM_image_MBS(data_slice, res, s, dims);
% I1 = Ctot_raw1;
% 
% data_slice = alignment_data2.alldata_raw(movies); 
% [scrap Ctot_raw2]=generate_STORM_image_MBS(data_slice, res, s, dims);
% I2 = Ctot_raw2;
% 
% % I1_f= (imfilter(mask.*I1, g_filt));
% % I2_f= (imfilter(mask.*I2, g_filt));
% 
% I1_f= mask.*(imfilter(mask.*I1, g_filt))./imfilter(double(mask), g_filt);
% I1_f(isnan(I1_f(:))) = 0;
% I1_fudgefact = mean(I1_f(logical(mask(:))))/mean(I1(logical(mask(:))));
% I1_densfact = mean(I1_f(logical(mask(:))));
% 
% I2_f= mask.*(imfilter(mask.*I2, g_filt))./imfilter(double(mask), g_filt);
% I2_f(isnan(I2_f(:))) = 0;
% I2_fudgefact = mean(I2_f(logical(mask(:))))/mean(I2(logical(mask(:))));
% I2_densfact = mean(I2_f(logical(mask(:))));
% 
% L1 = size(I1_f, 1);
% L2 = size(I1_f, 2);
% 
% N4 = (1/I2_densfact)*(1/I1_densfact)*fftshift(ifft2(fft2(I1_f, 2*L1+1, 2*L2+1).*conj(fft2(I2_f, 2*L1+1, 2*L2+1))))*I1_fudgefact.*I2_fudgefact;
% 
% normalize = sum(sum(mask));
% 
% NP_norm = N4/normalize;
% 
% [~, rad_ave_NP_edge]  = radial_average_EM_int(NP_norm,ceil(rmax/res),(rbinsize/res),0);
% XCor_norm = rad_ave_NP_edge(1:end-1);
% 
