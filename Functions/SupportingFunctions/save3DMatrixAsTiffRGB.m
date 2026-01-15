function save3DMatrixAsTiffRGB(data, filename, clim, cmap, outSize3, method)
% Save a 3D matrix as an RGB TIFF stack with desired colormap, scale limits,
% and optional 3D resize to exact output size.
%
% data     - 3D matrix (rows x cols x slices)
% filename - output filename (e.g., 'stack.tif')
% clim     - [min max] intensity limits for scaling (like caxis)
% cmap     - colormap matrix (Nx3), e.g., parula(256)
% outSize3 - [rows cols slices] exact output size for imresize3 (optional)
% method   - 'nearest' | 'linear' | 'cubic' (default 'linear')
%
% If outSize3 is empty, size is unchanged.

if nargin < 3 || isempty(clim)
    clim = [min(data(:)) max(data(:))];
end
if nargin < 4 || isempty(cmap)
    cmap = parula(256);
end
if nargin < 5
    outSize3 = []; % no resize
end
if nargin < 6 || isempty(method)
    method = 'linear';
end

% --- Normalize to [0,1] using clim (clip outside range) ---
dataNorm = (data - clim(1)) / max(eps, (clim(2) - clim(1)));
dataNorm = min(max(dataNorm, 0), 1);

% --- 3D resize (if requested) ---
if ~isempty(outSize3)
    vol = imresize3(dataNorm, outSize3, method);
else
    vol = dataNorm;
end

% --- Convert each slice to RGB using the colormap and write TIFF stack ---
nSlices = size(vol, 3);
nColors = size(cmap, 1);
for k = 1:nSlices
    slice = vol(:,:,k);
    sliceRGB = ind2rgb(gray2ind(slice, nColors), cmap);

    if k == 1
        imwrite(sliceRGB, filename, 'tif', 'Compression', 'none', 'WriteMode', 'overwrite');
    else
        imwrite(sliceRGB, filename, 'tif', 'Compression', 'none', 'WriteMode', 'append');
    end
end
end