function [tomShifted, cumShiftPx, phaseMatShifted] = RegisterBScanBulkMotion(tom, options)
  % RegisterBScanBulkMotion corrects axial and lateral rigid bulk motion
  % between adjacents B-scans (along third dimension) of an OCT tomogram
  %
  % [tomShifted, cumShiftPx] = RegisterBScanBulkMotion(tom, options);
  %
  % Example for options structure: options = struct('upsamplFac', 40,...
  % 'logInt', true, 'onlyAxial', false, 'highPassFilt', false,...
  % 'filtOutliers', false, 'verbosity', false, 'isComplex', true,...
  % 'polChIndx', 4);
  %
  % Inputs:
  % tom: Complex tomogram.
  % options: struct with optinal inputs:
  %   upsamplFac: Oversampling factor for sub-pixel B-scan registration.
  %   logInt: Register B-scans using logarithmic intensity?
  %   onlyAxial: Correct only axial (and not lateral) motion?
  %   latPxWise: Force lateral correction to be pixel (not subpixel) wise
  %   highPassFilt: Size of high-pass-filter to filter low-frequency motion
  %                 and prevent to flatten tissue
  %   filtOutliers: Replace values of atypical values for a mean value.
  %   verbosity: Print current B-scan being registered.
  %   isComplex: Whether the input is complex or not.
  %   polChIndx: Whether the input has polarization channels, correct both
  %              equally.
  %   anchorRef: Whether to anchor the first shift to be the same after
  %              highPassFilt or not. Useful when correcting a volume in a
  %              for loop for maintaining continuity of the shifts.
  %
  % Outputs:
  %   tomShifted: Tomogram after motion correction.
  %   cumShiftPx: Axial and lateral shifts for each B-scan, such that axial
  %               shifts are thisShiftPx(1, 1, :) and lateral shifts are
  %               thisShiftPx(2, 1, :).
  %
  %
  % Authors:  Sebastian Ruiz-Lopera {1*}, Rene Restrepo {1}, B. E. Bouma {2}
  %           and Nestor Uribe-Patarroyo {2}.
  % SRL - RR:
  %	 1. Applied Optics Group, Universidad EAFIT, Carrera 49 # 7 Sur-50,
  %     Medellin, Colombia.
  %
  % BEB - NUP:
  %  2. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA.
  %
  %	* <sruizlo2@eafit.edu.co>
  %
  % MGH-HMS-EAFIT OCT Postprocessing Project
  %
  % Changelog:
  %
  % V1.0 (2020-11-04): Initial version released
  %
  % Copyright Sebastian Ruiz Lopera (2020)
  %
  
  if nargin == 0
    error('At the least one input is required');
  elseif nargin > 1
    % If options are provided, unpack them
    StructToVars(options);
  end
  % First pad tom to avoid wrapping of content due to motion correction
  if ~exist('padding', 'var') || isempty(padding)
    padding = [100, 100];
  end
  
  if (~exist('noiseFloor', 'var') || isempty(noiseFloor)) && any(padding ~= 0)
    error('Need to provide noisefloor for padding')
  end
  
  [nZ, nX, nY, nChs] = size(tom);
  
  if ~exist('upsamplFac', 'var') || isempty(upsamplFac)
    upsamplFac = 20;
  end
  
  if ~exist('logInt', 'var') || isempty(logInt)
    logInt = false;
  end
  
  if ~exist('onlyAxial', 'var') || isempty(onlyAxial)
    onlyAxial = false;
  end
  
  if ~exist('latPxWise', 'var') || isempty(latPxWise)
    latPxWise = false;
  end
  
  if ~exist('highPassFilt', 'var') || isempty(highPassFilt)
    highPassFilt = false;
  end
  
  if ~exist('filtOutliers', 'var') || isempty(filtOutliers)
    filtOutliers = false;
  end
  
  if ~exist('verbosity', 'var') || isempty(verbosity)
    verbosity = false;
  end
  
  if ~exist('isComplex', 'var') || isempty(isComplex)
    isComplex = ~isreal(tom);
  end
  
  if ~exist('polChIndx', 'var') || isempty(polChIndx)
    polChIndx = ndims(tom) + 1;
  end
  
  if ~exist('anchorRef', 'var') || isempty(anchorRef)
    anchorRef = false;
  end
  
  if ~exist('zROI', 'var') || isempty(zROI)
    zROI = 1:nZ;
  end
  
  if ~exist('xROI', 'var') || isempty(xROI)
    xROI = 1:nX;
  end
  
  if ~exist('cumShiftZPx', 'var') || isempty(cumShiftZPx)
    cumShiftZPx = [];
  end
  
  if ~exist('cumShiftXPx', 'var') || isempty(cumShiftXPx)
    cumShiftXPx = [];
  end
  
  if ~exist('phaseMat', 'var') || isempty(phaseMat)
    doCorPhaseMat = false;
  else
    if all(size(phaseMat, [1 2 3]) == size(tom, [1 2 3]))
      doCorPhaseMat = true;
    else
      error("not the right phaseMat size for current tom")
    end
  end
  
  % Calculate FT of intensity tomogram in Z and X and sum along
  % polarization channel dimension
  if isComplex
    if logInt
      tomIntFT = fft(fft(10*log10(sum(abs(tom(zROI, xROI, :, :)).^2, polChIndx)), [], 1), [], 2);
    else
      tomIntFT = fft(fft(sum(abs(tom(zROI, xROI, :, :)).^2, polChIndx), [], 1), [], 2);
    end
  else
    if logInt
      tomIntFT = fft(fft(10*log10(sum(tom(zROI, xROI, :, :), polChIndx)), [], 1), [], 2);
    else
      tomIntFT = fft(fft(sum(tom(zROI, xROI, :, :), polChIndx), [], 1), [], 2);
    end
  end
  
  if isempty(cumShiftZPx) || isempty(cumShiftXPx)
    % Get shift in Z and X with sub-pixel image registration
    shiftZPx = zeros(1, 1, nY);
    shiftXPx = zeros(1, 1, nY);
    if verbosity
      fprintf('Processing Bscan (of %d):\n', nY)
    end
    for thisY = 2:nY
      if verbosity
        fprintf('%d,', thisY)
      end
      [thisShiftPx] = dftregistration(tomIntFT(:, :, thisY - 1), tomIntFT(:, :, thisY), upsamplFac);
      shiftZPx(:, :, thisY) = thisShiftPx(3);
      if onlyAxial
        shiftXPx(:, :, thisY) = 0;
      else
        if latPxWise
          shiftXPx(:, :, thisY) = round(thisShiftPx(4));
        else
          shiftXPx(:, :, thisY) = thisShiftPx(4);
        end
      end
    end
    
    if filtOutliers
      shiftZPx(abs(shiftZPx) > abs(mean(shiftZPx, 3)) + 3*std(shiftZPx, [], 3)) = mean(shiftZPx, 3);
      shiftXPx(abs(shiftXPx) > abs(mean(shiftXPx, 3)) + 3*std(shiftXPx, [], 3)) = mean(shiftXPx, 3);
    end
    
    % Perform accumulative sum to correct with respecto to first BScan
    cumShiftZPx = cumsum(shiftZPx, 3);
    cumShiftXPx = cumsum(shiftXPx, 3);
    % Reference to anchor shifts after highPassFilt
    cumShiftPxRef = cat(1, cumShiftZPx(1), cumShiftXPx(1));
    if any(highPassFilt > -1)
      if numel(highPassFilt) == 1
        highPassFilt = [highPassFilt, highPassFilt];
      end
      if highPassFilt(1) >= 0
        cumShiftZPx = RealHighPassFilter(cumShiftZPx, highPassFilt(1), 3, 'hat');
      end
      if highPassFilt(2) >= 0
        cumShiftXPx = RealHighPassFilter(cumShiftXPx, highPassFilt(2), 3, 'hat');
      end
    end
  else
    % cumShiftZPx and cumShiftXPx provided, need to permute to apply them across
    % 3rd index
    cumShiftZPx = permute(cumShiftZPx, [1 3 2]);
    cumShiftXPx = permute(cumShiftXPx, [1 3 2]);
  end
  clear tomIntFT
  
  % Apply padding if desired
  if any(padding ~= 0) && isComplex
    if ischar(noiseFloor)
      tom = padarray(tom, padding, noiseFloor, 'both');
    else
      tom = padarray(tom, padding, sqrt(noiseFloor), 'both');
    end
    [nZ, nX, ~, ~] = size(tom);
    if doCorPhaseMat
      phaseMat = padarray(phaseMat, padding, 'replicate', 'both');
    end
  elseif any(padding ~= 0) && ~isComplex
    tom = padarray(tom, padding, noiseFloor, 'both');
    [nZ, nX, ~, ~] = size(tom);
    if doCorPhaseMat
      phaseMat = padarray(phaseMat, padding, 'replicate', 'both');
    end
  end
  
  % Create vectors to apply shift
  zVect = ifftshift(-fix(nZ / 2):ceil(nZ / 2) - 1)';
  xVect = ifftshift(-fix(nX / 2):ceil(nX / 2) - 1);
  
  % Anchor shift to reference
  if anchorRef
    cumShiftZPx = cumShiftZPx + cumShiftPxRef(1) - cumShiftZPx(1);
    cumShiftXPx = cumShiftXPx + cumShiftPxRef(2) - cumShiftXPx(1);
  end
  
  % Calculate FT of complex tomogram in Z and X
  if onlyAxial || ~isComplex
    tomFT = fft(tom, [], 1);
    if doCorPhaseMat
      phaseMatFT = fft(phaseMat, [], 1);
    end
  else
    tomFT = fft(fft(tom, [], 1), [], 2);
    if doCorPhaseMat
      phaseMatFT = fft(fft(phaseMat, [], 1), [], 2);
    end
  end
  clear tom
  
  if isComplex
    % Apply shift in Fourier domain as a linear phase term
    if onlyAxial
      tomShifted = ifft(tomFT .* exp(-2i * pi *(cumShiftZPx .* zVect / nZ)), [], 1);
      if doCorPhaseMat
        phaseMatCor = ifft(phaseMatFT .* exp(-2i * pi * (cumShiftZPx .* zVect / nZ)), [], 1, 'symmetric');
      end
    else
      tomShifted = ifft(ifft(tomFT .* exp(-2i * pi * (cumShiftZPx .* zVect / nZ +...
        cumShiftXPx .* xVect / nX)), [], 1), [], 2);
      if doCorPhaseMat
        phaseMatCor = ifft(ifft(phaseMatFT .* exp(-2i * pi * (cumShiftZPx .* zVect / nZ +...
          cumShiftXPx .* xVect / nX)), [], 1), [], 2, 'symmetric');
      end
    end
  else
    % Only in z
    % Apply shift in Fourier domain as a linear phase term
    tomShifted = ifft(tomFT .* exp(-2i * pi * (cumShiftZPx .* zVect / nZ)),...
      [], 1);
    if doCorPhaseMat
      phaseMatCor = ifft(phaseMatFT .* exp(-2i * pi * (cumShiftZPx .* zVect / nZ)),...
        [], 1, 'symmetric');
    end
  end
  
  % Remove padding
  tomShifted = tomShifted(padding(1) + 1:end - padding(1),...
    padding(2) + 1:end - padding(2), :, :);
  if doCorPhaseMat
    phaseMatCor = phaseMatCor(padding(1) + 1:end - padding(1),...
      padding(2) + 1:end - padding(2), :, :);
  end
  
  if nargout > 1
    cumShiftPx = squeeze(cat(1, cumShiftZPx, cumShiftXPx));
  end
  
  if nargout > 2
    if doCorPhaseMat
      phaseMatShifted = phaseMatCor;
    else
      error("did not supply phaseMat")
    end
  end
  
  if verbosity
    fprintf('\nDone!\n')
  end
end
