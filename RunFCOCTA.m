function [alphaA1C, maskSNR, alphaA1CSNR] = RunFCOCTA( ...
  tomInput, backgroundInput, nPolChs, nReps, zNoiseBackground, ...
  zEnsembleWindowHalfSize, zEnsembleWindowExp2Diameter, ...
  xEnsembleWindowHalfSize, xEnsembleWindowExp2Diameter, ...
  yEnsembleWindowHalfSize, yEnsembleWindowExp2Diameter, ...
  displayDaspect, angioMap, saveFolder, tauDiff, verbose)
% RUNFCOCTA process complex-valued reconstructed raw OCT tomogram
% (irrespective of phase stability) with fully coherent methods to produce
% OCTA volumes (and optionally export the volumes to TIFF stacks)
% ------------------------------------------------------------------------
% Authors:  Tianhui (Cindy) Jie
%
% TJ: 
% 1. Massachusetts Institute of Technology, 
%    Institute for Medical Engineering and Sciences, Cambridge, MA, USA;
% 2. Wellman Center for Photomedicine, Harvard Medical School, 
%    Massachusetts General Hospital, Boston, MA, USA;
% <tjie@mit.edu>
%
% FC-OCTA (v1.0)
%
% Changelog:
%
% 1.0 (2026-01-16): Initial version released
%
% Copyright Tianhui (Cindy) Jie (2026)

arguments (Input)
  tomInput {eitherFileLocOr4DTom}            % file location of raw tomogram, OR 4D matrix of raw tomogram
  backgroundInput {eitherFileLocOr4DTom}     % file location of background measurement, OR 4D matrix of background measurement
  nPolChs (1, 1) double                      % number of polarization channels
  nReps (1, 1) double                        % number of B-scan repetitions
  zNoiseBackground (1, :) double             % z ROI above tissue in raw tomogram
  zEnsembleWindowHalfSize (1, 1) double      % 2*(this)+1 = z averaging window size
  zEnsembleWindowExp2Diameter (1, 1) double  % exp(-2)-diameter for z kernel shape
  xEnsembleWindowHalfSize (1, 1) double      % 2*(this)+1 = x averaging window size
  xEnsembleWindowExp2Diameter (1, 1) double  % exp(-2)-diameter for x kernel shape
  yEnsembleWindowHalfSize (1, 1) double      % 2*(this)+1 = y averaging window size
  yEnsembleWindowExp2Diameter (1, 1) double  % exp(-2)-diameter for y kernel shape
  displayDaspect (1, 3) double               % daspect for export, size of each voxel
  angioMap (256, 3) double                   % colormap for OCTA
  saveFolder (1, 1) string                   % file location to save TIFF stacks
  tauDiff (1, 1) double = 1                  % B-scan step between repetitions used to calculate OCTA
  verbose (1, 1) logical = true              % (if true) print out progress markers 
end

arguments (Output)
  alphaA1C  % OCTA volume based on A1C: 4D matrix [nZ (depth), nX (slow axis), nY (fast axis), nReps-1 (all combos of tauDiff)]
  maskSNR   % SNR mask calculated in FC-OCTA that is used on output TIFF stacks: 3D matrix [nZ, nX, nY]
  alphaA1CSNR  % OCTA volume based on A1CSNR: 4D matrix [nZ (depth), nX (slow axis), nY (fast axis), nReps-1 (all combos of tauDiff)]
end
  
  % --------------------------------------------------------------------
  % Load tomogram
  % --------------------------------------------------------------------
  if verbose
    fprintf('Loading your dataset for FC-OCTA processing... \n')
  end
  
  if isstring(tomInput) || ischar(tomInput)
    tomFile = struct2cell(load(tomInput));
    tom = tomFile{1};
  elseif isnumeric(tomInput) && ndims(tomInput) == 4
    tom = tomInput;
  end
  [nZ, nX, nBscans, ~] = size(tom);
  nY = nBscans / nReps;
  if verbose 
    fprintf(['Your dataset is %dx%dx%d [ZxXxY] and has \n' ...
             '  %d repetition(s), %d polarization channel(s). \n'], ...
            nZ, nX, nY, nReps, nPolChs)
  end

  % --------------------------------------------------------------------
  % Determine depth-dependent and approximate noise floors
  % --------------------------------------------------------------------
  if verbose
    fprintf('Starting noise floor determination... \n')
  end
  
  if isstring(backgroundInput) || ischar(backgroundInput)
    backgroundFile = struct2cell(load(backgroundInput));
    background = backgroundFile{1};
  elseif isnumeric(backgroundInput) && ndims(backgroundInput) == 4
    background = backgroundInput;
  end
  backgroundInt = mean(abs(background) .^ 2, [2 3]);
  
  % Remove weird peaks
  backgroundIntFiltered = cat(2, ...
    single(medfilt1(double(backgroundInt(:, 1)), 91, 'truncate')), ...
    single(medfilt1(double(backgroundInt(:, 2)), 91, 'truncate')));
  
  % Determine scaling factors
  darkFactor(1) = mean(squeeze( ...
    mean(abs(tom(zNoiseBackground, :, 1:10, 1)) .^ 2, [2 3])) ...
    ./ backgroundInt(zNoiseBackground, 1));
  darkFactor(2) = mean(squeeze( ...
    mean(abs(tom(zNoiseBackground, :, 1:10, 2)) .^ 2, [2 3])) ...
    ./ backgroundInt(zNoiseBackground, 2));
  
  % Concatenate polarization channels, get ROI, and move polarization
  % states to the 5th index
  noiseFloorZ = cat(5, backgroundIntFiltered(:, 1) * darkFactor(1), ...
    backgroundIntFiltered(:, 2) * darkFactor(2));
  
  clear backgroundInt backgroundIntFiltered darkFactor
  % Check noiseFloorZ
  if verbose
    figure(1), clf
    set(gcf, 'Name', 'depth-dependent noise floor')
    title(sprintf('your depth-dependent noise floor'))
    hold on
    plot(10 * log10(FlattenArrayTo2D( ...
      mean(abs(tom(:, :, 1:10, :)) .^ 2, [2 3])))) 
    plot(10*log10(noiseFloorZ(:, :)))
    hold off
    legend('dataset reference: pol. ch. 1', 'dataset reference: pol. ch. 2', ...
      'noise floor: pol. ch. 1', 'noise floor: pol. ch. 2')
  end
  
  % Depth average for approximate noise floor
  noiseFloorApprox = mean(noiseFloorZ, 'all');
  if verbose
    fprintf(['Noise floor determination done! \n' ...
             '  Approximate noise floor for your dataset is %f [dB] \n'], ...
            10 * log10(noiseFloorApprox));
  end
  
  % --------------------------------------------------------------------
  % Pixel-Wise Bulk B-scan Motion Correction
  % --------------------------------------------------------------------
  if verbose
    fprintf('Starting bulk motion correction... \n')
  end
  
  bulkMCSettings = struct(...
    'upsamplFac', 40, 'logInt', true, 'onlyAxial', false, 'isComplex', true, ...
    'polChIndx', 4, 'highPassFilt', [1 1], 'latPxWise', true);
  [tomMCed, cummulativeShift] = RegisterBScanBulkMotion( ...
    tom, catstruct(bulkMCSettings, struct('noiseFloor', noiseFloorApprox)));
  clear tom
  
  if verbose
    fprintf('  "tom" motion corrected and cleared -> "tomMCed" \n')
  end
  
  % the depth-dependent noise floor is shifted accordingly for each Bscan
  noiseFloorZMCed = reshape(ApplyAxialShift( ...
    repmat(noiseFloorZ, [1, 1, nY * nReps]), ...
    permute(cummulativeShift(1, :), [1, 3, 2])), [nZ, nReps, 1, nY, nPolChs]);
  
  clear cummulativeShift noiseFloorZ
  if verbose
    fprintf('Bulk motion correction done! \n')
  end
  
  % --------------------------------------------------------------------
  % Temporal relative phase-noise correction
  % --------------------------------------------------------------------
  if verbose
    fprintf('Starting phase-noise correction along repetition... \n')
  end
  
  % move repetitions to the 2nd index
  tomMCedReshaped = permute(reshape(tomMCed, ...
    [nZ, nX, nReps, nY, nPolChs]), [1, 3, 2, 4, 5]);
  clear tomMCed
  
  % fix phase jitter along repetitions, both slopes and offset
  tomPNCed = ApplyGeneralAlinePhaseJitterCorrection( ...
    tomMCedReshaped, 1:nZ, 10 * log10(noiseFloorApprox), true, 5, false);
  clear tomMCedReshaped
  
  if verbose
    fprintf('  "tomMCed" phase-noise corrected and cleared -> "tomPNCed" \n')
    fprintf('Phase-noise correction along repetitions done! \n')
  end
  
  % --------------------------------------------------------------------
  % Calculate ensemble averaging kernel
  % --------------------------------------------------------------------
  % create Gaussian kernel in z, x, y directions
  zEnsembleWindow = AnisotropicGaussianExp2Diameter(...
    [1, 2 * zEnsembleWindowHalfSize + 1], ...
    0, zEnsembleWindowExp2Diameter);
  xEnsembleWindow = AnisotropicGaussianExp2Diameter(...
    [2 * xEnsembleWindowHalfSize + 1, 1], ...
    xEnsembleWindowExp2Diameter, 0);
  yEnsembleWindow = AnisotropicGaussianExp2Diameter(...
    [2 * yEnsembleWindowHalfSize + 1, 1], ...
    yEnsembleWindowExp2Diameter, 0);
  yEnsembleWindow = permute(yEnsembleWindow, [1 3 2]);
  
  % calculate kernels
  axialWindow = zEnsembleWindow;
  lateralWindow = xEnsembleWindow .* yEnsembleWindow;
  ensembleWindow = axialWindow .* lateralWindow;
  
  % check kernel sizes
  effectiveAxialSize = sum(axialWindow ./ ...
    max(axialWindow, [], 'all'), 'all');
  effectiveLateralSize = sum(lateralWindow ./ ...
    max(lateralWindow, [], 'all'), 'all');
  effectiveKernelSize = sum(ensembleWindow ./ ...
    max(ensembleWindow, [], 'all'), 'all');
  if verbose
    fprintf(['Effective # pixels in ensemble averaging kernel: %f \n' ...
             '  axial dimension: %f \n  lateral dimension: %f \n'], ...
            effectiveKernelSize, effectiveAxialSize, effectiveLateralSize)
  end
  
  % --------------------------------------------------------------------
  % Create SNR mask
  % --------------------------------------------------------------------
  if verbose
    fprintf('Starting SNR mask creation... \n')
  end
  
  maskSNR = CreateSNRMask(tomPNCed, ensembleWindow, noiseFloorApprox);
  
  if verbose
    fprintf('  "tomPNCed"         -> "maskSNR" \n')
    fprintf('SNR mask created. \n')
  end
  
  % --------------------------------------------------------------------
  % Calculate A1C and A1CSNR
  % --------------------------------------------------------------------
  if verbose
    fprintf('Starting A1C and A1CSNR calculation... \n')
  end
  
  [alphaA1C, alphaA1CSNR] = CalculateA1CAndA1CSNR( ...
    tomPNCed, ensembleWindow, mean(noiseFloorZMCed, 2));
  clear tomPNCed
  
  if verbose
    fprintf('  "tomPNCed" cleared -> "alphaA1C" + "alphaA1CSNR" \n')
    fprintf('Calculation of A1C and A1CSNR complete! \n')
  end

  % --------------------------------------------------------------------
  % Export A1C and A1CSNR as TIFF stacks (optional)
  % --------------------------------------------------------------------
  if ~(saveFolder == "")

    if verbose
      fprintf('Exporting A1C-based and A1CSNR-based OCTA volumes as TIFF stacks... \n')
    end

    alphaA1CBscanCutSaveLoc = ...
      fullfile(saveFolder, 'alphaA1C_t' + string(tauDiff) + '_SNRMaskedBscan.tif');
    save3DMatrixAsTiffRGB(...
      maskSNR .* alphaA1C(:, :, :, tauDiff), ...
      alphaA1CBscanCutSaveLoc, [0 1], angioMap, ...
      [nZ * displayDaspect(1) / displayDaspect(2), nX, nY], 'linear')
    if verbose
      fprintf('  saved %s \n', alphaA1CBscanCutSaveLoc)
    end
    alphaA1CEnfaceCutSaveLoc = ...
      fullfile(saveFolder, 'alphaA1C_t' + string(tauDiff) + '_SNRMaskedEnface.tif');
    save3DMatrixAsTiffRGB(...
      permute(maskSNR .* alphaA1C(:, :, :, tauDiff), [3 2 1]), ...
      alphaA1CEnfaceCutSaveLoc, [0 1], angioMap, ...
      [nY * displayDaspect(2) / displayDaspect(3), nX, nZ], 'linear')
    if verbose
      fprintf('  saved %s \n', alphaA1CEnfaceCutSaveLoc)
    end
    alphaA1CSNRBscanCutSaveLoc = ...
      fullfile(saveFolder, 'alphaA1CSNR_t' + string(tauDiff) + '_SNRMaskedBscan.tif');
    save3DMatrixAsTiffRGB(...
      maskSNR .* alphaA1CSNR(:, :, :, tauDiff), ...
      alphaA1CSNRBscanCutSaveLoc, [0 1], angioMap, ...
      [nZ * displayDaspect(1) / displayDaspect(2), nX, nY], 'linear')
    if verbose
      fprintf('  saved %s \n', alphaA1CSNRBscanCutSaveLoc)
    end
    alphaA1CSNREnfaceCutSaveLoc = ...
      fullfile(saveFolder, 'alphaA1CSNR_t' + string(tauDiff) + '_SNRMaskedEnface.tif');
    save3DMatrixAsTiffRGB(...
      permute(maskSNR .* alphaA1CSNR(:, :, :, tauDiff), [3 2 1]), ...
      alphaA1CSNREnfaceCutSaveLoc, [0 1], angioMap, ...
      [nY * displayDaspect(2) / displayDaspect(3), nX, nZ], 'linear')
    if verbose
      fprintf('  saved %s \n', alphaA1CSNREnfaceCutSaveLoc)
      fprintf('A1C-based and A1CSNR-based OCTA volumes export complete! \n')
    end

  end

end

% enforce the input type that FC-OCTA will be running on
function eitherFileLocOr4DTom(input)
  if isstring(input) || ischar(input)
    return
  elseif isnumeric(input) && ndims(input) == 4
    return
  end
  error("Input must be a string/char OR a 4D numeric matrix.")
end