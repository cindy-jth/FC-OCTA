function [alphaA1C, alphaA1CSNR] = CalculateA1CAndA1CSNR(tom, ensembleWindow, noiseFloorZ, tauVec)
% CALCULATEA1CANDA1CSNR calculate OCTA volumes using A1C and A1CSNR metrics
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
  tom (:, :, :, :, :)                     % OCT tomogram [nZ, nReps, nX, nY, nPolChs]
  ensembleWindow (:,:,:)                  % averaging kernel [zWindowSize, xWindowSize, yWindowSize]
  noiseFloorZ (:, 1, 1, :, :)             % for each B-scan, average noise floor (intensity scale) at each depth [nZ, 1, 1, nY, nPolChs]
  tauVec (1, :) = 1 : (size(tom, 2) - 1)  % (default all: nReps-1) combos of tauDiff (B-scan step between repetitions used to calculate OCTA)
end

arguments (Output)
  alphaA1C  % OCTA volume based on A1C: 4D matrix [nZ (depth), nX (slow axis), nY (fast axis), nReps-1 (all combos of tauDiff)]
  alphaA1CSNR  % OCTA volume based on A1CSNR: 4D matrix [nZ (depth), nX (slow axis), nY (fast axis), nReps-1 (all combos of tauDiff)]
end

  [nZ, ~, nX, nY, ~] = size(tom);
  nTaus = numel(tauVec);

  [fracG1, ~, fracG1SNR] = CalcArbitraryEnsembleAveragedG1(tom, ...
    tauVec, ensembleWindow, noiseFloorZ, false, 5);
  
  alphaA1C = max(0, real(sqrt(1 - abs(fracG1))));
  alphaA1C = permute(reshape(alphaA1C, [nZ, nTaus, nX, nY]), [1 3 4 2]);

  alphaA1CSNR = max(0, real(sqrt(1 - abs(fracG1SNR))));
  alphaA1CSNR = permute(reshape(alphaA1CSNR, [nZ, nTaus, nX, nY]), [1 3 4 2]);

end
