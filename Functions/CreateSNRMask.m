function mask = CreateSNRMask(tom, ensembleWindow, noiseFloorApprox)
% CREATESNRMASK calculate 3D SNR mask for OCTA volumes based on the OCT
% tomogram signal strength
% ------------------------------------------------------------------------
% Authors:  Tianhui (Cindy) Jie
%
% TJ: 
% 1. Massachusetts Institute of Technology, 
%    Institute for Medical Engineering and Sciences, Boston, MA, USA
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
  tom (:, :, :, :, :)      % OCT tomogram [nZ, nReps, nX, nY, nPolChs]
  ensembleWindow (:,:,:)   % averaging kernel [zWindowSize, xWindowSize, yWindowSize]
  noiseFloorApprox (1, 1)  % average noise floor (intensity scale) reference
end

arguments (Output)
  mask  % SNR mask: 3D matrix [nZ, nX, nY] (binary mask with sigmoid edge transition)
end

  [nZ, ~, nX, nY, ~] = size(tom);

  tomInt = RunningArbitraryAndLateralAve(abs(tom) .^ 2, ensembleWindow, [2 5]);
  tomInt = reshape(tomInt, [nZ nX nY]);
  
  % Create a mask to identify any SNR < -3 dB for which we will discard
  % data using sigmoid function
  dBSNR = 10 * log10(max(0, tomInt / noiseFloorApprox - 1));
  mask = 1 ./ (1 + exp(-2 .* (dBSNR - 1)));
  mask(mask < 0.3) = 0;

end
