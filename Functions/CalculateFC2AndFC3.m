function [alphaFC2, alphaFC3] = CalculateFC2AndFC3(tom, ensembleWindow, noiseFloorZ, tauVec)
arguments (Input)
  tom (:, :, :, :, :) 
  ensembleWindow (:,:,:) 
  noiseFloorZ (:, 1, 1, :, :)
  tauVec (1, :) = 1 : (size(tom, 2) - 1)
end

arguments (Output)
  alphaFC2
  alphaFC3
end

  [nZ, ~, nX, nY, ~] = size(tom);
  nTaus = numel(tauVec);

  [fracG1, ~, fracG1SNR] = CalcArbitraryEnsembleAveragedG1(tom, ...
    tauVec, ensembleWindow, noiseFloorZ, false, 5);
  
  alphaFC2 = max(0, real(sqrt(1 - abs(fracG1))));
  alphaFC2 = permute(reshape(alphaFC2, [nZ, nTaus, nX, nY]), [1 3 4 2]);

  alphaFC3 = max(0, real(sqrt(1 - abs(fracG1SNR))));
  alphaFC3 = permute(reshape(alphaFC3, [nZ, nTaus, nX, nY]), [1 3 4 2]);

end
