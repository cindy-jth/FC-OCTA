function mask = CreateSNRMask(tom, ensembleWindow, noiseFloorApprox)
arguments (Input)
  tom (:, :, :, :, :) 
  ensembleWindow (:,:,:) 
  noiseFloorApprox (1, 1)
end

arguments (Output)
  mask
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
