function [ average ] = RunningArbitraryAndLateralAve(array, ensembleWindow, dim)
  %UNTITLED Summary of this function goes here
  %   Detailed explanation goes here
  
  % Because we moved all dim > 2 by one or two, no need to change any dims here.
  
  % Just use imfilter for coherent 2D averaging along 1st and all other
  % spatial dimensions (indices >= 4)
  sumAx = imfilter(array, ensembleWindow, 'replicate');
  % Also coherent averaging along dim (can be a vector with multiple dims)
  sumLat = sum(sumAx, dim);
  % Count number of elements to do averaging. prod is required in case dim
  % is a vector
  average = sumLat / (prod(size(array, dim)) * sum(ensembleWindow(:)));
end

