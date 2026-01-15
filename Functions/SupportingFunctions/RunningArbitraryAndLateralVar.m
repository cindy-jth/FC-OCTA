function [ variance ] = RunningArbitraryAndLateralVar(x, ensembleWindow, norm, dim, varargin)
  %UNTITLED Summary of this function goes here
  %   Detailed explanation goes here
  
  if nargin > 4 && ~isempty(varargin{1})
    keepType = varargin{1};
  else
    keepType = false;
  end
  
  if keepType
    typeOrig = class(x);
  end
  x = double(x);
  % Avoid catastrophic cancellation
  variance = RunningArbitraryAndLateralAve(abs(bsxfun(@plus, x, -mean(x, dim))) .^ 2, ensembleWindow, dim)...
    - abs(RunningArbitraryAndLateralAve(bsxfun(@plus, x, -mean(x, dim)), ensembleWindow, dim)) .^ 2;
  if norm == 0
    variance = variance * (ensembleWindow * xCorrWindow * size(x, dim)) / (sum(ensembleWindow(:)) * size(x, dim) - 1);
  end
  
  if keepType
    variance = cast(variance, typeOrig);
  end
end

