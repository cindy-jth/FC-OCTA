function [tomFixed] = ApplyGeneralAlinePhaseJitterCorrection(...
    tomCorrupt, zROI, noiseFloorDb, polAverage, varargin)
% APPLYGENERALALINEPHASEJITTERCORRECTION remove phase noise differences
% between A-lines
% ------------------------------------------------------------------------
% Authors:  Tianhui (Cindy) Jie, Sebastián Ruiz-Lopera, Néstor Uribe-Patarroyo
%
% TJ: 
% 1. Massachusetts Institute of Technology, 
%    Institute for Medical Engineering and Sciences, Cambridge, MA, USA;
% 2. Wellman Center for Photomedicine, Harvard Medical School, 
%    Massachusetts General Hospital, Boston, MA, USA;
% <tjie@mit.edu>
% SRL:
% 1. Massachusetts Institute of Technology, 
%    Department of Electrical Engineering and Computer Science, Cambridge, MA, USA;
% 2. Wellman Center for Photomedicine, Harvard Medical School, 
%    Massachusetts General Hospital, Boston, MA, USA;
% NUP: 
% 1. Massachusetts Institute of Technology, 
%    Institute for Medical Engineering and Sciences, Cambridge, MA, USA;
% 2. Wellman Center for Photomedicine, Harvard Medical School, 
%    Massachusetts General Hospital, Boston, MA, USA;
%
% FC-OCTA (v1.0)
%
% Changelog:
%
% 1.0 (2026-01-16): Initial version released
%
% Copyright Tianhui (Cindy) Jie, Sebastián Ruiz-Lopera, Néstor Uribe-Patarroyo (2026)
  
  % polDim is an optional argument
  if nargin > 4
    polDim = varargin{1};
  end
  
  % onlyOffset is an optional argument
  if nargin > 5
    onlyOffset = varargin{2};
  end
  
  % If noiseFloorDb is a scalar make it a vector to be indexable
  if isscalar(noiseFloorDb)
    noiseFloorDb = repmat(noiseFloorDb, [size(tomCorrupt, 1), 1]);
  end
  
  % Get info on tom
  nDims = ndims(tomCorrupt);
  % Colon operator for all dims except 1st and 2nd
  colonOp = repmat({':'}, [1 nDims - 2]);
  
  dims = size(tomCorrupt);
  % We subtract 1 from dim 2 as we only consider phase differences
  dims(2) = dims(2) - 1;
  
  % We need numel(dims) >= 3
  if numel(dims) < 3
    dims(3) = 1;
  end
  
  % We will convert tom into a very long Bscan and then reshape it back
  % once the phase is corrected
  
  nZ = size(tomCorrupt, 1);
  nZFit = numel(zROI);
  
  % We will create a system of equations like M x = B where M is a two
  % column matrix like [1, z1; 1, z1;...], x is our unkown column vector
  % [offset; slope], and B is a column vector with the measured phase
  % differences [phase1; phase2;...].
  
  % We will expand this to cover an entire Bscan with a system of equations
  % like M X = B where M is a two column matrix like [1, z1; 1, z1;...], x
  % is our unkown 2-row N-columns matrix [offset1, offset2,...; slope1,
  % slope2...] for N different Aline pairs, and B is a matrix with the measured
  % phase differences [Aline1Phase1, Aline2Phase2;... Aline1Phase2, Aline2Phase2;...].
  
  mMat = [ones(nZFit, 1), (zROI(:) - 1) / (nZ - 1)]; % z should start at 0
  
  % Conjugate product of adyacents Alines
  ccProd = tomCorrupt(zROI, 2:end, colonOp{:})...
    .* conj(tomCorrupt(zROI, 1:end - 1, colonOp{:}));
  
  % We will add weights to different z locations based on the amplitude of
  % the signal. We subtract the the noise floor from the intensity to get the weight,
  % weights have to be the same for all Alines, they can only depend on z
  % and pol channel
  weights = max(0, 10* log10(mean(abs(ccProd), 2)) - noiseFloorDb(zROI, :, colonOp{:}));
  if polAverage
    % If we are averaging by polarization channel, we need this to get a 1D
    % vector of weights
    weights = sum(weights, polDim);
  end
  
  % If we have two polarization channels, take advantage of a phase
  % difference average in that dimension if desired
  if polAverage
    ccProd = sum(ccProd, polDim);
  end
  % Now make ccProd a 2D array
  ccProd = ccProd(:, :);
  % And weights a 1D array
  weights = sum(weights(:, :), 2) / (prod(dims(3:end)) / 2); % We sum pol channels, not average them
  
  % We remove the offset first to avoid wrapping artifacts
  phaseOffset = angle(sum(ccProd, 1));
  
  if ~onlyOffset
    ccProdOffset = ccProd .* exp(-1i .* phaseOffset);
    
    % Now create the system of equations with M -> mMat, x -> xMat, B ->
    % bMat.
    bMat = angle(ccProdOffset);
    
    % diagonalize weights
    diagWeights = diag(weights);
    
    % Now we solve in this way (T means transpose):
    %                        M x = B
    % We now add weights in W
    %                      W M x = W B
    %                 (MT W M) x = MT W B
    %   (MT W M)^{-1} (MT W M) x = (MT W M)^{-1} MT W B
    %                          x = (MT W M)^{-1} MT W B
    
    % find the weighted LSF
    xMat =  (mMat.' * diagWeights * mMat) \ mMat.' * diagWeights * bMat;
    
    % Next step, correct the phases based on this slope, this time
    % considering all z
    zVecAll = (0:nZ - 1).' / (nZ - 1);
    
    % Now get offset and add the original phase offset we subtracted
    offsetDiff = xMat(1,:) + phaseOffset;
    % Regain dimensions, we set polDim to 1 to not waste memory
    if polAverage
      dims(polDim) = 1;
    end
    offsetDiff = reshape(offsetDiff, [1 dims(2:end)]);
    % Add zeros at the beginning to set the correction to each first Aline as
    % nothing
    offsetDiff = padarray(offsetDiff, [0, 1], 0, 'pre');
    % Do cumsum to obtain cumulative correction only across 2nd dimension
    offset = mod(cumsum(offsetDiff, 2), 2 * pi);
    
    % Regain dimensions
    slopeDiff = reshape(xMat(2, :), [1 dims(2:end)]);
    % Add zeros at the beginning to set the correction to each first Aline as
    % nothing
    slopeDiff = padarray(slopeDiff, [0, 1], 0, 'pre');
    % Do cumsum to obtain cumulative correction only across 2nd dimension
    slope = cumsum(slopeDiff, 2);
    
    % Now calculate correction with offset and slope
    corrPhase = slope .* zVecAll + offset;
    
    % Apply correction
    tomFixed = exp(-1i * corrPhase) .* tomCorrupt;
  else
    offsetDiff = phaseOffset;
    % Regain dimensions, we set polDim to 1 to not waste memory
    if polAverage
      dims(polDim) = 1;
    end
    offsetDiff = reshape(offsetDiff, [1 dims(2:end)]);
    % Add zeros at the beginning to set the correction to each first Aline as
    % nothing
    offsetDiff = padarray(offsetDiff, [0, 1], 0, 'pre');
    % Do cumsum to obtain cumulative correction only across 2nd dimension
    offset = mod(cumsum(offsetDiff, 2), 2 * pi);
    
    % Now calculate correction with offset and slope
    corrPhase = offset;
    
    % Apply correction
    tomFixed = exp(-1i * corrPhase) .* tomCorrupt;
  end
end




