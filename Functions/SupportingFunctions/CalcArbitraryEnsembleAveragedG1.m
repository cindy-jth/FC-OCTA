function [corrG, varargout] = CalcArbitraryEnsembleAveragedG1(tom, tauVec,...
    ensembleWindow, noiseFloor, varargin)
  %CalcEnsembleAveragedG1 First-order autocorrelation function
  %   Calculates the autocorrelation along 2nd index, with optional
  %   inputs.
  %
  % Inputs:
  %   tomInt:     array with the signal
  %   tauVec:     vector with desired delays
  %   ensembleWindow: kernel for ensemble averaging in arbitrary dimensions.
  %   noiseFloor:      noise floor in linear scale, required to calculate
  %                    SNR-corrected ACs
  %   isSNR:           if SNR is supplied instead of noiseFloor
  %   addEnsembleDims: vector with list of dimensions for further ensemble averaging (dims are collapsed)
  
  %
  % Outputs:
  %   corrG: 1st order AC
  %   meanInt: mean intensity inside each correlation window
  %   corrGSNR: SNR-corrected 1st order AC
  %   
  % This script and its functions follow the coding style that can be
  % sumarized in:
  % * Variables have lower camel case
  % * Functions upper camel case
  % * Constants all upper case
  % * Spaces around operators
  %
  % Authors:  Néstor Uribe-Patarroyo
  %
  % NUP: 
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA;
  % <uribepatarroyo.nestor@mgh.harvard.edu>

  % MGH Flow Measurement project (v1.0)
  %
  % Changelog:
  %
  % V1.0 (2017-09-01): Initial version released
  %
  % Copyright Néstor Uribe-Patarroyo (2017)
  
  if nargin >= 6 && ~isempty(varargin{2})
    addEnsembleDims = varargin{2};
  else
    addEnsembleDims = [];
  end
  
  % Here I replace all the ensemble averages to include axial running averaging
  % too!
  if numel(ensembleWindow) > 1
    % Shift all dims >= 2 by two due to the similar shift done to the signal
    % below
    if ~iscolumn(ensembleWindow)
      ensembleWindow = shiftdim(ensembleWindow, -1);
      ensembleWindow = permute(ensembleWindow, [2, 1, 3:ndims(ensembleWindow)]);
    end
    meanEnsemble = @(x, dim) RunningArbitraryAndLateralAve(x, ensembleWindow, [dim, addEnsembleDims]);
    varEnsemble = @(x, norm, dim) RunningArbitraryAndLateralVar(x, ensembleWindow, norm, [dim, addEnsembleDims]);
  else
    meanEnsemble = @(x, dim) mean(x, [dim, addEnsembleDims]);
    varEnsemble = @(x, norm, dim) var(x, norm, [dim, addEnsembleDims]);
  end
  
  % Get info on tom
  nDims = ndims(tom);
  colonOp = repmat({':'}, [1 nDims - 2]);
  nZ = size(tom, 1);
  corrWindow = size(tom, 2);
  nOther = num2cell(size(tom));
  nOther(1:2) = [];
  % And make addEnsembleDims singleton
  addEnsembleDimsCell = num2cell(addEnsembleDims - 2);
  nOther([addEnsembleDimsCell{:}]) = deal({1});
  % Possibly fix tauVec size
  tauVec = unique(min(tauVec, corrWindow - 1));
  nTaus = numel(tauVec);
  
  % If noise floor is empty, then don't calculate SNR-corrected version
  if isempty(noiseFloor) || nargin < 2 || nargout <= 2
    calcSNRCorrected = false;
  else
    calcSNRCorrected = true;
  end
  
  % Now we create the complex-conjugate products maxtrices
  corrG = zeros(nZ, nTaus, nOther{:}, 'like', tom);
  if calcSNRCorrected
    corrGSNRCorrected = zeros(nZ, nTaus, nOther{:}, 'like', tom);
  end
  if nargout > 1
    calcMeanInt = true;
    meanInt = zeros(nZ, 1, nOther{:}, 'like', tom);
    if nargout > 3
      calcVariances = true;
      ccProdGVar = zeros(nZ, nTaus, nOther{:}, 'like', tom);
      ccProdGSNRCorrectedVar = zeros(nZ, nTaus, nOther{:}, 'like', tom);
    else
      calcVariances = false;
    end
  else
    calcMeanInt = false;
    calcVariances = false;
  end
  
  % See if we are supplying SNR directly and not noise floor
  if nargin > 4 && ~isempty(varargin{1})
    isSNR = varargin{1};
  else
    isSNR = false;
  end
  
  l = 0;
  for thisTau = tauVec
    l = l + 1;
    idx1 = 1:corrWindow - thisTau;
    idx2 = 1 + thisTau:corrWindow;
    
    % Get signals
    g1 = tom(:, idx1, colonOp{:});
    g2 = tom(:, idx2, colonOp{:});
    
    % Calculate complex-conjugate products with proper normalization.
    % Ensemble averages are in the 2nd index!
    ccProdGDenominator = sqrt(meanEnsemble(abs(g1) .^ 2, 2) .* meanEnsemble(abs(g2) .^ 2, 2));
    % Get meanEnsemble of ccProd to obtain correlation coefficient of ensemble
    corrG(:, l, colonOp{:}) = (meanEnsemble(conj(g1) .* g2, 2) ./ ccProdGDenominator);
    
    
    if calcSNRCorrected
      if isSNR
        snr1 = noiseFloor;
        snr2 = noiseFloor;
      else
        % Estimate SNR
        if 0
          % In practice we would estimate the noise floor independently, and
          % estimating a different SNR for each part of the signal
          snr1 = max(0, (meanEnsemble(abs(g1) .^ 2, 2) ./ noiseFloor) - 1);
          snr2 = max(0, (meanEnsemble(abs(g2) .^ 2, 2) ./ noiseFloor) - 1);
        else
          % In practice we would estimate the noise floor independently, and
          % using the same SNR for each part of the signal
          thisSignalInt = meanEnsemble(abs(tom) .^ 2, 2);
          % Yes, same estimation for all l's
          if ~isempty(addEnsembleDims)
            % If noiseFloor is non-singleton along addEnsembleDims, average
            % there too
            snr1 = max(0, (thisSignalInt ./ mean(noiseFloor, addEnsembleDims)) - 1);
          else
            snr1 = max(0, (thisSignalInt ./ noiseFloor) - 1);
          end
          snr2 = snr1;
        end
      end
      
      % From SNR, correct value of cc products
      corrFactor = sqrt((1 + 1 ./ snr1) .* (1 + 1 ./ snr2));
      % And get SNR-corrected value of correlation coeficient
      corrGSNRCorrected(:, l, colonOp{:}) = corrG(:, l, colonOp{:}) .* corrFactor;
    end
    
    if calcVariances
      % Calculate the normalized variance (variance divided by the meanEnsemble)
      ccProdGVar(:, l, colonOp{:}) = bsxfun(@rdivide, varEnsemble(conj(g1) .* g2, 2), ccProdGDenominator);
      if calcSNRCorrected
        ccProdGSNRCorrectedVar(:, l, colonOp{:}) = bsxfun(@times, bsxfun(@rdivide, varEnsemble(conj(g1) .* g2, 2), ccProdGDenominator), corrFactor);
      end
    end
    
  end
  if calcSNRCorrected
    corrGSNRCorrected(:, tauVec == 0, colonOp{:}) = 1; % Because after the correction we know we need this
  end
  
  % Now I know how to calculate the var (the denominator is the same!). But what's
  % the meaning of var and can it give me something interesting?
  
  if calcMeanInt
    meanInt(:, 1, colonOp{:}) = meanEnsemble(abs(tom) .^ 2, 2);
    varargout{1} = meanInt;
  end
  
  if calcSNRCorrected
    varargout{2} = corrGSNRCorrected;
  end
  if calcVariances
    varargout{3} = ccProdGVar;
    if calcSNRCorrected
      varargout{4} = ccProdGSNRCorrectedVar;
    end
  end
  
end

