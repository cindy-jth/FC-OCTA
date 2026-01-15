function [ signalHighPass ] = RealHighPassFilter(signal, filterHalfSize, ind, varargin)
  %RealHighPassFilter Summary of this function goes here
  %   Detailed explanation goes here
  
  if nargin > 3 && ~isempty(varargin{1})
    windowType = varargin{1};
  else
    windowType = 'hat';
  end
  
  % filterHalfSize is size IN EXCESS of the DC component. So filterHalfSize = 0
  % filters out DC.
  switch windowType
    case 'hat'
      window = ones(2 * filterHalfSize + 1, 1);
    case 'hanning'
      window = hanning(2 * filterHalfSize + 1);
    case 'blackman'
      window = blackman(2 * filterHalfSize + 1);
    otherwise
      error('unknown window type')
  end
  
  nSamples = size(signal, ind);
  window = padarray(1 - window, floor(nSamples / 2) - filterHalfSize, 1, 'pre');
  window = padarray(window, ceil(nSamples / 2) - filterHalfSize - 1, 1, 'post');
  window = shiftdim(window(:), 1 - ind);
  
  signalFT = fft(fftshift(signal, ind), [], ind);
  signalFT = signalFT .* fftshift(window, ind);
  signalHighPass = real(ifftshift(ifft(signalFT, [], ind, 'symmetric'), ind));
  
end

