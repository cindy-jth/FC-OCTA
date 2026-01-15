function tomShifted = ApplyAxialShift(tom, axialShift)
  
  nZ = size(tom, 1);
  
  % Create vectors to apply shift
  zVect = ifftshift(-fix(nZ / 2):ceil(nZ / 2) - 1).';
  
  % Calculate FT of complex tomogram in Z
  tomFT = fft(tom, [], 1);
  
  % Apply shift in Fourier domain as a linear phase term
  tomShifted = ifft(tomFT .* exp(-2i * pi * (axialShift .* zVect / nZ)), [], 1);
end