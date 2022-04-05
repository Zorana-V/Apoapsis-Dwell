function rvValuesECEF = eci2ecef_Visic_Zorana(tValues,rvValuesECI,OmegaE)
% --------------------------------------------------------------%
% ------------------------- ECI2ECEF.M -------------------------%
% --------------------------------------------------------------%
% Transform the values of the position along a spacecraft orbit %
% from Earth-centered inertial Cartesian coordinates to         %
% Earth-centered Earth-fixed Cartesian coordinates.             %
% represents the transformation matrix at the time t            %
% Inputs:                                                       %
% tValues: column matrix of values of time at which the         %
% transformation is to be performed                             %
% rvValuesECI: matrix of values of the position expressed       %
% in Earth-centered inertial Cartesian                          %
% coordinates and stored ROW-WISE                               %
% OmegaE: the rotation rate of the Earth                        %
% Outputs:                                                      %
% rvValuesECEF: matrix of values of the position expressed      %
% in Earth-centered Earth-fixed Cartesian                       %
% and stored ROW-WISE                                           %
% --------------------------------------------------------------%


%for loop to include each value of t and its corresponding position vector
for ii = 1:length(tValues)
 %creating the transformation matrix from I to E
 TI2E = dcmI2E_Visic_Zorana(tValues(ii),OmegaE);
 %vector multiplication to convert from ECI to ECEF
 vec = TI2E*(rvValuesECI(ii,:).');
 
 for jj = 1:3
 %transform back into row vector
 rvValuesECEF(ii,jj) = vec(jj);
end
end
