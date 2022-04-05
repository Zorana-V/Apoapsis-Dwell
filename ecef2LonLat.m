function [lonE,lat] = ecef2LonLat(rvValuesECEF)
% --------------------------------------------------------------%
% ----------------------- ECEF2LONLAT.M ------------------------%
% --------------------------------------------------------------%
% Compute the Earth-relative longitude and geocentric latitude  %
% from the Earth-centered Earth-fixed position.                 %
% Input:                                                        %
% rvValuesECEF: matrix of values of the position expressed      %
% in Earth-centered Earth-fixed Cartesian                       %
% and stored ROW-WISE                                           %
% Outputs:                                                      %
% lonE: a column matrix containing the values of the            %
% Earth-relative longitude (radians)                            %
% lat: a column matrix containing the values of the             %
% geocentric latitude (radians)                                 %
% --------------------------------------------------------------%

%for loop to consider every row of the vector
 for ii = 1:size(rvValuesECEF,1)
 %deconstructing the rvValuesECEF vector into its 3D components
 xe = rvValuesECEF(ii,1);
 ye = rvValuesECEF(ii,2);
 ze = rvValuesECEF(ii,3);
 
 %solving for longitude and latitude
 lonE(ii) = atan2(ye,xe);
 lat(ii) = atan2(ze,sqrt((ye^2)+(xe^2)));
 end
end
