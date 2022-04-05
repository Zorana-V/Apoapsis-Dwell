function TI2E = dcmI2E_Visic_Zorana(t,OmegaE)
% --------------------------------------------------------------%
% -------------------------- DCMI2E.M --------------------------%
% --------------------------------------------------------------%
% Compute the transformation matrix from Earth-centered inertial%
% Cartesian coordinates to Earth-centered Earth-fixed Cartesian %
% coordinates. The transformation matrix is of size 3 by 3 and  %
% represents the transformation matrix at the time t            %
% Inputs:                                                       %
% t: the time at which the transformation is desired            %
% OmegaE: the rotation rate of the Earth                        %
% Outputs:                                                      %
% TI2E: 3 by 3 matrix that transforms components of a           %
% vector expressed in Earth-centered inertial                   %
% Cartesian coordinates to Earth-centered                       %
% Earth-fixed Cartesian coordinates                             %
% --------------------------------------------------------------%
% NOTE: the units of the inputs T and OMEGAE must be consistent %
% and the angular rate must be in RADIANS PER TIME UNIT         %
% --------------------------------------------------------------%

%first we built the transformation matrix using the given information
TI2E = [ cos(OmegaE*t) sin(OmegaE*t) 0; (-sin(OmegaE*t)) cos(OmegaE*t) 0; 0 0 1]; 
%transformation matrix
end
