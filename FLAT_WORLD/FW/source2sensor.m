function GSS=source2sensors(S,SS)
%GSS=source2sensors(S,SS) - gain from sources to sensors
%
%   GSS: gain matrix (nsources X nsensors)
%
%   S: source positions (nsources X 2)
%   SS: sensor positions (nsensors X2)

GSS=1./sqdist(S,SS);

