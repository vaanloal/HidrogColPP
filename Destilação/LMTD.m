function [deltaT] = LMTD(Thin, Thout, Tcin, Tcout)
delta1 = Thin - Tcout;
delta2 = Thout - Tcin;


deltaT = (delta1 - delta2) / log(delta1 / delta2);
end

