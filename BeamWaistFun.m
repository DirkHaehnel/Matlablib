function err = BeamWaistFun(w0,w1,w2,d,lambda);

err = (abs(w0*(sqrt(w2^2-w0^2) - sqrt(w1^2-w0^2))) - lambda*d/pi)^2

% beam waist position w0*sqrt(w1^2-w0^2)*pi/lambda