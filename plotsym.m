function H = plotsym(X,Y,Z,symbol)
% PLOTSYM generates a plot using symbols.
%
% PLOTSYM(X,Y) behaves exactly as plot(x,y).
%
% PLOTSYM(X,SYMBOL) produces a plot using the symbol
% defined by SYMBOL.  If SYMBOL is a single element vector,
% then is assumed to be the ASCII value of the symbol, and
% it will be converted to a string using setstr.
%
% PLOTSYM(X,Y,SYMBOL) produces a plot using the symbol
% defined by SYMBOL.  If SYMBOL is a single element vector,
% then is assumed to be the ASCII value of the symbol, and
% it will be converted to a string using setstr.
%
% PLOTSYM(X,Y,Z) behaves exactly as plot3(x,y,z) except if
% Z is a scalar.  When Z is a scalar, it is assumed to be the
% ASCII value of the symbol, and it will be converted to a
% string using setstr.
%
% PLOTSYM(X,Y,Z,SYMBOL) produces a 3-D plot using the symbol
% defined in SYMBOL.
%
% H = PLOTSYM3(...) returns the handle(s) to the plot.
%
% PLOTSYM uses the text command to produce the plot when SYMBOL
% is defined.
%
% EXAMPLES:
%
%     x = 0:pi/10:2*pi;
%     y = sin(x);
%     z = x;
%     plotsym(x,y,'A')
%     plotsym(x,y,z,154)  % Some random ASCII value 154

% Written by John L. Galenski III - May 26, 1994
% Copyright (c) by The MathWorks, Inc.

% Parse the input
if nargin < 1                              % Nothing given
  error('Requires at least one input.')
elseif (nargin == 1)                       % X given
  if (~isstr(X)) & (length(X) > 1)
    h = plot(X);
    return;
  else
    error('Input must be a vector.')
  end
elseif nargin == 2                         % X and Y or SYMBOL
  if isstr(Y) | length(Y)==1               % given
    symbol = Y;
    Y = X;
    X = 1:length(Y);
        Z = zeros(size(X));
    if ~isstr(symbol), symbol = setstr(symbol); end
  else
    h = plot(X,Y);
        return;
  end
elseif nargin == 3                         % X, Y, and Z or SYMBOL
  if isstr(Z) | length(Z)==1               % given
    if isstr(Z), symbol = Z;
    else, symbol = setstr(Z);
    end
    Z = zeros(size(X));
  else
    h = plot3(X,Y,Z);
    return;
  end
elseif nargin == 4                         % Everything given
  if ~isstr(symbol), symbol = setstr(symbol); end
else
  error('Invalid number of inputs.')
end

% Check to see if HOLD is on
if ~ishold
  newplot;
  % Set the correct limits.
  tmp = plot3(X,Y,Z);
  if all(Z == 0)
    view(2);
    set(gca,'Box','on')
  end
  lim = axis;
  delete(tmp);
  axis(lim);
end

% Generate the plot
h = text(X',Y',Z',symbol,'VerticalAlignment','middle', ...
    'HorizontalAlignment','center');

% Output?
if nargout
  H = h;
end