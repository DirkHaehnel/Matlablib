% getopt -- parse command options.
%
%      [opt, arg] = getopt(spec, ...)
%
% The 'getopt' command parses options from an argument list.
%
% First argument SPEC is a structure array where each element describes
% an option.  An option structure has these fields:
%
% name
%      The name of the option (a string) preceded by a hyphen.
% type
%      The type of the option (a string).  This field determines whether
%      the option requires an argument or not.  The option type is either
%      one of the built-in option types (see below), or a Matlab class
%      name.
% index
%      The structure array index of the option.  This field is only
%      mandatory for alias options.
% choices
% values
%      The choices field contains a set of valid option arguments.  If the
%      choices field is empty, any argument of the correct type is accepted.
%      If the choices *and* values field is not empty, the option argument
%      is replaced by the corresponding choice value.
% value
%      The value of the option.   For options requiring an argument, the
%      option value is set to the option argument.
%
% The 'getopt' command returns the modified options structure array OPT,
% and the remaining arguments ARG.  Options and arguments are permuted while
% being parsed.  The special argument '--' terminates parsing in all cases.
%
% An option argument can be provided in two ways.  Either as a separate
% argument immediately following the option it belongs to, or together with
% the option separated from it by an equal sign.
%
% If the argument SPEC is empty, 'getopt' does not accept any options and
% the return value OPT is an empty option structure array.
%
% Built-in option types:
%
% alias
%      The option is an alternative name for another option.  The index
%      field is the structure array index of the original option.  If the
%      value field is not empty, an alias option does not take an argument
%      but sets the value of the original option to that value.
% color
%      The option requires a Matlab color specification as its argument.
%      Option argument can be a short or long Matlab color name, a RGB
%      intensity tripple, or a gray-scale intensity value.
% cursor
%      The option requires a Matlab cursor name as its argument.
% line
%      The option requires a Matlab line style as its argument.
%      An argument of 'solid', 'dashed', 'dotted', and 'dashdot' is
%      accepted as an alternative name for the corresponding Matlab
%      line style.
% marker
%      The option requires a Matlab marker symbol as its argument.
%      An argument of 'plus', 'circle', 'asterisk', 'star', 'point',
%      'cross', 'x-mark', 'square', 'diamond', 'up', 'down', 'left',
%      'right', 'pentagram', and 'hexagram' is accepted as an alternative
%      name for the corresponding Matlab marker symbol.
% flag
%      The option is a simple "has seen" option.  Each time the option
%      occurs the value is incremented by one.  The option value defaults
%      to zero.
% symbol
%      The option argument must be a valid Matlab symbol name.
% toggle
%      The option is an on/off option.  Each time the option occurs, the
%      value is logically negated.  The option value defaults to zero.
%
% Example:
%
%      function foo(varargin)
%
%        % Initialize options.
%        for ind = 1:inf
%          switch ind
%           case 1
%            opt(ind).name = '-flag';
%           case 2
%            opt(ind).name = '-mode';
%            opt(ind).type = 'char';
%            opt(ind).choices = {'box', 'line'};
%            opt(ind).values = [4, 2];
%            opt(ind).value = 4;
%           case 3
%            opt(ind).name = '-color';
%            opt(ind).type = 'color';
%           case 4
%            opt(ind).name = '-cursor';
%            opt(ind).type = 'cursor';
%           case 5
%            opt(ind).name = '-fullcross';
%            opt(ind).type = 'alias';
%            opt(ind).index = 4;
%            opt(ind).value = 'fullcrosshair';
%           otherwise
%            break;
%          end
%        end
%
%        % Parse options from argument list.
%        [opt, arg] = getopt(opt, varargin{:});

%% getopt.m --- the preceding comment is the documentation string.

% Copyright (C) 2004 Ralph Schleicher

% Author: Ralph Schleicher <rs@nunatak.allgaeu.org>
% Time-stamp: <2004-07-15 22:36:45 CEST>

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation; either version 2,
% or (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program; see the file COPYING.  If not, write to
% the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA 02111-1307, USA.

% As a special exception, Ralph Schleicher gives permission to link
% the code of this program with MATLAB from The Mathworks, Inc. (or
% with modified versions of MATLAB that use the same license as
% MATLAB), and distribute linked combinations including the two.
% You must obey the GNU General Public License in all respects for
% all of the code used other than with MATLAB.  If you modify this
% file, you may extend this exception to your version of the file,
% but you are not obligated to do so.  If you do not wish to do so,
% delete this exception statement from your version.

%% Commentary:

% Uncomment the following line to run a local test.
% function local_test, getopt_test

%% Code:


% Program entry point.
function [opt, rest] = getopt(spec, varargin)

  % Check number of arguments.
  error(nargchk(1, inf, nargin));

  % Create empty structure array.
  if isempty(spec)
    spec = repmat(struct('name', [], ...
                         'type', [], ...
                         'index', [], ...
                         'choices', [], ...
                         'values', [], ...
                         'value', []), ...
                  0, 0);
  end

  % Check options structure.
  opt = spec;

  % Add missing fields.
  fields = fieldnames(opt);
  switch 'name'
   case fields
   otherwise
    error('Invalid options structure');
  end
  switch 'type'
   case fields
   otherwise
    [opt.type] = deal([]);
  end
  switch 'index'
   case fields
   otherwise
    [opt.index] = deal([]);
  end
  switch 'choices'
   case fields
   otherwise
    [opt.choices] = deal([]);
  end
  switch 'values'
   case fields
   otherwise
    [opt.values] = deal([]);
  end
  switch 'value'
   case fields
   otherwise
    [opt.value] = deal([]);
  end

  % Check field values.
  for ind = 1:numel(opt)
    if ~ getopt_option_p(opt(ind).name)
      error(sprintf('Option %d has an invalid option name', ind));
    end
    if isempty(opt(ind).type)
      opt(ind).type = 'flag';
    end
    switch opt(ind).type
     case 'alias'
      if isempty(opt(ind).index)
        error(sprintf('Option %d requires an option index', ind));
      end
     case 'flag'
      opt(ind).arg = 0;
      if isempty(opt(ind).value)
        opt(ind).value = 0;
      end
     case 'toggle'
      opt(ind).arg = 0;
      if isempty(opt(ind).value)
        opt(ind).value = 0;
      end
     otherwise
      opt(ind).arg = 1;
    end
    if ~ isempty(opt(ind).index)
      if all(opt(ind).index ~= 1:numel(opt))
        error(sprintf('Option %d has an invalid option index', ind));
      end
    end
    if ~ isempty(opt(ind).choices) & ~ isempty(opt(ind).values)
      if numel(opt(ind).choices) ~= numel(opt(ind).values)
        error(sprintf('Option %d has a non-matching number of choice/value pairs', ind));
      end
    end
  end

  % Check alias loops.
  for ind = 1:numel(opt)
    switch opt(ind).type
     case 'alias'
      next = ind;
      while strcmp(opt(next).type, 'alias')
        next = opt(next).index;
        if next == ind
          error(sprintf('Option %d refers to itself', ind));
        end
      end
    end
  end

  % Known option names.
  options = {opt.name};

  % Index of the next element of the argument list to be processed.
  optind = 1;

  % While scanning the argument list, options and option arguments
  % are removed from the list.  At the end, the list contains all
  % the remaining arguments.
  while optind <= numel(varargin)

    % Next element of the argument list.
    optopt = varargin{optind};

    % The special argument '--' terminates parsing in all cases.
    if ischar(optopt) & strcmp(optopt, '--')
      varargin(optind) = [];
      break;
    end

    % Skip over non-option arguments.
    if ~ getopt_option_p(optopt)
      optind = optind + 1;
      continue;
    end

    % Remove option from argument list.
    varargin(optind) = [];

    % Check for embedded option argument.
    has_arg = strfind(optopt, '=');
    if has_arg
      start = has_arg(1);
      % Extract option argument (may be empty, take care).
      optarg = optopt(start:end);
      optarg(1) = [];
      % Evaluate simple (constant) arguments.
      optarg = getopt_eval(optarg);
      % Delete option argument from option.
      optopt(start:end) = [];
    end

    % Only process known options.
    switch optopt
     case options
     otherwise
      error(sprintf('Unrecognized option ''%s''', optopt));
    end

    % Determine option structure index (can't fail).
    for ind = 1:numel(opt)
      if strcmp(optopt, opt(ind).name)
        break;
      end
    end

    % Resolve aliases.
    switch opt(ind).type
     case 'alias'
      optval = opt(ind).value;
      while strcmp(opt(ind).type, 'alias')
        ind = opt(ind).index;
        if isempty(optval)
          optval = opt(ind).value;
        end
      end
     otherwise
      optval = [];
    end

    % Update the option value.
    if opt(ind).arg == 0
      % Option has no argument.
      if has_arg
        error(sprintf('Option ''%s'' does not take an argument', optopt));
      end
      if isempty(optval)
        switch opt(ind).type
         case 'flag'
          optval = opt(ind).value + 1;
         case 'toggle'
          optval = ~ opt(ind).value;
         otherwise
          error(sprintf('Option %d has an unknown option type', ind));
        end
      end
    else
      % Option has an argument.  The argument is either provided implicitly
      % by the option value of an alias option or explicitly by an option
      % argument.
      if has_arg
        if ~ isempty(optval)
          error(sprintf('Option ''%s'' does not take an argument', optopt));
        end
      else
        if ~ isempty(optval)
          optarg = optval;
        else
          if numel(varargin) < optind
            error(sprintf('Option ''%s'' requires an argument', optopt));
          end
          % Get option argument from the argument list.
          optarg = varargin{optind};
          varargin(optind) = [];
          if getopt_option_p(optarg)
            error(sprintf('Option ''%s'' requires an argument', optopt));
          end
        end
      end
      switch opt(ind).type
       case 'color'
        [true, optval] = getopt_color_p(optarg);
        if ~ true
          error(sprintf('Option ''%s'' requires a Matlab color specification', optopt));
        end
       case 'cursor'
        [true, optval] = getopt_cursor_p(optarg)
        if ~ true
          error(sprintf('Option ''%s'' requires a Matlab cursor name', optopt));
        end
       case 'line'
        [true, optval] = getopt_line_p(optarg);
        if ~ true
          error(sprintf('Option ''%s'' requires a Matlab line style', optopt));
        end
       case 'marker'
        [true, optval] = getopt_marker_p(optarg)
        if ~ true
          error(sprintf('Option ''%s'' requires a Matlab marker symbol', optopt));
        end
       case 'symbol'
        optval = optarg;
        if ~ isvarname(optarg)
          error(sprintf('Option ''%s'' requires a Matlab symbol name', optopt));
        end
       otherwise
        optval = optarg;
        if ~ isa(optarg, opt(ind).type)
          error(sprintf('Option ''%s'' requires an argument of class ''%s''', optopt, opt(ind).type));
        end
        if ~ isempty(opt(ind).choices)
          tem = 0;
          for i = 1:numel(opt(ind).choices)
            if iscell(opt(ind).choices)
              elem = opt(ind).choices{i};
            else
              elem = opt(ind).choices(i);
            end
            if isequal(optarg, elem)
              tem = i;
              break;
            end
          end
          if tem == 0
            error(sprintf('Invalid argument for option ''%s''', optopt));
          end
          if ~ isempty(opt(ind).values)
            if iscell(class(opt(ind).values))
              optval = opt(ind).values{i};
            else
              optval = opt(ind).values(i);
            end
          end
        end
      end
    end

    % Assign option value.
    opt(ind).value = optval;
  end

  % Remove undocumented option structure fields.
  opt = rmfield(opt, 'arg');

  % Remaining arguments.
  rest = varargin;


% Attempt to evaluate argument ARG.
function val = getopt_eval(arg)

  if exist(arg)
    val = arg;
  else
    val = eval(arg, 'arg');
  end


% Return non-zero if ARG looks like an option.
function p = getopt_option_p(arg)

  p = (ischar(arg) & numel(arg) > 1 & arg(1) == '-');


% Return non-zero if ARG is a valid color.
function [p, val] = getopt_color_p(arg)

  val = [];

  switch class(arg)
   case 'char'
    switch arg
     case {'y', 'yellow', ...
           'm', 'magenta', ...
           'c', 'cyan', ...
           'r', 'red', ...
           'g', 'green', ...
           'b', 'blue', ...
           'w', 'white', ...
           'k', 'black'};
      val = arg;
    end
   case 'double'
    switch numel(arg)
     case 1
      if arg >= 0 & arg <= 1
        val = repmat(arg, 1, 3);
      end
     case 3
      if all(arg >= 0 & arg <= 1)
        val = arg;
      end
    end
  end

  p = ~ isempty(val);


% Return non-zero if ARG is a valid cursor.
function [p, val] = getopt_cursor_p(arg)

  val = [];

  switch class(arg)
   case 'char'
    switch arg
     case {'crosshair', 'arrow', 'watch', 'topl', 'topr', 'botl', 'botr', ...
           'circle', 'cross', 'fleur', 'left', 'right', 'top', 'bottom', ...
           'fullcrosshair', 'ibeam'};
      val = arg;
    end
  end

  p = ~ isempty(val);


% Return non-zero if ARG is a valid line style.
function [p, val] = getopt_line_p(arg)

  val = [];

  if isempty(arg)
    arg = 'none';
  end

  switch class(arg)
   case 'char'
    switch arg
     case {'none', ...
           '-', '--', ':', '-.'}
      val = arg;
     case 'solid'
      val = '-';
     case 'dashed'
      val = '--';
     case 'dotted'
      val = ':';
     case 'dashdot'
      val = '-.';
    end
  end

  p = ~ isempty(val);


% Return non-zero if ARG is a valid marker.
function [p, val] = getopt_marker_p(arg)

  val = [];

  if isempty(arg)
    arg = 'none';
  end

  switch class(arg)
   case 'char'
    switch arg
     case {'none', ...
           '+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'}
      val = arg;
     case 'plus'
      val = '+';
     case 'circle'
      val = 'o';
     case {'asterisk', 'star'}
      val = '*';
     case 'point'
      val = '.';
     case {'cross', 'x-mark'}
      val = 'x';
     case 'square'
      val = 's';
     case 'diamond'
      val = 'd';
     case 'up'
      val = '^';
     case 'down'
      val = 'v';
     case 'left'
      val = '<';
     case 'right'
      val = '>';
     case 'pentagram'
      val = 'p';
     case 'hexagram'
      val = 'h';
    end
  end

  p = ~ isempty(val);


% Local test procedure.
function getopt_test

  for ind = 1:inf
    switch ind
     case 1
      opt(ind).name = '-flag';
     case 2
      opt(ind).name = '-mode';
      opt(ind).type = 'char';
      opt(ind).choices = {'box', 'line'};
      opt(ind).values = [4, 2];
      opt(ind).value = 4;
     case 3
      opt(ind).name = '-color';
      opt(ind).type = 'color';
     case 4
      opt(ind).name = '-cursor';
      opt(ind).type = 'cursor';
     case 5
      opt(ind).name = '-full';
      opt(ind).type = 'alias';
      opt(ind).index = 4;
      opt(ind).value = 'fullcrosshair';
     case 6
      opt(ind).name = '-fullcross';
      opt(ind).type = 'alias';
      opt(ind).index = 5;
     otherwise
      break;
    end
  end

  [opt, arg] = getopt(opt, ...
                      'foo', ...
                      '-flag', ...
                      '-mode=line', ...
                      magic(3), ...
                      '-color', [1 1 1], ...
                      '-fullcross', ...
                      'bar');

  opt(1)
  opt(2)
  opt(3)
  opt(4)
  arg


% More fun with aliases.
function getopt_test_cc

  % Initialize options.
  for ind = 1:inf
    switch ind
     case 1
      % Warning level.
      opt(ind).name = '-W';
     case 2
      % Turn on all warnings.
      opt(ind).name = '-Wall';
      opt(ind).type = 'alias';
      opt(ind).index = 1;
      opt(ind).value = inf;
     case 3
      % No warnings.
      opt(ind).name = '-w';
      opt(ind).type = 'alias';
      opt(ind).index = 1;
      opt(ind).value = 0;
     case 4
      % Enable 'asm' keyword.
      opt(ind).name = '-fasm';
     case 5
      % Disable 'asm' keyword.
      opt(ind).name = '-fno-asm';
      opt(ind).type = 'alias';
      opt(ind).index = 4;
      opt(ind).value = 0;
     otherwise
      break;
    end
  end


% local variables:
% time-stamp-line-limit: 200
% time-stamp-format: "%:y-%02m-%02d %02H:%02M:%02S %Z"
% end:

%% getopt.m ends here
