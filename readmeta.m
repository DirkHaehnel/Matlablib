function output = readmeta(filename, cmdno, argno)
%READMETA	Read old matlab meta file.	
%	READMETA can produce a listing of the commands and arguments for the
%	named meta file or return a specified argument to one of the plot
% 	commands found in the meta file.  If readmeta has only one argument
%	(the meta file name) then an annotated listing of the commands and
%	their data is produced.  The command number found on the annotated
%	list can be used with an argument number to return as an output the
%	desired data.
%	readmeta('matlab.met')              - produces listing of matlab.met
%	out = readmeta('matlab.met',i,j)    - returns j'th argument to i'th
%	                                      command found in matlab.met

if (nargout > 1) | ((nargout == 1) & (nargin == 1))
	error('Too many outputs for READMETA'), 
end
if (nargin ~= 1) & (nargin ~= 3), 
	error('Incorrect number of inputs for READMETA'), 
end

if ~isstr(filename),
	error('First argument to READMETA must be the meta file name'),
end
if nargin > 1
	if any(size(cmdno) ~= 1)
		error('Second argument to READMETA must be a scalar'),
	end
	if any(size(argno) ~= 1)
		error('Third argument to READMETA must be a scalar'),
	end
end

% old matlab function table
fnames = [
'plot    ';'loglog  ';'semilogx';'semilogy';'polar   ';'mesh    ';
'title   ';'xlabel  ';'ylabel  ';'text    ';'axis    ';'grid    ';
'prtscr  ';'hold    ';'subplot ';'shg     ';'cla     ';'clg     ';
'hog     ';'sha     ';'aspect  ';'contour ';'polyline';'ginput  ';
'polymark'
];

fnum = [600:604 606:615 621:630];

lastfun = 1000;

% open meta file
fp = fopen(filename,'r');
if fp == -1, error('Can''t open desired meta file'), end

% loop over meta commands
cmdcount = 0;
while 1
	[ia,count] = fread(fp, 5, 'int');
	if count ~= 5, break, end
	if ia(1) < 600  |  ia(1) >= 1000
		% read in axis limits
		[lim,count] = fread(fp, 4, 'double');
		if count ~= 4, break, end

	else
		% display command name and data
		cmdcount = cmdcount + 1;
		if nargin == 1
			fprintf(1,'%d - %s\n',cmdcount,fnames(fnum==ia(1),:));
		end

		% loop on command arguments
		for i=1:ia(2)
			[mathead,count] = fread(fp, 4, 'long');
			if count ~= 4, fclose(fp), return, end
			m = mathead(1);
			n = mathead(2);
			str_flag = mathead(3);
			cmplx_flag = mathead(4);
			[rdata,count] = fread(fp, [m n], 'double');
			if count ~= m*n, fclose(fp), return, end
			if cmplx_flag
				[idata,count] = fread(fp, [m n], 'double');
				if count ~= m*n, fclose(fp), return, end
			end
			if nargin == 1
				fprintf(1,'   arg - %d - ',i);
				if str_flag
					fprintf(1,'string - ''%s''\n',setstr(rdata'));
				else
					fprintf('numeric - %d x %d, ',m,n);
					if cmplx_flag
						fprintf(1,'complex\n');
					else
						fprintf(1,'\n');
					end
				end
			elseif (cmdcount == cmdno)  &  (argno == i)
				output = rdata;
				if cmplx_flag, output = output + idata*sqrt(-1); end
				if str_flag, output = setstr(output); end
				fclose(fp);
				return;
			end
		end
		if nargin == 1, fprintf(1,'\n'); end
	end
end

fclose(fp);

if nargin > 1,
	error('READMETA could not find the desired command/argument'),
end


------ Message Header Follows ------
Received: from windsurf.UUCP by MathWorks.COM with UUCP id AA23683
  (5.65c/IDA-1.5-dp for loren); Tue, 30 Nov 1993 18:46:32 -0500
Received: by windsurf.UUCP (4.1/SMI-4.1)
	id AA21132; Tue, 30 Nov 93 15:33:48 EST
Date: Tue, 30 Nov 93 15:33:48 EST
From: mathw@windsurf.mathworks.com (Steve Bangert)
Message-Id: <9311302033.AA21132@windsurf.UUCP>
To: loren@MathWorks.COM
Subject: meta files





