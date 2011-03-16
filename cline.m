%
% Draw a color-coded line by using the edge of a patch with no facecolor
%
% SYNTAX
% ======
% h = cline(x, y [, z, cdata])
%
% INPUT
% =====
% x                     vector with x-values
% y                     vector with y-values
% z (opt.)              vector with z-values
% cdata (opt.)          vector with color-data
%
% 2 input arguments =>  cdata = y; z=0      % s. Example 1
% 3 input arguments =>  cdata = z           % s. Example 2
% 4 i.a. & z = []   =>  cdata = y; z=0      % s. Example 4
%
% OUPUT
% =====
% h                 Handle to line (i.e. patch-object !!!)
%
% Examples
% ========
% t = 2*pi:.1:8*pi;
%
% cline(sqrt(t).*sin(t), sqrt(t).*cos(t)); view(3)                       % Example 1
% cline(sqrt(t).*sin(t), sqrt(t).*cos(t), t); view(3)                    % Example 2
% cline(sqrt(t).*sin(t), sqrt(t).*cos(t), t, rand(size(t))); view(3)     % Example 3
% cline(sqrt(t).*sin(t), sqrt(t).*cos(t), [], rand(size(t))); view(3)	 % Example 4
% 
%
% Author & Version
% ================
% S. Hölz, TU-Berlin, seppel_mit_ppATweb.de
% V 1.0, 16.4.2007
% Created using Matlab 7.0.4 (SP2)
%

% Info
% ====
% This function uses the edges of a patch to represent the colored 2D/3D-line. The marker-related
% properties (i.e. 'maker','markersize','markeredgecolor','markerfacecolor') can be used as with a
% regular line. 
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% The line-related properties (i.e. 'linestyle','linewidth') WILL HAVE NO EFFECT 
% while displaying the line on screen, but will change the output when printing to file !!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!
%
% This is not a flaw of this function but rather the way Matlab interprets the interpolated edges of
% a patch on screen and when printing. I do not know, if this behavior is consistent in other
% versions of Matlab.
%

% License 
% =======
% Copyright (c) 2009, Sebastian Hölz
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

function h = cline(x, y, z, cdata)

    % Check input arguments
    error(nargchk(2, 4, nargin))
    if ~isnumeric(x) || ~isnumeric(y) || ~isvector(x) || ~isvector(y) || length(x)~=length(y); 
        error('x and y must be numeric and conforming vectors'); 
    end
    if (nargin == 3 && (~isnumeric(z) || ~isvector(z) || length(x)~=length(z))) || ...
       (nargin == 4 && ~isempty(z) && (~isnumeric(z) || ~isvector(z) || length(x)~=length(z))) ...
       (nargin == 4 && (~isnumeric(cdata) || ~isvector(cdata) || length(x)~=length(cdata)))
        error('z (and cdata) must be a numeric vector and conforming to x and y'); 
    end
    
    % Draw line as patch
    if nargin == 2
        p = patch([x(:)' nan], [y(:)' nan], 0);
        cdata = [y(:)' nan];
    elseif nargin == 3
        p = patch([x(:)' nan], [y(:)' nan], [z(:)' nan], 0);
        cdata = [z(:)' nan];
    elseif nargin == 4
        if isempty(z); z = zeros(size(x)); end
        p = patch([x(:)' nan], [y(:)' nan], [z(:)' nan], 0);
        cdata = [cdata(:)' nan];
    end
    
    set(p,'cdata', cdata, 'edgecolor','interp','facecolor','none')
    
    
    % Create output
    if nargout == 1; h = p; end