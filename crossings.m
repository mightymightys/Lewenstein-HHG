function [ci,ti] = crossings(x,t,val,opt)
%
% CROSSINGS returns the crossings of a given value in a vector. Note that
% it does not return "occurrences" of that value, unless the value is
% actually crossed.
%
% SYNTAX
%   [I,Ti] = CROSSINGS(X,T,VALUE,'option').
%   Input arguments T, VALUE and 'option' can be defaulted by assigning 
%   them the [] value. For instance, 
%   >> [I,Ti] = CROSSINGS(X,[],[],'option'), 
%   will assign the default values to T and VALUE
%
% INPUTS 
%   X        vector whose crossings are sought
%   T        time coordinates of X. Obviously length(T) must equal 
%            length(X). Default is T = 1:length(X).
%   VALUE:   arbitrary number whose crossings are sought. Default is
%            VALUE = 0 (zero-crossings). 
%   'option' when there is an actual VALUE-crossing (and not just a VALUE
%            value), 'option' determines to which element (the one before
%            or the one after the crossing) of X the crossing is assigned
%            to. Possible values are:  
%               'dis': the VALUE-crossing is assigned to the element which
%                      is closer to VALUE ('d'istance criterion). If
%                      abs(X(n)-VALUE)==abs(X(n+1)-VALUE), then the
%                      VALUE-crossing is assigned to n+1. 
%               'pre': the VALUE-crossing is assigned to the element
%                      before the actual VALUE-crossing (the 'p'receding sample)
%               'fol': the VALUE-crossing is assigned to the element
%                      after the actual VALUE-crossing (the 'f'ollowing sample)
%               'int': the VALUE-crossing is assigned to a number after
%                      interpolating the element number before and after 
%                      the VALUE-crossing. 
%                      NOTE: when using the 'int' option you must ask for
%                      both output arguments (see description of I in the
%                      'OUTPUT' section) 
%               EXAMPLE: assume that VALUE=0, X(4)=-2 and X(5)=3. Then
%                        'dis' will assign the 0-crossing to i=4
%                        'pre' will assign the 0-crossing to i=4
%                        'fol' will assign the 0-crossing to i=5
%                        'int' will assign the 0-crossing to t=4.4
%             Default is 'dis'.  
%
% OUTPUTS
%   I   indices where the crossings occur (only integer values). If the
%       'int' option is used, then I will contain the indices as if the
%       'dis' option were used and the fractional numbers will be contained
%       in Ti.
%   Ti  the respective time instants of the crossings.
%
% EXAMPLES:
% Assume that the VALUE = 0.
%   X = [1 -1 0 -1]: CROSSING will return the 0-crossing from 1 to -1. 
%                    The sequence -1, 0, -1 is not considered a crossing (0
%                    is only reached but not crossed). 
%   X = [-1 0 0 -1]: CROSSING will not return any 0-crossings (0 is only
%                    reached but not crossed).  
%   X = [1 0 0 -1] : CROSSING will return only one 0-crossing at element 3
%                    (from 0 to -1)
% 
% Occurrences of VALUE at the beginning or end of X are ignored (since no
% actual crossings occur).
%
% Christos Saragiotis, 2007-08-29
% Copyright (c) Christos Saragiotis, 2007
% christos dot saragiotis at gmail dot com
%
% Modifications
%  2007/10/05: Corrected a concatenation error when input data were column
%              vectors 
%  



% check the number of input arguments
narginchk(1,4);

if nargin == 4
    if strcmp(opt,'dis') && strcmp(opt,'pre') && strcmp(opt,'fol') && strcmp(opt,'int') 
        opt = 'dis';
    end
else if nargin < 4
        opt = 'dis';
        if nargin < 3
            val = 0;
            if nargin < 2
                t = 1:length(x);
            end
        end
    end
end
        

if isempty(opt), opt = 'd'; end
if isempty(val), val = 0; end
if isempty(t), t = 1:length(x); end

if length(t) ~= length(x)
    error('X and T must have the same length!');
end

iscolumn = ( size(x,1)>1 );
x = x(:)'; t = t(:)';

x = x - val;

i_0 = find(x == 0);                     % indices of 0's in x

% Ignore consecutive 0's if they occur at the beginning of x
N = 1;
while ~isempty(i_0) && (x(N) == 0) 
    i_0(1) = []; 
    N = N + 1;
end 

% Ignore consecutive 0's if they occur at the end of s
N = length(x);
while ~isempty(i_0) && (x(N) == 0)
    i_0(end) = [];
    N = N-1;
end


ci_0 = [];
if ~isempty(i_0)
    prec_0_sgn = sign(x(i_0-1));        % sign of elements preceding a 0
    prec_0s_i = find(prec_0_sgn == 0);  % indices of elements that precede 
                                        % a 0 and are 0 also 

    i_0( prec_0s_i-1 ) = [];            % elimination of 0's preceding also 
    prec_0_sgn(prec_0s_i) = [];         % 0's

    foll_0_sgn = sign(x(i_0+1));        % There are no zero values in 
                                        % foll_sgn due to previous 4 lines

    temp = prec_0_sgn.*foll_0_sgn;     
    ci_0 = i_0( find( temp<0 ) );       % crossing indices;
end


ci_1 = find( (x(1:end-1).* x(2:end)) < 0 );

switch opt
    case 'dis'
        absDiffX = abs(x(ci_1+1))-abs(x(ci_1));
        dist = find(absDiffX <= 0);
        ci_1(dist) = ci_1(dist) + 1;
        t_c = t(ci_1);
    case 'int'
        dt = t(ci_1+1) - t(ci_1);
        dx = x(ci_1+1) - x(ci_1);
        t_c = t(ci_1) - x(ci_1).*dt./dx;

        absDiffX = abs(x(ci_1+1))-abs(x(ci_1));
        dist = find(absDiffX <= 0);
        ci_1(dist) = ci_1(dist) + 1;
    case 'fol'
        ci_1 = ci_1 + 1;
        t_c = t(ci_1);
    case 'pre'
        t_c = t(ci_1);
end

ci = sort([ci_0 ci_1]);
ti = sort([t(ci_0) t_c]);

if iscolumn
    ci = ci'; ti = ti';
end

