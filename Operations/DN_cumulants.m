function out = DN_cumulants(y,thecum)
% Very simple function that uses the skewness and kurtosis functions in Matlab's Statistics Toolbox to calculate these higher order moments of input time series, y
% The appropriate higher order moment is specified with the second input, thecum
% Ben Fulcher, 2008

if nargin < 2 || isempty(thecum)
    error('You must specify the cumulant type as the second input!')
end

switch thecum
	case 'skew1' % skewness
		out = skewness(y);
    case 'skew2' % corrects for bias
        out = skewness(y,0);
	case 'kurt1' % kurtosis
		out = kurtosis(y);
    case 'kurt2' % corrects for bias
        out = kurtosis(y,0);
    otherwise
        error(sprintf('Unknown cumulant specifier: %s',thecum))
end

end