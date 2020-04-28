function [ nt ] = bateman( n0, rate, time )
% solve one path, first order decay reactions with Bateman's equations

% n0 - initial number of first isotope in decay chain
% k - 1d vector of decay chain rate constants
% t - 1d vector of time values

%check for column vectors   
if ~iscolumn(rate)
    rate = rate';
end
    
if ~iscolumn(time)
    time = time';
end

d = length(rate);
t = length(time);

% rate constant-weighted, time-dependent decay terms
rate_time = repmat(rate', t, 1);
decay = n0 * rate_time .* exp(-rate_time .* repmat(time, 1, d));

% 2d matrix of rate constant differences
rate_diff = repmat(rate, 1, d) - repmat(rate', d, 1);

% 2d matrix of scaled rate constants
% eps needed to prevent division by zero
rate_scale = repmat(rate, 1, d) ./ (rate_diff + diag(eps));

% now remove large diagonal values due to dviding by eps
rate_scale = rate_scale - diag(diag(rate_scale));

% 2d matrix of scaled rate cumulative products
rate_prod = cumprod((rate_scale + eye(d)), 1);

% rate factors to multiply by decay
rate_facs = (rate_prod - triu(rate_prod, 1)) ./ repmat(rate, 1, d);

% 2d matrix of time-dependent number of isotopes
nt = decay * rate_facs';
end