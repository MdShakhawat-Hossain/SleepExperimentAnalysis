function out = OXY_HRF(s,x,Fr)
    % OXY_HRF -- Use Analytic Methods to compute HRF functions of the observed
    %       changes in oxygen signals based on the assumption that the
    %       oxygen response to locomotion is a linear, time-invariant system.
    %
    % USAGE:
    %       h = OXY_HRF(s,x,Fs);
    %
    % INPUTS:
    %       s: binarized locomotion signal, sample by 1
    %       x: oxygen time series, sample by 1
    %       Fr: sample rate, Hz, currently we only analyze 3 Hz data
    %
    % OUTPUTS:
    %       HRF: impulse response
    %       timeShift: shift input data by timeShift seconds to prevent
    %                  boundary artifacts
    %       IRLength: the length of impulse response, in samples
    %
    % Last Modified: 05-18-2017
    % Modified By: Qingguang Zhang (qingguang.zhang@gmail.com)
    %
    % See also:
    %       ImpResp
    %       HRFfeatures

    % Reference:
    % Bingxing Huo et al., Quantitative spearation of arterial and venous
    % cerebral blood volume increases during voluntary locomotion. NeuroImage
    % 105 (2015): 369-379.

    % remove trend
    x = detrend(x);
    
    % make sure x is a column vector
    if size(x,1) < size(x,2), x = x(:); end
    
    % Make sure there is no NaN data in binary speed data
    idx = find(isnan(s));     s(idx) = 0;     clear idx;
    % make sure s is a column vector
    if size(s,1) < size(s,2), s = s(:);  end
    s = detrend(s); % won't affect the shape of HRF
    
    timeShift = 10; % shift stimulus backward by 10 seconds to avoid boundary effects
    IRLength = 70*Fr; % length of impulse response fuction, 70 seconds

    %% Compute HRF directly
    HRF = ImpResp(s,x,Fr,timeShift,IRLength); % HRF contains timeshift data
    HRF = HRF(1:IRLength);
    HRF_time = -timeShift:1/Fr:(IRLength-timeShift*Fr-1)/Fr;
    HRF_time = HRF_time(:);


%     %% Get HRF features: peak, AUC and onsettime
%     out1 = HRFfeatures(HRF_time, HRF, Fr, timeShift, IRLength);
%     
%     out = out1; 
    out.timeShift = timeShift;
    out.IRLength = IRLength;
    out.HRF = HRF;
    out.HRF_time = HRF_time;
end


function H = ImpResp(datain,dataout,Fr, timeShift, IRLength)
% ImpResp: compute impulse response using numerically analysis, 
%          this script works for single-input-single-output (SISO) system
%
% INPUTS:
%       datain: stimulus data (input) for the model
%       dataout: response data (output) for the model, sample by 1
%       Fr: sample rate for the time series
%       timeShift: time shift for the computing model, in seconds
%       IRLength: the length of the impulse response, in samples
%
% OUTPUTS:
%       H: impulse response, samples by 1
%
% USUAGE:
%       H = ImpResp(dataStimulus, dataResponse, timeShift, IRLength);
%
% Last Modified: 05-18-2017 
% Modified By: Qingguang Zhang (qingguang.zhang@gmail.com)

if nargin < 2
    error('HRFNumerical -> Not enough inputs');
elseif size(datain,1)~=size(dataout,1)
    error('Input and output should be of same lengths.')
end

if nargin < 3
    Fr = 3; % default frame rate, 3 Hz
    timeShift = 10 * Fr; % default time shift, 10 seconds, convert to samples
    IRLength = round((length(datain) - timeShift)/5); % Length of the impulse response
elseif nargin < 4
    timeShift = 10 * Fr; % shift 10 seconds data, convert to samples
    IRLength = round((length(datain) - timeShift)/5); % Length of the impulse response
else
    timeShift = timeShift * Fr;
end

datain = datain(timeShift+1:end); % shift input time backward, assuming non-causal system
dataout=dataout(1:end-timeShift,:); % truncate output to keep same lengths
% Note: defination of non-causal system: the current value does not
% depend on previous values

N = length(datain); % length of processed input data

% Generate a numerical matrix for the input data
c = zeros(2*N-1,1);
c(1:N) = datain; % this is the first column of the Toeplitz Matrix

r = zeros(IRLength,1);
r(1) = datain(1); % this is the first column of the Toeplitz Matrix

TM = toeplitz(c,r); % generate Toeplitz Matrix (TM)

S = [ones(2*N-1,1), TM]; % the numerical matrix for input data, s

% Generate a numerical matrix for the output data
X = zeros(2*N-1, 1);
X(1:N) = dataout(1:N);

% % Method 1:
% H = inv(S'*S)*S'*X; % To avoid rank deficiency, try to use this part of
%                     % the code; otherwise, use Method 2, which is
%                     % computational efficient

% Method 2:
H = S\X; % this is more acurate and faster than H = inv(S'*S)*S'*X;
end

function out1 = HRFfeatures(HRF_time, HRF, Fr, timeShift, IRLength)
% HRFfeatures: get HRF features, e.g., AUC, peak, onsettime
%
% INPUTS:
%       HRF_time: time stamp for HRF
%       HRF: oxygen HRF
%       Fr: sample rate for the time series
%       timeShift: time shift for the computing model, in seconds
%       IRLength: the length of the impulse response, in samples
%
% OUTPUTS:
%       out1: structure contains all the HRF features
%
% USUAGE:
%       out1 = HRFfeatures(HRF_time, HRF, Fr, timeShift, IRLength)
%
% Last Modified: 05-18-2017 
% Modified By: Qingguang Zhang (qingguang.zhang@gmail.com)


% Interpolate HRF to get a resolution of 1 ms
step = 1/1000;

HRF_time_interp = -timeShift:step:(IRLength-timeShift*Fr-1)/Fr;
HRF_time_interp = HRF_time_interp(:);

% turn off warning to prevent MATLAB throwing messages
warning('off','MATLAB:interp1:EmptyMethod');
HRF_interp = interp1(HRF_time, HRF, HRF_time_interp,'spline');

% output
out1.HRF_time_interp = HRF_time_interp;
out1.HRF_interp = HRF_interp;

% 1. find peak value of HRF
pos_peak_st = 0; % start time for detecting positive peak, in seconds
pos_peak_fin = 10; % end time for detecting positive peak, in seconds
neg_peak_st = 0; % start time for detecting negative peak, in seconds
neg_peak_fin = 10; % end time for detecting negative peak, in seconds

% 1.1 find positive peak
st = find(HRF_time_interp >= pos_peak_st, 1, 'first');
fin = find(HRF_time_interp <= pos_peak_fin, 1, 'last');
[~,HRFpeakIdxtmp] = max(abs(HRF_interp(st:fin))); % peak value
HRFpeakIdx.pos = HRFpeakIdxtmp + st-1;
if HRF_interp(HRFpeakIdx.pos) > 0 % make sure it is positive
    HRFpeaktime.pos = HRF_time_interp(HRFpeakIdx.pos); % peak time
    HRFpeak.pos = HRF_interp(HRFpeakIdx.pos); % peak amplitude
else % not positive
    HRFpeakIdx.pos = NaN;
    HRFpeaktime.pos = NaN;
    HRFpeak.pos = NaN;
end

% 1.2 find negative peak
neg_peak_fin = min([neg_peak_fin, HRFpeaktime.pos]); % pick the smaller one
st = find(HRF_time_interp >= neg_peak_st, 1, 'first');
fin = find(HRF_time_interp <= neg_peak_fin, 1, 'last');
if min(HRF_interp(st:fin)) >= 0 % there is no negative peak
    HRFpeakIdx.neg = NaN;
    HRFpeaktime.neg = NaN; % peak time
    HRFpeak.neg = NaN; % peak amplitude
else
    [~,HRFpeakIdxtmp] = min(HRF_interp(st:fin)); % peak value
    HRFpeakIdx.neg = HRFpeakIdxtmp + st-1;
    HRFpeaktime.neg = HRF_time_interp(HRFpeakIdx.neg); % peak time
    HRFpeak.neg = HRF_interp(HRFpeakIdx.neg); % peak amplitude
end

% 1.3 find zero-crossing point
% zero-crossing between positive and negative peaks
if isnan(HRFpeakIdx.neg) || isnan(HRFpeakIdx.pos) % no negative or positive peaks, no crossing
    HRFcxIdx.btwn = NaN;
    HRFcxtime.btwn = NaN;    
else
    HRF_seg = HRF_interp(HRFpeakIdx.neg:HRFpeakIdx.pos);
    HRF_seg_time = HRF_time_interp(HRFpeakIdx.neg:HRFpeakIdx.pos);
    tmp = find(HRF_seg >= 0, 1, 'first');
    HRFcxIdx.btwn = tmp+HRFpeakIdx.neg-1; % zero crossing index between peaks
    HRFcxtime.btwn =  HRF_seg_time(tmp); % zero crossing time between peaks
end

% zero-crossing after positive peaks
if ~isnan(HRFpeakIdx.pos) % positive peak exist
    HRF_seg = HRF_interp(HRFpeakIdx.pos:end);
    HRF_seg_time = HRF_time_interp(HRFpeakIdx.pos:end);
    tmp = find(HRF_seg <= 0, 1, 'first');
    HRFcxIdx.post = tmp+HRFpeakIdx.pos-1; % zero crossing index after positive peak
    HRFcxtime.post =  HRF_seg_time(tmp); % zero crossing time after positive peak
elseif ~isnan(HRFpeakIdx.neg) % no positive peak, negative peak exist
    HRFcxIdx.post
end

% real time 0
HRFcxIdx.first = find(HRF_time_interp == 0);

% output
out1.HRFpeakIdx.neg = HRFpeakIdx.neg;
out1.HRFpeaktime.neg = HRFpeaktime.neg;
out1.HRFpeak.neg = HRFpeak.neg;

out1.HRFpeakIdx.pos = HRFpeakIdx.pos;
out1.HRFpeaktime.pos = HRFpeaktime.pos;
out1.HRFpeak.pos = HRFpeak.pos;

out1.HRFcxIdx.first = HRFcxIdx.first;
out1.HRFcxtime.first = 0;

out1.HRFcxIdx.btwn = HRFcxIdx.btwn;
out1.HRFcxtime.btwn = HRFcxtime.btwn;

out1.HRFcxIdx.post = HRFcxIdx.post;
out1.HRFcxtime.post = HRFcxtime.post;


% 2. find area under the curve (AUC) for the negative and positive peak
if HRFcxIdx.first == HRFcxIdx.btwn || isnan(HRFcxIdx.btwn)
    HRF_AUC.neg = 0;
else
    HRF_AUC.neg = trapz(HRF_time_interp(HRFcxIdx.first : HRFcxIdx.btwn),... 
                        HRF_interp(HRFcxIdx.first : HRFcxIdx.btwn));
end

if HRFcxIdx.btwn == HRFcxIdx.post || isnan(HRFcxIdx.btwn) || isnan(HRFcxIdx.post)
    HRF_AUC.pos = 0;
else
    HRF_AUC.pos = trapz(HRF_time_interp(HRFcxIdx.btwn : HRFcxIdx.post),...
                        HRF_interp(HRFcxIdx.btwn : HRFcxIdx.post));
end

% output
out1.HRF_AUC.neg = HRF_AUC.neg;
out1.HRF_AUC.pos = HRF_AUC.pos;

% 3. find onset time for the negative and positive peak
[onsettime.pos, slope.pos, intercept.pos] = ...
    Oxy_findOnsetTime(HRF_time_interp(HRFcxIdx.btwn : HRFpeakIdx.pos), HRF_interp(HRFcxIdx.btwn : HRFpeakIdx.pos));

if ~isnan(HRFpeakIdx.neg)
    [onsettime.neg, slope.neg, intercept.neg] = ...
        Oxy_findOnsetTime(HRF_time_interp(HRFcxIdx.first : HRFpeakIdx.neg), HRF_interp(HRFcxIdx.first : HRFpeakIdx.neg));
else
    onsettime.neg = NaN;
    slope.neg = NaN;
    intercept.neg = NaN;
end

% output
out1.onsettime.pos = onsettime.pos;
out1.slope.pos = slope.pos;
out1.intercept.pos = intercept.pos;

out1.onsettime.neg = onsettime.neg;
out1.slope.neg = slope.neg;
out1.intercept.neg = intercept.neg;

end

function [onsettime, slope, intercept] = Oxy_findOnsetTime(time, data)
% time: time points for rising segment
% data: data points for rising segment
if ~isempty(data)
    peakamp = data(end);
    peak80 = peakamp*0.8; % 80% peak value
    peak20 = peakamp*0.2; % 20% peak value

    if peakamp <= 0
        idx20 = find(data <= peak20, 1, 'first');
        idx80 = find(data <= peak80, 1, 'first');
    else
        idx20 = find(data >= peak20, 1, 'first');
        idx80 = find(data >= peak80, 1, 'first');
    end

    % compute onset time
    X = [ones(length(time(idx20:idx80)),1) time(idx20:idx80)];
    Y = data(idx20:idx80);
    b = X\Y; % b has two elements, first one is intercept, second one is slope

    onsettime = -b(1)/b(2);
    slope = b(2);
    intercept = b(1);
else
    onsettime = NaN;
    slope = NaN;
    intercept = NaN;
end
end