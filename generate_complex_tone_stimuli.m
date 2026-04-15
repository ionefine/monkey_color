function [stimuli, t] = generate_complex_tone_stimuli(config)
%GENERATE_COMPLEX_TONE_STIMULI Create one or two spatialized complex-tone stimuli.
%
% [stimuli, t] = generate_complex_tone_stimuli(config)
%
% This function generates one or two 200 ms complex-tone stimuli with:
%   - Fundamental frequency selectable (e.g., 600 Hz or 5000 Hz)
%   - Pitch trajectory: 'constant', 'rising', or 'descending'
%   - Timbre: 'flute' or 'oboe'
%   - Spatial position at x = -135, 0, or +135 deg (y = 0, z = 1 m) via ITD + ILD
%   - Raised-cosine onset/offset ramps of 25 ms
%
% INPUT
%   config: struct with fields:
%       .fs                Sampling rate in Hz (default 48000)
%       .duration          Duration in seconds (default 0.200)
%       .rampDuration      Ramp duration in seconds (default 0.025)
%       .nStimuli          Number of stimuli, 1 or 2 (default 2)
%       .stim              1xN struct array (N = nStimuli) with fields:
%       .snd               Backward-compatible alias for .stim
%           .timbre        'flute' or 'oboe'
%           .f0            Fundamental frequency in Hz
%           .trajectory    'constant' | 'rising' | 'descending'
%           .glideRatio    End/start F0 ratio for glides (default 1.2)
%           .loc           Source azimuth: +135 (left), -135 (right), 0 (center)
%           .azimuthDeg    Backward-compatible alias for .loc
%           .distanceM     Source distance in meters (default 1)
%           .level         Linear gain scalar (default 0.2)
%           .phaseMode     'random' or 'sine' (default 'random')
%
% OUTPUT
%   stimuli: 1xN struct array with fields:
%       .mono              Nx1 mono waveform
%       .stereo            Nx2 binaural waveform [L R]
%       .meta              Copy of parameters used for this stimulus
%   t: Nx1 time vector in seconds
%
% EXAMPLE
%   cfg.nStimuli = 2;
%   cfg.stim(1) = struct('timbre','flute','f0',600,'trajectory','constant', ...
%                        'loc',135,'distanceM',1,'level',0.25);
%   cfg.stim(2) = struct('timbre','oboe','f0',5000,'trajectory','descending', ...
%                        'glideRatio',1.3,'loc',0,'distanceM',1,'level',0.2);
%   [stimuli, t] = generate_complex_tone_stimuli(cfg);
%   sound(stimuli(1).stereo, 48000);
%
% NOTE
%   The spatialization uses a lightweight naturalistic approximation:
%   - ITD from a spherical-head model
%   - ILD via frequency-independent head-shadow approximation by azimuth

if nargin < 1
    config = struct();
end

cfg = applyDefaults(config);

if ~(cfg.nStimuli == 1 || cfg.nStimuli == 2)
    error('config.nStimuli must be 1 or 2.');
end

nSamples = round(cfg.duration * cfg.fs);
t = (0:nSamples-1)' / cfg.fs;
env = raisedCosineEnvelope(nSamples, round(cfg.rampDuration * cfg.fs));

stimuli = repmat(struct('mono', [], 'stereo', [], 'meta', []), 1, cfg.nStimuli);

for k = 1:cfg.nStimuli
    p = cfg.stim(k);

    validateStim(p, k);

    f0_t = instantaneousF0(p.f0, p.trajectory, p.glideRatio, t);

    x = synthComplexTone(f0_t, cfg.fs, p.timbre, p.phaseMode);
    x = x .* env;
    x = p.level * x;

    azimuthDeg = stimulusAzimuth(p);
    [left, right] = binauralSpatialize(x, cfg.fs, azimuthDeg, p.distanceM);

    % Hard limit to avoid accidental clipping if multiple tones are mixed later.
    peak = max(abs([left; right]));
    if peak > 0.999
        left = left / peak * 0.999;
        right = right / peak * 0.999;
    end

    stimuli(k).mono = x;
    stimuli(k).stereo = [left, right];
    stimuli(k).meta = p;
end

end

% ------------------------------ Helpers ---------------------------------

function cfg = applyDefaults(cfg)
if ~isfield(cfg, 'fs'),           cfg.fs = 48000;      end
if ~isfield(cfg, 'duration'),     cfg.duration = 0.200; end
if ~isfield(cfg, 'rampDuration'), cfg.rampDuration = 0.025; end
if ~isfield(cfg, 'nStimuli'),     cfg.nStimuli = 2;     end
if ~isfield(cfg, 'stim') && isfield(cfg, 'snd')
    cfg.stim = cfg.snd; % Backward-compatible alias.
end

if ~isfield(cfg, 'stim') || isempty(cfg.stim)
    cfg.stim = repmat(defaultStim(), 1, cfg.nStimuli);
else
    for i = 1:cfg.nStimuli
        d = defaultStim();
        if i <= numel(cfg.stim)
            u = cfg.stim(i);
            fn = fieldnames(d);
            for j = 1:numel(fn)
                if ~isfield(u, fn{j}) || isempty(u.(fn{j}))
                    u.(fn{j}) = d.(fn{j});
                end
            end
            cfg.stim(i) = u;
        else
            cfg.stim(i) = d;
        end
    end
end
end

function s = defaultStim()
s = struct( ...
    'timbre', 'flute', ...
    'f0', 600, ...
    'trajectory', 'constant', ...
    'glideRatio', 1.2, ...
    'loc', 135, ...
    'azimuthDeg', 135, ...
    'distanceM', 1, ...
    'level', 0.2, ...
    'phaseMode', 'random');
end

function validateStim(p, idx)
if ~ismember(lower(p.timbre), {'flute','oboe'})
    error('stim(%d).timbre must be ''flute'' or ''oboe''.', idx);
end
if ~ismember(lower(p.trajectory), {'constant','rising','descending'})
    error('stim(%d).trajectory must be ''constant'', ''rising'', or ''descending''.', idx);
end
if ~(isscalar(p.f0) && p.f0 > 0)
    error('stim(%d).f0 must be a positive scalar.', idx);
end
azimuthDeg = stimulusAzimuth(p);
if ~(isscalar(azimuthDeg) && any(abs(azimuthDeg - [-135, 0, 135]) < 1e-6))
    error('stim(%d).loc must be one of -135, 0, or +135 degrees.', idx);
end
if ~(isscalar(p.distanceM) && p.distanceM > 0)
    error('stim(%d).distanceM must be a positive scalar.', idx);
end
if ~(isscalar(p.level) && p.level > 0)
    error('stim(%d).level must be a positive scalar.', idx);
end
if ~ismember(lower(p.phaseMode), {'random','sine'})
    error('stim(%d).phaseMode must be ''random'' or ''sine''.', idx);
end
end

function f0_t = instantaneousF0(f0, trajectory, glideRatio, t)
switch lower(trajectory)
    case 'constant'
        f0_t = f0 * ones(size(t));
    case 'rising'
        f0_t = f0 * (glideRatio .^ (t / max(t(end), eps)));
    case 'descending'
        f0_t = f0 * ((1/glideRatio) .^ (t / max(t(end), eps)));
    otherwise
        error('Unknown trajectory: %s', trajectory);
end
end

function azimuthDeg = stimulusAzimuth(p)
% Prefer .loc if present, otherwise use legacy .azimuthDeg.
hasLoc = isfield(p, 'loc') && ~isempty(p.loc);
hasAz = isfield(p, 'azimuthDeg') && ~isempty(p.azimuthDeg);

if hasLoc
    azimuthDeg = p.loc;
elseif hasAz
    azimuthDeg = p.azimuthDeg;
else
    azimuthDeg = 135;
end
end

function x = synthComplexTone(f0_t, fs, timbre, phaseMode)
n = numel(f0_t);
t = (0:n-1)' / fs;

fMax = 0.48 * fs;
maxHarm = max(1, floor(fMax / max(f0_t)));

weights = harmonicWeights(maxHarm, timbre);
weights = weights(:) / (sum(weights) + eps);

if strcmpi(phaseMode, 'random')
    phi = 2*pi*rand(maxHarm,1);
else
    phi = zeros(maxHarm,1);
end

% Integrate instantaneous angular frequency to phase per harmonic.
basePhase = 2*pi*cumsum(f0_t)/fs;

x = zeros(n,1);
for h = 1:maxHarm
    x = x + weights(h) * sin(h * basePhase + phi(h));
end

% Remove DC and normalize.
x = x - mean(x);
x = x / (max(abs(x)) + eps);

% Timbre post-shaping to increase perceptual contrast.
if strcmpi(timbre, 'flute')
    % Smooth and suppress high-frequency buzz (flute-like purity).
    x = filter([1 1 1 1]/4, 1, x);
    x = filter([1 1 1 1]/4, 1, x);
    x = x + 0.01*randn(size(x)); % faint airy component
else
    % Brighter and reedy: boost high-frequency detail + soft saturation.
    x = x + 0.16 * [0; diff(x)];
    x = tanh(1.6 * x);
end

x = x / (max(abs(x)) + eps);
end

function w = harmonicWeights(maxHarm, timbre)
h = (1:maxHarm)';

switch lower(timbre)
    case 'flute'
        % Dominant low harmonics, very fast rolloff for a pure tone color.
        w = exp(-1.5*(h-1));
        if maxHarm >= 2, w(2) = w(2) * 0.65; end
        if maxHarm >= 3, w(3) = w(3) * 0.45; end
        if maxHarm >= 4, w(4:end) = w(4:end) * 0.3; end

    case 'oboe'
        % Dense harmonic stack with slower decay and odd-harmonic emphasis.
        w = 1 ./ (h .^ 0.55);
        oddIdx = mod(h,2)==1;
        w(oddIdx) = w(oddIdx) * 1.7;
        % Add strong resonance-like bumps (reedy nasal character).
        for center = [5, 9, 14]
            if center <= maxHarm
                idx = max(1,center-1):min(maxHarm,center+1);
                bump = [0.85; 1.7; 0.85];
                w(idx) = w(idx) .* bump(1:numel(idx));
            end
        end

    otherwise
        error('Unknown timbre: %s', timbre);
end
end

function env = raisedCosineEnvelope(nSamples, rampSamples)
if rampSamples*2 > nSamples
    error('Ramp is too long for selected duration.');
end

env = ones(nSamples,1);
if rampSamples > 0
    ramp = 0.5*(1 - cos(pi*(0:rampSamples-1)'/rampSamples));
    env(1:rampSamples) = ramp;
    env(end-rampSamples+1:end) = flipud(ramp);
end
end

function [L, R] = binauralSpatialize(x, fs, azimuthDeg, distanceM)
% Naturalistic ITD/ILD approximation for horizontal plane.
% +135 = left, -135 = right, 0 = center.

c = 343;         % m/s
headRadius = 0.0875; % ~17.5 cm diameter

theta = deg2rad(azimuthDeg);

% Woodworth-style ITD (signed):
% Positive => left ear leads for leftward source.
itd = (headRadius/c) * (theta + sin(theta));

% ILD approximation by azimuth magnitude (frequency-independent).
% Max around 20 dB at extreme lateral positions.
ildDBMax = 20;
ildDB = ildDBMax * sin(abs(theta));

if azimuthDeg > 0
    gL = 1;
    gR = 10^(-ildDB/20);
elseif azimuthDeg < 0
    gR = 1;
    gL = 10^(-ildDB/20);
else
    % Center location: symmetric ears.
    gL = 1;
    gR = 1;
end

% Distance attenuation (inverse law), referenced to 1 m.
distGain = 1 / max(distanceM, eps);
gL = gL * distGain;
gR = gR * distGain;

% Fractional delays for each ear from signed ITD.
if itd >= 0
    dL = 0;
    dR = itd;
else
    dL = -itd;
    dR = 0;
end

L = fractionalDelay(x, dL * fs) * gL;
R = fractionalDelay(x, dR * fs) * gR;

% Keep same length.
L = L(1:numel(x));
R = R(1:numel(x));
end

function y = fractionalDelay(x, delaySamples)
% Linear interpolation fractional-delay line.
n = (0:numel(x)-1)';
y = interp1(n, x, n - delaySamples, 'linear', 0);
end
