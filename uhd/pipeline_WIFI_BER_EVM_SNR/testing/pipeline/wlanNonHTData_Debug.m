function [y,afterMap] = wlanNonHTData_Debug(PSDU,cfgNonHT,varargin)
%wlanNonHTData Non-HT Data field processing of the PSDU
%
%   Y = wlanNonHTData(PSDU,CFGNONHT) generates the non-HT format Data
%   field time-domain waveform for the input PLCP Service Data Unit (PSDU). 
%
%   Y is the time-domain non-HT Data field signal. It is a complex matrix
%   of size Ns-by-Nt, where Ns represents the number of time-domain samples
%   and Nt represents the number of transmit antennas.
%
%   PSDU is the PHY service data unit input to the PHY. It is a double
%   or int8 typed column vector of length CFGNONHT.PSDULength*8, with each
%   element representing a bit.
%
%   CFGNONHT is the format configuration object of type <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a> which
%   specifies the parameters for the non-HT format. Only OFDM modulation
%   type is supported.
%
%   Y = wlanNonHTData(...,SCRAMINIT) optionally allows specification of the
%   scrambler initial state, or initial pseudorandom scrambler sequence
%   bits SCRAMINIT.
%
%   When bandwidth signaling is not used (CFGNONHT.SignalChannelBandwidth
%   is false) SCRAMINIT is the initial state of the scrambler. SCRAMINIT
%   must be a scalar between 1 and 127 inclusive, or a corresponding double
%   or int8-typed binary 7-by-1 column vector. When not specified, 93 is
%   used.
%
%   When bandwidth signaling is used (CFGNONHT.SignalChannelBandwidth is
%   true), SCRAMINIT is the initial pseudorandom scrambler sequence as
%   described in IEEE 802.11-2016 Table 17-7. The valid range depends on
%   the value of CFGNONHT.BandwidthOperation and CFGNONHT.ChannelBandwidth.
%   SCRAMINIT must be a scalar between MIN and MAX inclusive or a double or
%   int8-typed binary column vector of length NB. When not specified,
%   SCRAMINIT defaults to the NB most significant bits of 93. The values of
%   MIN, MAX, and NB are specified in the table below.
%   
%   |  BandwidthOperation   | ChannelBandwidth | MAX | MIN | NB |
%   |       'Absent'        |      'CBW20'     | 31  |  1  | 5  |
%   |       'Absent'        |         -        | 31  |  0  | 5  |
%   | 'Static' or 'Dynamic' |      'CBW20'     | 15  |  1  | 4  |
%   | 'Static' or 'Dynamic' |         -        | 15  |  0  | 4  |
%
%   Y = wlanNonHTData(...,'OversamplingFactor',OSF) generates the
%   NonHT-Data oversampled by a factor OSF. OSF must be >=1. The resultant
%   cyclic prefix length in samples must be integer-valued for all symbols.
%   The default is 1.
%
%   % Example 1:
%   %  Generate the signal for a 20 MHz non-HT OFDM data field for 36 Mbps.
%
%     cfgNonHT = wlanNonHTConfig('MCS',5);                % Configuration
%     psdu = randi([0 1],cfgNonHT.PSDULength*8,1,'int8'); % PSDU in bits
%     tx = wlanNonHTData(psdu,cfgNonHT);
%
%   % Example 2:
%   %  Generate signal for 80 MHz non-HT duplicate data field with random
%   %  PSDU and scrambler initial sequence bits.
%
%     cfgNonHT = wlanNonHTConfig;
%     cfgNonHT.ChannelBandwidth = 'CBW80';
%     cfgNonHT.SignalChannelBandwidth = true;
%     cfgNonHT.BandwidthOperation = 'Static';
%     
%     % Random PSDU bits
%     psdu = randi([0 1],cfgNonHT.PSDULength*8,1,'int8');
%     % Random initial scrambler sequence bits
%     scramInit = randi(scramblerRange(cfgNonHT));
%     % Data field
%     tx = wlanNonHTData(psdu,cfgNonHT,scramInit);
%   
%   See also wlanNonHTConfig, wlanNonHTOFDMDemodulate,
%   wlanNonHTDataBitRecover, wlanInterpretScramblerState,
%   wlanWaveformGenerator, wirelessWaveformGenerator.

%   Copyright 2015-2021 The MathWorks, Inc.

%#codegen

narginchk(2,5);

% Validate inputs
% Validate the format configuration object
validateattributes(cfgNonHT, {'wlanNonHTConfig'}, ...
    {'scalar'}, mfilename, 'format configuration object');
% Only applicable for OFDM and DUP-OFDM modulations
coder.internal.errorIf( ~strcmp(cfgNonHT.Modulation, 'OFDM'), ...
                        'wlan:wlanNonHTData:InvalidModulation');
s = validateConfig(cfgNonHT); 
validateattributes(PSDU, {'double', 'int8'},...
    {'real', 'binary', 'size', [cfgNonHT.PSDULength*8 1]}, ...
    mfilename, 'PSDU input');
[scramInitBits,osf] = processVarargin(cfgNonHT,varargin{:});

% Determine number of symbols and pad length
numSym = s.NumDataSymbols;
numPad = s.NumPadBits;

mcsTable = wlan.internal.getRateTable(cfgNonHT);
Nservice = 16;
Ntail = 6;

% Scramble padded data
% [service; psdu; tail; pad] processing
paddedData = [zeros(Nservice,1,'int8'); PSDU; zeros(Ntail,1,'int8'); zeros(numPad, 1,'int8')];
scrambData = wlanScramble(paddedData, scramInitBits);

% Zero-out the tail bits again for encoding
scrambData(16+length(PSDU) + (1:Ntail)) = zeros(Ntail,1);

% BCC Encoding
encodedData = wlanBCCEncode(scrambData, mcsTable.Rate);

% BCC Interleaving
interleavedData = wlanBCCInterleave(encodedData, 'Non-HT', mcsTable.NCBPS);

% Constellation mapping
mappedData = wlanConstellationMap(interleavedData, mcsTable.NBPSCS);
afterMap = mappedData;
% Non-HT pilots, from IEEE Std 802.11-2012, Eqn 18-22
% Reshape to form OFDM symbols
mappedData = reshape(mappedData, mcsTable.NCBPS/mcsTable.NBPSCS, numSym);

% Non-HT pilots, from IEEE Std 802.11-2012, Eqn 18-22
z = 1; % Offset by 1 to account for HT-SIG pilot symbol
pilotValues = wlan.internal.nonHTPilots(numSym,z);
        
% Data packing with pilot insertion, replicate over 20 MHz subchannels
ofdm = wlan.internal.vhtOFDMInfo('NonHT-Data', cfgNonHT.ChannelBandwidth, 1);
packedData = complex(zeros(ofdm.FFTLength, numSym));
packedData(ofdm.ActiveFFTIndices(ofdm.DataIndices), :) = repmat(mappedData, ofdm.NumSubchannels, 1);
packedData(ofdm.ActiveFFTIndices(ofdm.PilotIndices), :) = repmat(pilotValues, ofdm.NumSubchannels, 1);

% Apply gamma rotation, replicate over antennas and apply cyclic shifts
[data,scalingFactor] = wlan.internal.legacyFieldMap(packedData, ofdm.NumTones, cfgNonHT);

% OFDM modulate
y = wlan.internal.wlanOFDMModulate(data, ofdm.CPLength, osf)*scalingFactor;

end

% Validate and return scrambler initialization and get OSF
function [scramInitBits,osf] = processVarargin(cfgNonHT,varargin)
    % Default options
    osf = 1;
    % Use most significant bits of default as initial scrambler sequence
    userScramInitBits = int8([1; 0; 1; 1; 1; 0; 1]); % Default is 93 
    bandwidthSignaling = any(strcmp(cfgNonHT.ChannelBandwidth,{'CBW20','CBW40','CBW80','CBW160'})) && cfgNonHT.SignalChannelBandwidth;

    % Validate scrambler initialization
    [range,numReqScramBits] = scramblerRange(cfgNonHT);
    minVal = range(1);
    maxVal = range(2);
    if nargin>1
        if ~(ischar(varargin{1}) || isstring(varargin{1}))
            scramInit = varargin{1};
            % Validate scrambler init
            validateattributes(scramInit,{'double','int8'},{'real','integer','nonempty'},mfilename,'Scrambler initialization');
    
            if isscalar(scramInit)
                % Check for correct range
                scramInitErrorIf(any((scramInit<minVal) | (scramInit>maxVal)));
                userScramInitBits = int8(int2bit(scramInit,numReqScramBits));
            else
                % Check for binary vector of correct size and range
                scramInitErrorIf(~iscolumn(scramInit) || numel(scramInit)~=numReqScramBits  || ...
                  any((scramInit~=0) & (scramInit~=1)) || all(bit2int(double(scramInit(:)),numel(scramInit))<minVal)); % Cast to double for codegen
                userScramInitBits = int8(scramInit);
            end
            osf = wlan.internal.parseOSF(varargin{2:end});
        else
            osf = wlan.internal.parseOSF(varargin{:});
        end
    end
    if bandwidthSignaling
        scramInitBits = bandwidthSignalScramblerInit(userScramInitBits(1:numReqScramBits,1),cfgNonHT.ChannelBandwidth,cfgNonHT.BandwidthOperation);
    else
        scramInitBits = userScramInitBits(1:numReqScramBits,1); % Index for codegen
    end
    
    function scramInitErrorIf(condition)
        if bandwidthSignaling
            coder.internal.errorIf(condition,'wlan:wlanNonHTData:InvalidScramInitBWSignaling',minVal,maxVal,numReqScramBits);
        else
            coder.internal.errorIf(condition,'wlan:wlanNonHTData:InvalidScramInit',minVal,maxVal,numReqScramBits);
        end
    end
end

% Get scrambler initial bits
function scramInitBits = bandwidthSignalScramblerInit(userScramInitBits,chanBW,bandwidthOperation)
    scramInitBits = zeros(7,1,'int8');
    numReqScramBits = size(userScramInitBits,1);

    % Signal channel bandwidth - IEEE 802.11-2016 Table 17-7, 17-8, and 17-10
    switch chanBW
      case 'CBW20'
          bwBits = [0; 0]; % [B5 B6]
      case 'CBW40'
          bwBits = [1; 0]; % [B5 B6]
      case 'CBW80'
          bwBits = [0; 1]; % [B5 B6]
      otherwise % 'CBW160'
          bwBits = [1; 1]; % [B5 B6]
    end
    scramInitBits(6:7) = bwBits;

    % Signal dynamic bandwidth capability
    switch bandwidthOperation
        case 'Static'
            scramInitBits(5) = 0;
        case 'Dynamic'
            scramInitBits(5) = 1;
      otherwise % 'Absent'
            assert(numReqScramBits==5);
    end
    
    % Use only required number of pseudorandom bits
    scramInitBits(1:numReqScramBits) = userScramInitBits(:,1); % Index for codegen

    % Create scrambler initialization from initial scrambler sequence
    % (equivalent of feeding bits into scrambler then flicking the switch)
    scramInitBits = wlan.internal.scramblerInitialState(scramInitBits);
end