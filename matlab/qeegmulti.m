function [sf] = qeegmulti(varargin) 

    %% parse input                  
    parser = inputParser;
    parser.addRequired('Input', @iscellstr);
    parser.addParameter('SamplingRate', 500, @isscalar); 
    %parser.addParameter('Band', [1 40], @(x) isvector(x) && length(x)==2); 
    parser.addParameter('Channels', [], @isnumeric); 
    parser.addOptional('Legend', {}, @iscellstr);
    parser.parse(varargin{:});

    % display results of command line parse...
    parser.Results
    disp(['Input filename   : ' parser.Results.Input]);
    disp(['SamplingRate(1/s): ' num2str(parser.Results.SamplingRate)]);
    %disp(['Freq Band(Hz)    : ' mat2str(parser.Results.Band)]);
    disp(['Channels         : ' mat2str(parser.Results.Channels)]);

    l = parser.Results.Legend;  % might be default {}

    % PlotBands is a cell array {m, 2}, where there is one row per band.
    % PlotBands{m, 1} is a string label, used in plots
    % PlotBands{m, 2} is a 2-element array, [band_min band_max]. 
    % All freq bins whose values fall in this range are counted
    %
    % bands are "Theta"6-10Hz, Beta 15-30 LowGamma 30-50
    PlotBands = { 'Theta (6-10)', [6 10] ; 'Beta (15-30)', [15 30] ; 'LowGamma (30-50)', [30 50]};
    bandSums = cell(1, length(parser.Results.Input));


    
    %% Loop over input files. For each, generate plots and spectra.
    
    % cell matrix for spectra results
    % sf{i, 1} is spectra, 1xn. 
    % sf{i, 2} is the frequencies where spectra given, 1xn
    sf=cell(length(parser.Results.Input), 2);

    % cell matrix for power in the bands specified in PlotBands.     
    % row = input file, columns correspond to bands
    bandSums=cell(length(parser.Results.Input), size(PlotBands, 1));    
    for i=[1:length(parser.Results.Input)]
        f = parser.Results.Input(i);
        
        c = parser.Results.Channels(1);
        if length(parser.Results.Channels)>1
            c = parser.Results.Channels(i);
        end
        %[spectra,freqs] = anaEEGFile(f{:}, c);
        [sf{i, 1}, sf{i, 2}] = anaEEGFile(f{:}, c);

        % if legend not set on command line, use filenames
        if any(ismember(parser.UsingDefaults, 'Legend'))
            % parse filename and use basename in legend
            [filepath, filebase, fileext]=fileparts(f{:});    
            l{i}=filebase;
        end

        for j=[1:size(PlotBands, 1)]
            bandSums{i, j} = [bandSums{i, j} ...
                sum(sf{i, 1}(find(sf{i, 2}>=PlotBands{j, 2}(1) & sf{i, 2}<=PlotBands{j, 2}(2)))) ];
        end

        
    end

    % now plot power spectra. While we're at it, generate the summed
    % power over bands for the band plot.

    
    figure('Name', 'Power Spectra');
    for i=[1:length(parser.Results.Input)]
        hold on;
        plot(sf{i, 2}', sf{i, 1});
        
    end
    legend(l);
    xlim([0 50]);

    % generate band plot
    % matrix for plot has the ROWS = each group of bars, one value per
    % input file.
    % Band sum values for a single input file should go along COLUMNS. 
    %
    % note: cell2mat(bandSums) gives a matrix with per-file band sums along
    % the rows. The transpose of that is what we give to the plot function.

    figure('Name', 'Band Power');
    bar(cell2mat(bandSums)');
    set(gca,'xticklabel', {PlotBands{1, 1}, PlotBands{2, 1}, PlotBands{3, 1}});
    leg = legend(parser.Results.Legend);
%    leg.Location='best';
    leg.Location='northwest';
end

function [spectra, freqs] = anaEEGFile(f, channel)

    [filepath, filebase, fileext]=fileparts(f);
    NEdata = [];
    NEdataRefFilPack = [];


%% load data

    disp('Loading data ...');
    if strcmp(fileext, '.easy')
        eeg = pop_easy(f, 0, 0, num2str(channel));
    elseif strcmp(fileext, '.edf')
        error('Cannot load edf files, use easy file.');
    else
      error('Unrecognized file extension')
    end


%% first plot - raw data

    h=figure('Name', filebase);
    subplot(311);
    plot(eeg.times, eeg.data);
    xlabel('Time(ms)');
    ylabel('uV');
    title(['raw(channel ' num2str(channel) ')']);

    subplot(312);
    plot(eeg.times, detrend(eeg.data));
    xlabel('Time(ms)');
    ylabel('uV');
    title(['detrend(channel ' num2str(channel) ')']);

%% spectrum
    [spectra,freqs] = spectopo(eeg.data(1,:), 0, eeg.srate, 'plot', 'off');

%% spectogram

    nfft = 256;
    numoverlap = 128;
    window = hanning(nfft);

    %# spectrogram: make it look like specgram
    subplot(313);
    [S,F,T,P] = spectrogram(eeg.data(1, :), window, numoverlap, nfft, eeg.srate);
    imagesc(T, F, 20*log10(P));
    axis xy, colormap(jet), ylabel('Frequency');
    title('Spectrogram');

return;
end