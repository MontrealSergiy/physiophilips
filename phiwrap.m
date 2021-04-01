function phiwrap(save_dir, varargin)   
%% A command line wrapper for the main entry function 
% tapas_physio_main_create_regressors of the Physio (for Philips equipment)
% This version target mostly Philips equipment 
% and enables compilation and cbrain integration of physio.
% While few tool parameters such as input files becames command line
% parameters, the bulk of one hunderd of parameters and options are to be
% provided in json config file. User may dump his existing matlab json 
% physio structure into such file. Parameters can be set from
% an exisiting example of config, provided by toolset. The command line
% parameters have highest priority, than ones from the json file, and
% 
% for documentation of the parameters, see also tapas_physio_new (e.g., via
% edit tapas_physio_new). However parameters related to 
%
%
%
% NOTE: All inputs in physio-structure have to be specified previous to
%       running this function.
%
% IN
%   save_dir         the folder where results will be stored, the only
%                    positional parameter
%   physlogfile      the file with physio logs
%   fmrimovefile     optional physio log file (Philips uses same file for all data)
%   presetdefaults   optional preset defaults choice, takes values  
%                    'Philips_ECG3T_v1', in the latter 
%                    case the parameters will get default values from a corresponding 
%                    example script of the toolset;
%   configfile       optional custom file with a json dump of misc parameters, it does
%                    overwrite the positional arguments. Note that input or 
%                    output files occur in the config dump, thouse will be ignored.
%                    either preset defaults or configfile are recommended
%                    otherwise execution will fail
%
% EXAMPLES
%
%   phiwrap('results_dir',  ...
%           'physlogfile',  'SCANPHYSLOG.log', ...
%           'fmrimovefile', 'rp_fMRI.txt', ...
%           'configfile',   'myparams.json') 
%
% REFERENCES
%
% CBRAIN        www.cbrain.com 
% RETROICOR     regressor creation based on Glover et al. 2000, MRM 44 and
%               Josephs et al. 1997, ISMRM 5, p. 1682
%               default model order based on Harvey et al. 2008, JMRI 28
% RVT           (Respiratory volume per time) Birn et al. 2008, NI 40
% HRV           (Heart-rate  variability) regressor creation based on
%               Chang et al2009, NI 44
%
% See also tapas_physio_new

% Author:    Serge Boroday
% Created:   2021-03-16
% Copyright: McGill University
%
% The original tool is by Institute for Biomedical Engineering, 
%               University of Zurich and ETH Zurich.
%
% This file is a wrapper for TAPAS PhysIO Toolbox, Phillips equipment only


%% Parse arguments

p = inputParser;
validDefaultConfig = @(x) ischar(x) && (strcmpi(x,'PHILIPS_ECG3T_V1') || strcmpi(x,'standard'));
addRequired(p, 'save_dir', @isvarname);  % prevent polluting parent folders 
addParameter(p,'physlogfile', '', @isfile);   % all in file for Philips
addParameter(p,'configfile',  '', @isfile);
addParameter(p,'fmrimovefile','', @isfile)
addParameter(p,'presetdefaults', 'standard', validDefaultConfig)
parse(p, save_dir, varargin{:});

physlogfile    = p.Results.physlogfile;
configfile     = p.Results.configfile;
presetdefaults = p.Results.presetdefaults;
fmrimovefile   = p.Results.fmrimovefile;
% hasmovefile    = p.Results.fmrimovefile;

%% Create default parameter structure with all fields
physio = tapas_physio_new();


%% Individual Parameter settings. Modify to your need and remove default settings
physio.save_dir = save_dir;
physio.log_files.vendor = 'Philips'; 


physio = load_params(presetdefaults, physio);

if fmrimovefile
            physio.model.movement.include = true;
            physio.model.movement.file_realignment_parameters = {fmrimovefile};
else
            physio.model.movement.include = false;
end    

% read bulk of params from json file
if configfile
    uconfig = jsondecode(fileread(configfile));
    physio = update_physio(physio, uconfig);
    
end





 


% input params
physio.log_files.cardiac = {physlogfile};
physio.log_files.respiration = {physlogfile};
physio.log_files.scan_timing = {physlogfile};

% update results folder
physio.save_dir = save_dir;

mkdir save_dir
% save final input params
fid = fopen(strcat(save_dir, '/physio_params_in.json'), 'wt');
fprintf(fid, '%s', jsonencode(physio));
fclose(fid);



% postpone figure generation in first run - helps with compilation
if isfield(physio, 'verbose') && isfield(physio.verbose, 'level')
     verbose_level = physio.verbose.level;
     physio.verbose.level = 0;
     if isfield(physio, 'fig_output_file')
         fig_output_file = physio.verbose.fig_output_file;
     else
         fig_output_file = 'PhysIO_output.png'; 
     end    
else
  verbose_level = 0;
end 


%% Run physiological recording preprocessing and noise modeling

physio = tapas_physio_main_create_regressors(physio);


%% Build figures
if verbose_level
  physio.verbose.fig_output_file = fig_output_file; % has to reset, the old value is distorted
  physio.verbose.level = verbose_level;
  tapas_physio_review(physio);
end


end

function [physio] = update_physio(physio_1, physio_2)
%% update major fields of nested physio config struct
% only these fields and their subfields are overwriten

f = {'preproc', 'scan_timing', 'model', 'verbose', 'on_secs'}; 

for i = 1:length(f)
    fieldname = f{i};
    if isfield(physio_1, fieldname) && isfield(physio_2, fieldname)
        physio_1.(fieldname) = merge_struct(physio_1.(fieldname), physio_2.(fieldname));
    end     
end
physio = physio_1;
end

function [S] = merge_struct(S_1, S_2)
% update the first struct with values and keys of the second and returns the result
% deep update, merges substructrues recursively, the values from the first
% coinside

f = fieldnames(S_2);

for i = 1:length(f)
    if isfield(S_1, f{i}) && isstruct(S_1.(f{i})) && isstruct(S_2.(f{i}))
        S_1.(f{i}) = merge_struct(S_1.(f{i}), S_2.(f{i}));
    else   
        S_1.(f{i}) = S_2.(f{i});
    end        
end
S = S_1;
end

function [physio] = load_params(byexample, physio)
    
    if strcmpi(byexample, "PHILIPS_ECG3T_V1")
        %cut and paste of an example from https://www.tapas.tnu-zurich.com/examples_v4.0.0.zip, 
        %save initalization and execution
     
        physio.model.movement.include = true;
        physio.model.movement.file_realignment_parameters = {'rp_fMRI.txt'};
        physio.model.movement.order = 6;      
        physio.model.movement.censoring_threshold = [3 Inf];
        physio.model.movement.censoring_method = 'MAXVAL';

        physio.log_files.vendor = 'Philips';
        physio.log_files.relative_start_acquisition = 0;
        physio.log_files.align_scan = 'last';
        physio.scan_timing.sqpar.Nslices = 37;
        physio.scan_timing.sqpar.TR = 2.5;
        physio.scan_timing.sqpar.Ndummies = 3;
        physio.scan_timing.sqpar.Nscans = 495;
        physio.scan_timing.sqpar.onset_slice = 19;
        physio.scan_timing.sync.method = 'gradient_log';
        physio.scan_timing.sync.grad_direction = 'y';
        physio.scan_timing.sync.zero = 0.4;
        physio.scan_timing.sync.slice = 0.45;
        physio.preproc.cardiac.modality = 'ECG';
        physio.preproc.cardiac.filter.include = false;
        physio.preproc.cardiac.filter.type = 'butter';
        physio.preproc.cardiac.filter.passband = [0.3 9];
        physio.preproc.cardiac.initial_cpulse_select.method = 'load_from_logfile';
        physio.preproc.cardiac.initial_cpulse_select.max_heart_rate_bpm = 90;
        physio.preproc.cardiac.initial_cpulse_select.file = 'initial_cpulse_kRpeakfile.mat';
        physio.preproc.cardiac.initial_cpulse_select.min = 0.4;
        physio.preproc.cardiac.posthoc_cpulse_select.method = 'off';
        physio.preproc.cardiac.posthoc_cpulse_select.percentile = 80;
        physio.preproc.cardiac.posthoc_cpulse_select.upper_thresh = 60;
        physio.preproc.cardiac.posthoc_cpulse_select.lower_thresh = 60;
        physio.model.orthogonalise = 'none';
        physio.model.censor_unreliable_recording_intervals = false;
        physio.model.output_multiple_regressors = 'multiple_regressors.txt';
        physio.model.output_physio = 'physio.mat';
        physio.model.retroicor.include = true;
        physio.model.retroicor.order.c = 3;
        physio.model.retroicor.order.r = 4;
        physio.model.retroicor.order.cr = 1;
        physio.model.rvt.include = false;
        physio.model.rvt.delays = 0;
        physio.model.hrv.include = false;
        physio.model.hrv.delays = 0;
        physio.model.noise_rois.include = false;
        physio.model.noise_rois.thresholds = 0.9;
        physio.model.noise_rois.n_voxel_crop = 0;
        physio.model.noise_rois.n_components = 1;
        physio.model.noise_rois.force_coregister = 1;

        physio.model.other.include = false;
        physio.verbose.level = 2;
        physio.verbose.process_log = cell(0, 1);
        physio.verbose.fig_handles = zeros(1, 0);
        physio.verbose.fig_output_file = 'PhysIO_output.png';
        physio.verbose.show_figs = false;
        physio.verbose.save_figs = true;
        physio.verbose.close_figs = true;
        physio.ons_secs.c_scaling = 1;
        physio.ons_secs.r_scaling = 1;
           
    end
end