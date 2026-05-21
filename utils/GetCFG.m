function cfg = GetCFG(phantom_case)

datapath = '/home/trulssv/Documents/IodineGadolinium/md_sim';


% --------------------------

% Note that Linux version cannot handle .pp phantoms.
% Convertion to .ppm does not work.
cfg.seed = 42; % set seed to get same noise for all phantoms


switch phantom_case
    case 1
        cal_date = '241121_155419'; % data of calibration acquisitions
        calpath = fullfile(datapath, ['calibration_' cal_date]);
        cal_baseline = fullfile(calpath,['cal_baseline_' cal_date '.mat']);
        load(cal_baseline, 'cfg');
    case 2
        ca.l_date = '241121_155419';
    case 3
        cal_date = '241120_170809';
        calpath = fullfile(datapath, ['calibration_' cal_date]);
        cal_baseline = fullfile(calpath,['cal_baseline_' cal_date '.mat']);
        load(cal_baseline, 'cfg');
        cfg.recon_fov = 500;
    case 4
        data_path = '/media/data/trulssv/Antropomorphic_Phantoms/calibration_250326_155807/phantoms';
        cal_dir = 'vox_case_00100';
        cal_filename = 'recon_for_250329_110221_250329_110708.mat';
        
        cal_baseline = fullfile(data_path, cal_dir, cal_filename);
        load(cal_baseline, 'cfg');
    case 5 
        data_path = '/media/data/trulssv/Antropomorphic_Phantoms/calibration_250326_155807/phantoms';
        cal_dir = 'vox_case_00100';
        cal_filename = 'recon_for_250329_110221_250329_110708.mat';
        
        cal_baseline = fullfile(data_path, cal_dir, cal_filename);
        load(cal_baseline, 'cfg');
           case 9
        data_path = '/media/data/trulssv/Antropomorphic_Phantoms/calibration_250326_155807/phantoms';
        cal_dir = 'vox_case_00100';
        cal_filename = 'recon_for_250329_110221_250329_110708.mat';
        
        cal_baseline = fullfile(data_path, cal_dir, cal_filename);
        load(cal_baseline, 'cfg');
end



% Add save_transmission flag to cfg struct
if exist('save_transmission','var') 
    cfg.this_is_a_transmissionscan = save_transmission;
    cfg.copy_rows = 1; % don't use MD anyway when saving transmission/pathlength
    clear save_transmission
end    


end