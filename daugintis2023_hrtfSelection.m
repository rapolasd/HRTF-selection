function [varargout] = daugintis2023_hrtfSelection(subjects, hrtf_dir, varargin)
%DAUGINTIS2023_HRTFSELECTION select a good and a bad mathcing HRTF based on
%barumerli2023 prediction
%   
%   Requires Auditory Modelling Toolbox (tested with v1.4.1-dev)
%
%
%   The function outputs the error table and saves it as a .csv to be used
%   with the python script of the same name. The python function can be run
%   directly with this function (if python environment is set up) by adding
%   'run_python' flag in the function arguments. However, it has caused
%   matlab to crash before. Alternativelly, the python code can be run
%   within the same folder after this function finishes running.
%   
%   Input parameters:
%       subjects : a cell array of file names within the folder, for which
%                  the good and the bad HRTFs will be selected.
%       hrtf_dir : relative or absolute path to the directory that contains
%                  all HRTFs from which the selection must be made. NB:
%                  subject HRTFs also have to be in this folder.
%
%   Optional flags:
%       'redo_fast'     : runs model prediction with one iteration only and 
%                         doesn't save the error table (for testing purposes)
%       'save_matrices' : explicitly save modelled loacalisation matrices 
%                         as a .mat file. NB: the files can be very big!
%       'run_python'    : run python script directly from the matlab
%                         function. Requires pyenv to be set and can crash
%                         MATLAB. Alternatively, the python script should
%                         be run after the MATLAB script is finished.
%
%   Output:
%       dir_err_t : matlab table, which contains the aggregated errors for
%                   each direction within +-11.5 deg elevation and +-30 deg
%                   azimuth for each subject (template HRTF) with each 
%                   target HRTFs. This tables is also automatically saved
%                   as a dir_err_table.csv.
%       err_dir : original unflattened 4D array that contains the same 
%                 aggregated errors as dir_err_t table.
%       
%
%   Example use:
%       dir_err_t = daugintis2023_hrtfSelection({'P0001_Windowed_48kHz.sofa','P0019_Windowed_48kHz.sofa'},'./test_hrtfs');
%
%       After this finishes, run "python daugintis2023_hrtfSelection.py"
%       from the Terminal (or "python daugintis2023_hrtfSelection.py no_plot" 
%       to prevent saving distribution and selection plots, useful when a
%       lot of subjects are run). The python script will output a
%       selected_HRTFs.csv table, with a good and bad HRTF for each
%       subject. If the target HRTF set is small there might not be a good
%       HRTF selected but that should generally not be the case with a big
%       enough HRTF dataset.



    definput.flags.redo = {'redo','redo_fast'};
    definput.flags.save = {'no_save_matices','save_matrices'};
    definput.flags.python = {'no_run_python','run_python'};
    
    [flags,kv]  = ltfatarghelper({},definput,varargin);
    
    if flags.do_save_matrices
        save_loc_matrices = true;
    else
        save_loc_matrices = false;
    end

    if ~iscell(subjects) && size(subjects,1) == 1
        subjects = {subjects};
    end
    num_subjects = length(subjects);

    % Average model parameters
    sigma_itd = 0.569;
    sigma_ild = 0.75;
    sigma_mon = 4.3;
    sigma_motor = 13.45;
    sigma_prior = 11.5;
    
    % Quick run flag for testing purposes
    if flags.do_redo_fast
        num_exp = 1;
    else
        num_exp = 300;  
    end

    % Extracting the HRTF filenames from the directory
    hrtf_list = dir(fullfile(hrtf_dir,'*.sofa'));    
    num_hrtf = length(hrtf_list);

    
    % Defining limited set of directions at +-30 azimuth +-11.5 elevation
    % (Should be +-30 lateral +-11.5 polar angles but the result is the same)
    dirs = amt_load('barumerli2023','dirs.mat');
    dirs = dirs.cache.value;
    dirs_sph = zeros(size(dirs));
    [dirs_sph(:,1),dirs_sph(:,2),dirs_sph(:,3)] = cart2sph(dirs(:,1), dirs(:,2), dirs(:,3));
    dirs_sph(:,1:2) = rad2deg(dirs_sph(:,1:2));
    dirs_sph_lim = dirs_sph(abs(dirs_sph(:,1)) <= 30 & abs(dirs_sph(:,2)) <= 11.5,:);
    

    clear template_par target_par err_tab err_dir
    amt_disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    amt_disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    amt_disp('Reading SOFA files and extracting features:');

    % Extracting features from the HRTFs
    n_s=0;
    for s = 1:num_hrtf % HRTFs
        sofa_fn = fullfile(hrtf_list(s).folder, hrtf_list(s).name);
        sofa = SOFAload(sofa_fn);
        % Convert to DTF
        sofa = SOFAhrtf2dtf(sofa);
        if any(strcmp(subjects, hrtf_list(s).name))
            amt_disp(sprintf('HRTF (Template and Target): %s', hrtf_list(s).name))
            n_s = n_s+1;
            [template_par(n_s), target_par(s)] = barumerli2023_featureextraction(sofa, ...
                                                    'dtf', ...
                                                    'targ_az', dirs_sph_lim(:,1), ...
                                                    'targ_el', dirs_sph_lim(:,2), ...
                                                    'fs', sofa.Data.SamplingRate);
                                                    
            subject_list{n_s} = hrtf_list(s).name;
        else
            amt_disp(sprintf('HRTF (Target): %s', hrtf_list(s).name)) 
            [~, target_par(s)] = barumerli2023_featureextraction(sofa, 'dtf', ...
                                            'targ_az', dirs_sph_lim(:,1), ...
                                            'targ_el', dirs_sph_lim(:,2), ...
                                            'fs', sofa.Data.SamplingRate);
                                             
%             target_par(s) = barumerli2023_featureextraction(sofa, 'dtf', ...
%                                             'targ_az', dirs_sph_lim(:,1), ...
%                                             'targ_el', dirs_sph_lim(:,2), ...
%                                             'fs', sofa.Data.SamplingRate, ...
%                                             'target');      
        end
    end

    if num_subjects ~= n_s
        warning('subjects list may contain HRTFs, missing from the directory.');
    end
   
    amt_disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    amt_disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    amt_disp('Estimating localisation performance')

    num_dir = size([target_par.monaural],1);
    num_t = num_dir*num_exp;
    if save_loc_matrices
        m = zeros(num_t,8,n_s,num_hrtf);
    end

    % 4D array to store aggregated directional errors
    % number of directions x 6 (target az, el, lat, pol, rmsP, querr) x
    % number of subjects x number of HRTFs
    err_dir = zeros(num_dir, 6, n_s, num_hrtf);
    for s = 1:n_s % looping over subjects
        amt_disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        amt_disp(sprintf('SUBJECT: %s', subject_list{s}));

        template_par_s = template_par(s);
        parfor n = 1:num_hrtf % target HRTFs
            m_temp = barumerli2023('template', template_par_s, ...
                            'target', target_par(n), ...
                            'num_exp', num_exp, ...
                            'sigma_itd', sigma_ild, ...
                            'sigma_ild', sigma_itd, ...
                            'sigma_spectral', sigma_mon, ...
                            'sigma_motor', sigma_motor,...
                            'sigma_prior', sigma_prior,...
                            'MAP');
%             err_tab(s,n) = barumerli2022_metrics(m_temp,'middle_metrics');
            fprintf('Target HRTF %s finished\n', hrtf_list(n).name)

            % Aggregate and store errors for each direction
            for d = 1:num_dir
                [rmsPlocal, querr] = rmsPandQuerr(m_temp(d:num_dir:num_t,:));
                err_dir(d,:, s, n) = [m_temp(d,1), m_temp(d,2), m_temp(d,5), m_temp(d,6), rmsPlocal, querr];
            end
            if save_loc_matrices
                m(:,:,s,n)  = m_temp;
            end
        end

        
    end

    err_dir = err_dir(~isnan(err_dir(:,5,1,1)),:,:,:);
    nd = length(err_dir);
    dir_err_t = array2table(reshape(reshape(permute(err_dir,[1 4 3 2]),[],2),[],size(err_dir,2)), ...
        'VariableNames', {'azi_target', 'ele_target', 'lat_target', 'pol_target', 'rmsP', 'querr'});

    subj_t = array2table([repelem(subjects', num_hrtf*nd, 1), repmat(repelem({hrtf_list.name}', nd, 1), num_subjects, 1)], ...
        'VariableNames', {'template_HRTF', 'target_HRTF'});

    dir_err_t = [subj_t dir_err_t];

    varargout{1} =  dir_err_t;
    varargout{2} =  err_dir;
    
    if ~flags.do_redo_fast
        if flags.do_save_matrices
            save('m.mat', 'm', '-v7.3');
            amt_disp('Localisation matrix saved.');
        end        
        writetable(dir_err_t, 'dir_err_table.csv')
    end     

    if flags.do_run_python
        if pyenv().Version == ""
            warning("Python environment not set up in MATLAB. Please run daugintis2023_hrtfSelection.py manually.")
        else
            pyrunfile("daugintis2023_hrtfSelection.py 'plot'")
        end
    end
    
end

function [rmsPlocal, querr] = rmsPandQuerr(m)
      idx=find(abs(m(:,5))<=30);
      m=m(idx,:);
      idx=find((abs(mynpi2pi(m(:,8)-m(:,6))))<=90);  
      rmsPlocal=sqrt(sum((mynpi2pi(m(idx,8)-m(idx,6))).^2)/size(idx,1));
      querr=100-size(idx,1)/size(m,1)*100;
end

function out_deg=mynpi2pi(ang_deg)
ang=ang_deg/180*pi;
out_rad=sign(ang).*((abs(ang)/pi)-2*ceil(((abs(ang)/pi)-1)/2))*pi;
out_deg=out_rad*180/pi;
end