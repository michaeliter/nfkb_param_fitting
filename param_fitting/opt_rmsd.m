function rmsd = optimize(p, sheet, row_num)
%% OPTIMIZE Simulate with passed parameters
% Set up simulation
    warning('off');
    lb = [0.01775, 0.1*log(2)/9, 0.05188, 1E-7, 0.0004, 7, 0.04];
    ub = [1.775, 10*log(2)/9, 5.188, 1E-5, 0.04, 28, 0.3];
    ref = [0.1775, log(2)/9, 0.5188, 1E-6, 0.004, 14, 0.08];
    params = constrain(p,lb,ub);
    p_mod = [6 3 params(1); 15 1 params(2); 16 1 params(2); 66 1 params(3); 68 1 params(4); 75 1 params(5); 6 4 params(6)];
    init_mod = {'NFkB', params(7)};
    doses = [100];
    dose_scale = 1/1500; % Convert to uM. Pam3CSK molecular weight: 1.5KDa
    names = {'NFkBn', 'IkBaNFkBn'};
    options = struct;
    options.DEBUG = 0;
    options.SIM_TIME = 8*60 + 5;
    % Simulate (only need to equilibrate on first iteration)
    output = [];
    for i = 1:length(doses)
        if isempty(output)
            try
                [t,x,simdata] = nfkbSimulate({'Pam3CSK',doses(i)*dose_scale},names, p_mod, init_mod,options);
            catch
                warning('on');
                warning('nfkbSimulate failed. Proceeding to set rmsd to inf.');
                warning('off');
                rmsd = inf;
                return
            end
        else
            options.STEADY_STATE = simdata.STEADY_STATE;
            [~,x] = nfkbSimulate({'Pam3CSK',doses(i)*dose_scale},names, p_mod, {},options);
        end
        output = cat(3,output,x);
    end
    sim_nfkb_curve = squeeze(output(:,strcmp(names,'NFkBn'),:));
    sim_ikba_curve = squeeze(output(:,strcmp(names,'IkBaNFkBn'),:));
    sim_nfkb_curve = sim_nfkb_curve + sim_ikba_curve;
    sim_nfkb_curve = sim_nfkb_curve(1 : 5 : end);
    % check for NaN values in sim
    if any(isnan(sim_nfkb_curve),'all')
        warning('on');
        warning('NaN values detected in simulation. Proceeding to set rmsd to inf.');
        warning('off');
        rmsd = inf;
        return
    end
%% Get experimental data
    row = 'A' + string(row_num) + ':CT' + string(row_num);
    experimental = readmatrix(sheet, 'Range', row).';
    experimental = experimental.*0.0313;
    % bring sim to same zero point as experimental data
    zero = sim_nfkb_curve(1,1) - experimental(1,1);
    sim_nfkb_curve = sim_nfkb_curve - zero;
    dif = sim_nfkb_curve - experimental;
    dif(1:48,1) = dif(1:48,1).*sqrt(2);
    ratio = params./ref;
    ratio(1,1:5) = log10(ratio(1,1:5));
    ratio(1,6:7) = log2(ratio(1,6:7));
    rmsd = rms(dif, 'omitnan') + 0.1*rms(ratio);
end
