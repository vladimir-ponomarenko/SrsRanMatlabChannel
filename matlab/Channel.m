
clear global; clear classes; close all; clc;
fprintf('Starting MATLAB script v1.0.0 ...\n');
warning('on','all');

jeromq_jar_path = '/home/user/jeromq/target/jeromq-0.6.0.jar'; 
quadriga_base_path = '/home/user/Quadriga/quadriga_src'; 

if ~isfile(jeromq_jar_path), error('jeromq JAR file not found: %s', jeromq_jar_path); end
if ~isfolder(quadriga_base_path), error('Quadriga source folder not found: %s', quadriga_base_path); end

javaaddpath(jeromq_jar_path);
try
    addpath(genpath(quadriga_base_path));
    rehash toolboxcache;
    if ~exist('qd_layout', 'class'), error('QuaDRiGa class "qd_layout" not found.'); end
    disp('QuaDRiGa paths added successfully.');
catch ME_path
    error('Error adding QuaDRiGa path: %s', ME_path.message);
end

try
    import org.zeromq.ZMQ.*;
    import org.zeromq.*;
    disp('ZeroMQ Java classes imported successfully.');
catch ME_zmq_import
     error('Failed to import ZeroMQ classes. Check JAR path "%s". Error: %s', jeromq_jar_path, ME_zmq_import.message);
end

port_api = 2111;
context = [];
socket_api_proxy = [];
recv_timeout_ms = 200;

try
    context = ZMQ.context(1);
    socket_api_proxy = context.socket(ZMQ.REP);
    socket_api_proxy.setLinger(0);
    socket_api_proxy.setReceiveTimeOut(recv_timeout_ms);
    socket_api_proxy.setRcvHWM(2000);
    socket_api_proxy.setSndHWM(2000);
    socket_api_proxy.bind(sprintf('tcp://*:%d', port_api));
    fprintf('ZMQ REP server bound to tcp://*:%d (Recv Timeout: %d ms)\n', port_api, recv_timeout_ms);
catch ME_zmq_bind
    fprintf(2,'FATAL ERROR binding ZMQ socket: %s\n', ME_zmq_bind.message);
    if ~isempty(context) && ismethod(context,'term'), try context.term(); catch, end, end
    return;
end

global pauseFlag stopFlag updateChannelFlag channelParams;
pauseFlag = false;
stopFlag = false;
updateChannelFlag = false;
channelParams = struct();

disp('Setting up simulation parameters...');

channelParams.modelingMode = 'Manual'; 

channelParams.center_frequency = 2.56e9;
channelParams.scenario = '3GPP_38.901_UMa_NLOS';
channelParams.tx_pos = [0; 0; 25];
channelParams.rx_pos = [200; 0; 1.5];
channelParams.numClusters = 16;
channelParams.bandwidth = 23.04e6;
channelParams.N_fft_qdrg_calc = 131072; 
fprintf('Setting QuaDRiGa channel calculation N_fft (N_fft_qdrg_calc) = %d\n', channelParams.N_fft_qdrg_calc);
channelParams.H_freq_normalized = complex(ones(channelParams.N_fft_qdrg_calc, 1));

channelParams.manualPL_model = '3GPP_UMa_NLOS';
channelParams.h_tx = channelParams.tx_pos(3);
channelParams.h_rx = channelParams.rx_pos(3);
channelParams.current_pl_dB = NaN;
channelParams.pg_lin_sqrt_applied = 1.0;

noise_power_dBm = -190;
noise_power_watts = 10^((noise_power_dBm - 30) / 10);
noise_sigma = sqrt(noise_power_watts / 2);
fprintf('AWGN settings: Noise Power = %.1f dBm, Sigma = %.2e\n', noise_power_dBm, noise_sigma);

disp('Calculating initial channel model...');
updateChannelModel(true);

if isempty(channelParams.H_freq_normalized) || isnan(channelParams.pg_lin_sqrt_applied)
    fprintf(2,'Initial channel model calculation failed. Exiting.\n');
    cleanup(socket_api_proxy, context, []);
    return;
end
disp('Initial channel model calculation complete.');

fig = figure('Name', 'QuaDRiGa/Manual PL Control (v1.0.0)', ...
             'NumberTitle', 'off', ...
             'Position', [100 50 550 450], ...
             'CloseRequestFcn', @closeFigureCallback, ...
             'Visible', 'on');

y_pos = 400; spacing = 35; label_width = 150; edit_width = 100; popup_width = 240;

uicontrol('Style', 'pushbutton', 'String', 'Pause/Resume', 'Position', [20, y_pos, 100, 30], 'Callback', @togglePause);
uicontrol('Style', 'pushbutton', 'String', 'STOP Script', 'Position', [130, y_pos, 100, 30], 'Callback', @stopScript);

y_pos = y_pos - spacing - 5;
uicontrol('Style','text','String','Modeling Mode:','Position',[20 y_pos label_width 20],'HorizontalAlignment','left');
modePopup = uicontrol('Style','popupmenu','String',{'Manual PL + QDR SSF','Full QuaDRiGa'}, ...
                      'Position',[label_width+30 y_pos popup_width-40 25], ...
                      'Callback', @modeChangedCallback);

y_pos = y_pos - spacing;
uicontrol('Style','text','String','2D Distance (m):','Position',[20 y_pos label_width-50 20],'HorizontalAlignment','left');
distEdit = uicontrol('Style','edit','String',sprintf('%.1f', norm(channelParams.tx_pos(1:2) - channelParams.rx_pos(1:2))),'Position',[label_width+30 y_pos edit_width 25]);

y_pos = y_pos - spacing;
uicontrol('Style','text','String','QDR Scenario:','Position',[20 y_pos label_width 20],'HorizontalAlignment','left');

scenarioList = { ...
    '3GPP_38.901_UMi_NLOS','3GPP_38.901_UMi_LOS', ...
    '3GPP_38.901_UMa_NLOS','3GPP_38.901_UMa_LOS', ...
    '3GPP_38.901_RMa_NLOS','3GPP_38.901_RMa_LOS', ... 
    '3GPP_38.901_Indoor_Mixed_Office_NLOS', '3GPP_38.901_Indoor_Mixed_Office_LOS',... 
    '3GPP_38.901_Indoor_Open_Office_NLOS', '3GPP_38.901_Indoor_Open_Office_LOS', ... 
    '3GPP_38.901_InF_SL_NLOS', '3GPP_38.901_InF_SL_LOS', ...  
    '3GPP_38.901_InF_DL_NLOS', '3GPP_38.901_InF_DL_LOS', ...  
    '3GPP_38.901_InF_SH_NLOS', '3GPP_38.901_InF_SH_LOS', ...  
    '3GPP_38.901_InF_DH_NLOS', '3GPP_38.901_InF_DH_LOS', ...  
    '5G-ALLSTAR_DenseUrban_NLOS','5G-ALLSTAR_DenseUrban_LOS', ... 
    '5G-ALLSTAR_Urban_NLOS','5G-ALLSTAR_Urban_LOS', ... 
    '5G-ALLSTAR_Suburban_NLOS','5G-ALLSTAR_Suburban_LOS', ... 
    '5G-ALLSTAR_Rural_NLOS','5G-ALLSTAR_Rural_LOS', ... 
    'mmMAGIC_UMi_NLOS','mmMAGIC_UMi_LOS', ... 
    'mmMAGIC_Indoor_NLOS','mmMAGIC_Indoor_LOS', ... 
    'QuaDRiGa_Industrial_NLOS','QuaDRiGa_Industrial_LOS', ... 
    'QuaDRiGa_UD2D_NLOS','QuaDRiGa_UD2D_LOS', ... 
    'LOSonly', 'None (Flat SSF)'};
scenarioPopup = uicontrol('Style','popupmenu','String',scenarioList,'Position',[label_width+30 y_pos popup_width 25]);
initialScenarioIndex = find(strcmp(scenarioPopup.String, channelParams.scenario), 1);
if isempty(initialScenarioIndex)
    warning('Initial scenario "%s" not found in list. Defaulting to first item.', channelParams.scenario);
    initialScenarioIndex = 1;
    channelParams.scenario = scenarioList{1};
end
set(scenarioPopup, 'Value', initialScenarioIndex);

y_pos = y_pos - spacing;
uicontrol('Style','text','String','#Clusters (QDR):','Position',[20 y_pos label_width 20],'HorizontalAlignment','left');
clustersEdit = uicontrol('Style','edit','String',num2str(channelParams.numClusters),'Position',[label_width+30 y_pos edit_width 25]);

% --- Элементы для Ручного Режима ---
y_pos = y_pos - spacing - 5;
manualPLText = uicontrol('Style','text','String','Manual PL Model:','Position',[20 y_pos label_width 20],'HorizontalAlignment','left');
manualPLList = {'3GPP_UMa_NLOS', '3GPP_UMa_LOS', ...
                  '3GPP_UMi_NLOS', '3GPP_UMi_LOS', ...
                  '3GPP_RMa_NLOS', '3GPP_RMa_LOS', ...
                  'FSPL'};
manualPLPopup = uicontrol('Style','popupmenu','String',manualPLList,'Position',[label_width+30 y_pos popup_width-60 25]);
initialPLIndex = find(strcmp(manualPLPopup.String, channelParams.manualPL_model), 1);
if isempty(initialPLIndex), initialPLIndex = 1; end
set(manualPLPopup, 'Value', initialPLIndex);

y_pos = y_pos - spacing;
txHeightText = uicontrol('Style','text','String','Tx Height (m):','Position',[20 y_pos label_width-50 20],'HorizontalAlignment','left');
txHeightEdit = uicontrol('Style','edit','String',sprintf('%.1f', channelParams.h_tx),'Position',[label_width-40 y_pos edit_width-20 25]);
rxHeightText = uicontrol('Style','text','String','Rx Height (m):','Position',[label_width+edit_width-20 y_pos label_width-50 20],'HorizontalAlignment','right');
rxHeightEdit = uicontrol('Style','edit','String',sprintf('%.1f', channelParams.h_rx),'Position',[label_width+edit_width+label_width-60 y_pos edit_width-20 25]);


y_pos = y_pos - spacing - 10;
applyButton = uicontrol('Style', 'pushbutton', 'String', 'Apply & Update Channel', 'Position', [20, y_pos, 220, 30], ...
          'Callback', {@applyParametersCallback, modePopup, distEdit, scenarioPopup, clustersEdit, manualPLPopup, txHeightEdit, rxHeightEdit});

y_pos = y_pos - spacing - 5;
statusText = uicontrol('Style', 'text', 'String', 'Status: Running', 'Position', [20, y_pos, 510, 30], 'HorizontalAlignment', 'left','FontSize',10);

y_pos = y_pos - spacing + 10;
currentGainText = uicontrol('Style', 'text', ...
                          'String', sprintf('PL: %.2f dB | Sqrt Gain: %.3e', channelParams.current_pl_dB, channelParams.pg_lin_sqrt_applied), ...
                          'Position', [20, y_pos, 510, 25], 'HorizontalAlignment', 'left','FontSize',9);

handles.manualPLText = manualPLText;
handles.manualPLPopup = manualPLPopup;
handles.txHeightText = txHeightText;
handles.txHeightEdit = txHeightEdit;
handles.rxHeightText = rxHeightText;
handles.rxHeightEdit = rxHeightEdit;
guidata(fig, handles);

modeChangedCallback(modePopup);

disp('MATLAB GUI is ready. Waiting for ZMQ requests...');


packet_count = 0;
last_status_update = tic;
last_drawnow_time = tic;
error_count = 0;
max_errors = 10;
avg_in_power_lin = 0;
avg_out_power_lin = 0;
avg_packet_count = 0;

while ~stopFlag
    current_time_loop = tic;

    if ~ishandle(fig)
        disp('Figure closed externally. Stopping script.');
        stopFlag = true;
        continue;
    end
    if toc(last_drawnow_time) > 0.1
        drawnow limitrate;
        last_drawnow_time = current_time_loop;
        if ~ishandle(fig)
             disp('Figure closed during drawnow. Stopping script.');
             stopFlag = true;
             continue;
        end
    end

    if updateChannelFlag
        updateChannelModel(false);
        updateChannelFlag = false;
        if ~ishandle(fig), stopFlag = true; continue; end
    end

    if pauseFlag
        if toc(last_status_update) > 0.5
             set(statusText, 'String', sprintf('Status: Paused (Pkts processed: %d)', packet_count));
             last_status_update = current_time_loop;
        end
        pause(0.01);
        continue;
    end

    if toc(last_status_update) > 1.0
        status_str = sprintf('Status: Running (Pkts: %d)', packet_count);
        if avg_packet_count > 0
            in_pwr_db = 10*log10(avg_in_power_lin / avg_packet_count);
            out_pwr_db = 10*log10(avg_out_power_lin / avg_packet_count);
            in_pwr_db_str = sprintf('%.1f', in_pwr_db); if isinf(in_pwr_db), in_pwr_db_str = '-Inf'; end
            out_pwr_db_str = sprintf('%.1f', out_pwr_db); if isinf(out_pwr_db), out_pwr_db_str = '-Inf'; end
            status_str = sprintf('%s | Avg P_in: %s dB | Avg P_out: %s dB', status_str, in_pwr_db_str, out_pwr_db_str);
            avg_in_power_lin=0; avg_out_power_lin=0; avg_packet_count=0;
        else
            status_str = sprintf('%s | Waiting for data...', status_str);
        end
        set(statusText, 'String', status_str);
        last_status_update = current_time_loop;
    end

    msg_bytes = [];
    try
        msg_bytes = socket_api_proxy.recv(ZMQ.DONTWAIT);
        if isempty(msg_bytes)
            pause(0.001);
            continue;
        end
        error_count = 0;
    catch ME_recv
        err_id = ''; if isprop(ME_recv, 'identifier'), err_id = ME_recv.identifier; end
        if contains(err_id, 'timeout') || (isprop(ME_recv, 'message') && ...
           (strcmp(ME_recv.message, 'EAGAIN') || contains(ME_recv.message,'Resource temporarily unavailable')))
             pause(0.001);
             continue;
        else
             fprintf(2,'ZMQ Recv Error: %s (ID: %s)\n', ME_recv.message, err_id);
             error_count = error_count + 1;
             if error_count >= max_errors
                 fprintf(2,'Too many consecutive ZMQ errors. Stopping script.\n');
                 stopFlag = true;
             end
             pause(0.05);
             continue;
        end
    end

    processed_bytes = uint8([]);
    processed_data_valid = false;

    if ~isempty(msg_bytes) && length(msg_bytes) >= 8
        try
            complex_signal_in = bytesToComplex(msg_bytes);

            if ~isempty(complex_signal_in)
                original_length = length(complex_signal_in);
                N_fft_qdrg = length(channelParams.H_freq_normalized);
                H_freq_norm_current = channelParams.H_freq_normalized;
                pg_sqrt_current = channelParams.pg_lin_sqrt_applied;

                avg_in_power_lin = avg_in_power_lin + mean(abs(complex_signal_in).^2);


                signal_fft = fft(complex_signal_in, original_length);
                H_applied = complex(ones(original_length, 1));

                if ~isempty(H_freq_norm_current) && N_fft_qdrg > 0
                    if original_length == N_fft_qdrg
                        H_applied = H_freq_norm_current;
                    elseif N_fft_qdrg > 1
                        f_qdrg_norm = (0:N_fft_qdrg-1)' / N_fft_qdrg;
                        f_signal_norm = (0:original_length-1)' / original_length;
                        try
                             H_applied = interp1(f_qdrg_norm, H_freq_norm_current, f_signal_norm, 'linear', H_freq_norm_current(end));
                        catch ME_interp
                             warning('H_freq interpolation failed: %s. Using flat SSF (H=1).', ME_interp.message);
                             H_applied = complex(ones(original_length, 1));
                        end
                        nan_mask_interp = isnan(H_applied);
                        if any(nan_mask_interp)
                            warning('NaN found after H_freq interpolation. Replacing with 1.');
                            H_applied(nan_mask_interp) = 1;
                        end
                    else
                         warning('Stored H_freq length <= 1. Using flat SSF (H=1).');
                         H_applied = complex(ones(original_length, 1));
                    end
                else
                    warning('Stored H_freq is empty or invalid length. Using flat SSF (H=1).');
                    H_applied = complex(ones(original_length, 1));
                end

                processed_signal_fft = signal_fft(:) .* H_applied(:);
                signal_after_multipath = ifft(processed_signal_fft, original_length);
                signal_after_pl = signal_after_multipath * pg_sqrt_current;

                noise = noise_sigma * (randn(original_length, 1) + 1i * randn(original_length, 1));
                complex_signal_out = signal_after_pl + noise;

                avg_out_power_lin = avg_out_power_lin + mean(abs(complex_signal_out).^2);
                avg_packet_count = avg_packet_count + 1;

                processed_bytes = complexToBytes(complex_signal_out);
                if isempty(processed_bytes)
                    fprintf(2,'Encoding complex signal to bytes failed.\n');
                else
                    processed_data_valid = true;
                end
            else
                 fprintf(2,'Decoding received ZMQ bytes to complex signal failed.\n');
                 processed_bytes = uint8([]);
                 processed_data_valid = true;
            end
        catch ME_process
            fprintf(2,'Error processing signal: %s @ file %s, line %d\n', ME_process.message, ME_process.stack(1).file, ME_process.stack(1).line);
            processed_bytes = uint8([]);
            processed_data_valid = true;
        end
    else
         processed_bytes = uint8([]);
         processed_data_valid = true;
    end

    try
        socket_api_proxy.send(processed_bytes);
        if processed_data_valid
            packet_count = packet_count + 1;
            error_count = 0;
        end
    catch ME_send
         fprintf(2,'ZMQ Send Error: %s\n', ME_send.message);
         error_count = error_count + 1;
         if error_count >= max_errors
             fprintf(2,'Too many consecutive ZMQ errors. Stopping script.\n');
             stopFlag = true;
         end
         pause(0.05);
    end
end

cleanup(socket_api_proxy, context, fig);




function updateChannelModel(isInitialRun)
    global channelParams currentGainText statusText fig updateChannelFlag;

    if ~ishandle(fig)
        warning('GUI figure handle is invalid during update. Aborting update.');
        updateChannelFlag = false;
        return;
    end

    if ~isInitialRun
        set(statusText, 'String', 'Status: Updating Channel Model...');
        drawnow;
    end
    disp(' ');
    disp('<<< Starting Channel Model Update >>>');
    fprintf('Mode: %s\n', channelParams.modelingMode);

    H_freq_raw = [];
    pl_dB_calculated = Inf;
    pg_lin_sqrt = 1e-10;
    pl_description = '';

    try
        H_freq_raw = getQuadrigaHfreq(channelParams);

        if isempty(H_freq_raw)
            error('Failed to get raw H_freq from QuaDRiGa.');
        end

        if strcmp(channelParams.modelingMode, 'Manual')
            pl_dB_calculated = calculateManualPathLoss(channelParams);
            pl_description = sprintf('Manual %s', channelParams.manualPL_model);

            avg_power_ssf = mean(abs(H_freq_raw).^2);
            if avg_power_ssf > 1e-30
                H_freq_normalized = H_freq_raw ./ sqrt(avg_power_ssf);
            else
                warning('Raw SSF H_freq average power near zero (%.2e). Setting normalized H flat.', avg_power_ssf);
                H_freq_normalized = complex(ones(size(H_freq_raw)));
            end

        elseif strcmp(channelParams.modelingMode, 'FullQDR')
            avg_power_full = mean(abs(H_freq_raw).^2);
            if avg_power_full > 1e-30
                pl_dB_calculated = -10 * log10(avg_power_full);
                H_freq_normalized = H_freq_raw ./ sqrt(avg_power_full);
            else
                 warning('Full QuaDRiGa H_freq average power near zero (%.2e). Setting H flat and PL to Inf.', avg_power_full);
                 pl_dB_calculated = Inf;
                 H_freq_normalized = complex(ones(size(H_freq_raw)));
            end
            pl_description = sprintf('Full QDR (%s)', channelParams.scenario);
        else
            error('Unknown modelingMode: %s', channelParams.modelingMode);
        end

        if ~isnan(pl_dB_calculated) && ~isinf(pl_dB_calculated)
            pg_lin = 10^(-pl_dB_calculated / 10);
            pg_lin_sqrt = sqrt(pg_lin);
            if isnan(pg_lin_sqrt) || isinf(pg_lin_sqrt) || pg_lin_sqrt < 0
                warning('Sqrt(Path Gain) NaN/Inf/Negative (PL=%.1f dB). Setting sqrt(gain)=1e-10.', pl_dB_calculated);
                pg_lin_sqrt = 1e-10;
            elseif pg_lin_sqrt == 0
                 warning('Sqrt(Path Gain) is zero (PL=Inf?). Setting sqrt(gain)=1e-10.');
                 pg_lin_sqrt = 1e-10;
            end
        else
            warning('Calculated PL is NaN or Inf. Applying sqrt(gain)=1e-10.');
            pg_lin_sqrt = 1e-10;
            if isnan(pl_dB_calculated), pl_dB_calculated = Inf; end
        end


        if size(H_freq_normalized, 1) == channelParams.N_fft_qdrg_calc && size(H_freq_normalized, 2) == 1
             channelParams.H_freq_normalized = H_freq_normalized;
             fprintf('Successfully updated H_freq_normalized [Size: %d x %d]\n', size(H_freq_normalized,1), size(H_freq_normalized,2));
        else
              warning('Calculated norm H size [%d x %d] != target [%d x 1]. Keeping previous H_freq.',...
                   size(H_freq_normalized,1), size(H_freq_normalized,2), channelParams.N_fft_qdrg_calc);
              if isempty(channelParams.H_freq_normalized) || length(channelParams.H_freq_normalized) ~= channelParams.N_fft_qdrg_calc
                   channelParams.H_freq_normalized = complex(ones(channelParams.N_fft_qdrg_calc, 1));
                   warning('Initialized H_freq_normalized to flat channel.');
              end
        end
        channelParams.pg_lin_sqrt_applied = pg_lin_sqrt;
        channelParams.current_pl_dB = pl_dB_calculated;

        pl_disp_str = sprintf('%.2f', pl_dB_calculated);
        if isinf(pl_dB_calculated), pl_disp_str = 'Inf'; end
        statusString = sprintf('PL: %s dB (%s) | Sqrt Gain: %.3e', ...
            pl_disp_str, pl_description, pg_lin_sqrt);
        set(currentGainText, 'String', statusString);
        if ~isInitialRun, set(statusText, 'String', 'Status: Running'); end
        drawnow;
        disp('<<< Channel Model update process finished >>>');
        disp(statusString);

    catch ME_update
         fprintf(2,'ERROR during channel model update: %s\n', ME_update.message);
         if isempty(channelParams.H_freq_normalized) || length(channelParams.H_freq_normalized) ~= channelParams.N_fft_qdrg_calc
             channelParams.H_freq_normalized = complex(ones(channelParams.N_fft_qdrg_calc, 1));
         end
         channelParams.pg_lin_sqrt_applied = 1e-10;
         channelParams.current_pl_dB = Inf;
         set(currentGainText, 'String', 'PL: Inf dB (Update Error) | Sqrt Gain: 1.000e-10');
         if ~isInitialRun, set(statusText, 'String', 'Status: ERROR during update!'); end
         drawnow;
         disp('<<< Channel Model update FAILED >>>');
    end
end


function H_freq_raw_out = getQuadrigaHfreq(params)
    H_freq_raw_out = [];

    try
        enable_internal_pl = strcmp(params.modelingMode, 'FullQDR');

        fprintf('--> Calculating QuaDRiGa H_freq (Scenario=%s, Clusters=%d, N_fft=%d, Internal PL: %s)...\n', ...
             params.scenario, params.numClusters, params.N_fft_qdrg_calc, mat2str(enable_internal_pl));

        if strcmpi(params.scenario, 'None (Flat SSF)')
            fprintf('    QuaDRiGa scenario "None". Returning flat H_freq (mag 1).\n');
            H_freq_raw_out = complex(ones(params.N_fft_qdrg_calc, 1));
            return;
        end

        s = qd_simulation_parameters;
        s.center_frequency = params.center_frequency;
        s.use_absolute_delays = true;
        s.show_progress_bars = false;
        s.samples_per_meter = 1;

        l = qd_layout(s);
        l.tx_position = params.tx_pos;
        l.rx_position = params.rx_pos;
        tx_ant = qd_arrayant('omni'); tx_ant.center_frequency = s.center_frequency; l.tx_array = tx_ant;
        rx_ant = qd_arrayant('omni'); rx_ant.center_frequency = s.center_frequency; l.rx_array = rx_ant;
        try
            l.set_scenario(params.scenario);
        catch ME_set_scen
            error('Failed to set QuaDRiGa scenario "%s": %s', params.scenario, ME_set_scen.message);
        end

        b = l.init_builder;
        if isempty(b), error('Failed to init builder.'); end
        b.scenpar.NumClusters = params.numClusters;

        if ~enable_internal_pl
            b.plpar = [];
            fprintf('    Internal QuaDRiGa PL model explicitly disabled.\n');
        else
            fprintf('    Internal QuaDRiGa PL model enabled.\n');
        end

        b.gen_parameters;
        c = b.get_channels;

        if isempty(c), error('get_channels returned empty.'); end
        if numel(c) ~= 1, warning('get_channels returned %d objects, using first.', numel(c)); end
        channel_object = c(1,1);

        if isempty(channel_object.coeff) || channel_object.no_path == 0
            warning('QuaDRiGa generated no paths (coeff empty or no_path=0).');
        end

        H_freq_raw_tmp = channel_object.fr(params.bandwidth, params.N_fft_qdrg_calc);
        if isempty(H_freq_raw_tmp), error('channel.fr returned empty.'); end

        sz = size(H_freq_raw_tmp);
        if sz(1) == params.N_fft_qdrg_calc && sz(2)==1 && sz(3)==1 % [Freq, Rx, Tx]
             H_freq_raw_out = H_freq_raw_tmp(:,1,1);
        elseif sz(1)==1 && sz(2)==1 && sz(3) == params.N_fft_qdrg_calc % [Rx, Tx, Freq]
             H_freq_raw_out = squeeze(H_freq_raw_tmp(1,1,:));
        else
             error('Unexpected H_freq dimension from channel.fr: [%s]', num2str(sz));
        end

        if size(H_freq_raw_out, 1) == 1, H_freq_raw_out = H_freq_raw_out.'; end
        if size(H_freq_raw_out, 1) ~= params.N_fft_qdrg_calc
             error('Internal error: Raw H_freq size %d != requested N_fft %d', size(H_freq_raw_out,1), params.N_fft_qdrg_calc);
        end

        nan_inf_mask_raw = isnan(H_freq_raw_out) | isinf(H_freq_raw_out);
        if any(nan_inf_mask_raw)
             warning('NaN/Inf found in raw H_freq. Replacing with eps.');
             H_freq_raw_out(nan_inf_mask_raw) = complex(eps,0);
        end

        fprintf('--> QuaDRiGa H_freq calculation finished (Raw H size: %d x %d).\n', ...
            size(H_freq_raw_out,1), size(H_freq_raw_out,2));

    catch ME_qdr
        fprintf(2, 'ERROR calculating QuaDRiGa H_freq: %s @ %s:%d\n', ME_qdr.message, ME_qdr.stack(1).name, ME_qdr.stack(1).line);
        H_freq_raw_out = [];
    end
end


function pl_dB = calculateManualPathLoss(params)
    pl_dB = Inf;
    try
        dist_3d = norm(params.rx_pos - params.tx_pos);
        d_2d    = norm(params.rx_pos(1:2) - params.tx_pos(1:2));
        fc_GHz  = params.center_frequency / 1e9;
        h_tx    = params.h_tx;
        h_rx    = params.h_rx;
        c_light = physconst('LightSpeed');

        min_dist_UMa = 35; max_dist_UMa_LOS = 5000; max_dist_UMa_NLOS = 5000;
        min_dist_UMi = 10; max_dist_UMi_LOS = 500;  max_dist_UMi_NLOS = 2000;
        min_dist_RMa = 10; max_dist_RMa_LOS = 10000; max_dist_RMa_NLOS = 10000;

        fprintf('--> Calculating Manual PL: Model=%s, Dist3D=%.1fm, Dist2D=%.1fm, f=%.2fGHz, hTx=%.1fm, hRx=%.1fm\n', ...
                params.manualPL_model, dist_3d, d_2d, fc_GHz, h_tx, h_rx);

        dist_3d_safe = max(1.0, dist_3d);
        d_2d_safe    = max(1.0, d_2d);

        h_bs_eff = max(1e-3, h_tx - 1.0);
        h_ut_eff = max(1e-3, h_rx - 1.0);
        d_BP_prime = 4 * h_bs_eff * h_ut_eff * fc_GHz * 1e9 / c_light;

        switch params.manualPL_model
            case '3GPP_UMa_NLOS'
                 d_calc = max(min_dist_UMa, d_2d);
                 if d_calc > max_dist_UMa_NLOS, fprintf('Warning: UMa NLOS dist %.1f > %.1f km limit\n', d_calc/1000, max_dist_UMa_NLOS/1000); end
                 pl1 = 32.4 + 20*log10(dist_3d_safe) + 20*log10(fc_GHz);
                 pl2 = 13.54 + 39.08*log10(dist_3d_safe) + 20*log10(fc_GHz) - 0.6*(h_rx-1.5);
                 pl_dB = max(pl1, pl2);
            case '3GPP_UMa_LOS'
                 d_calc = max(min_dist_UMa, d_2d);
                 d_BP = max(min_dist_UMa, d_BP_prime);
                 if d_calc <= d_BP
                     pl_dB = 28.0 + 22 * log10(dist_3d_safe) + 20 * log10(fc_GHz);
                 else
                     if d_calc > max_dist_UMa_LOS, fprintf('Warning: UMa LOS dist %.1f > %.1f km limit\n', d_calc/1000, max_dist_UMa_LOS/1000); end
                     d_BP_safe = max(1e-9, d_BP);
                     pl_dB = 28.0 + 40 * log10(dist_3d_safe) + 20 * log10(fc_GHz) ...
                           - 9 * log10( d_BP_safe^2 + (h_tx-h_rx)^2 );
                 end
            case '3GPP_UMi_NLOS'
                 d_calc = max(min_dist_UMi, d_2d);
                 if d_calc > max_dist_UMi_NLOS, fprintf('Warning: UMi NLOS dist %.1f > %.1f km limit\n', d_calc/1000, max_dist_UMi_NLOS/1000); end
                 pl1 = 32.4 + 20*log10(dist_3d_safe) + 20*log10(fc_GHz);
                 pl2 = 35.3 * log10(dist_3d_safe) + 22.4 + 21.3 * log10(fc_GHz) - 0.3*(h_rx-1.5);
                 pl_dB = max(pl1, pl2);
            case '3GPP_UMi_LOS'
                 d_calc = max(min_dist_UMi, d_2d);
                 d_BP = max(min_dist_UMi, d_BP_prime);
                 if d_calc <= d_BP
                     pl_dB = 32.4 + 21 * log10(dist_3d_safe) + 20 * log10(fc_GHz);
                 else
                      if d_calc > max_dist_UMi_LOS, fprintf('Warning: UMi LOS dist %.1f > %.1f m limit\n', d_calc, max_dist_UMi_LOS); end
                      d_BP_safe = max(1e-9, d_BP);
                      pl_dB = 32.4 + 40 * log10(dist_3d_safe) + 20 * log10(fc_GHz) ...
                            - 9.5 * log10(d_BP_safe^2 + (h_tx-h_rx)^2);
                 end
             case '3GPP_RMa_NLOS'
                 d_calc = max(min_dist_RMa, d_2d);
                 if d_calc > max_dist_RMa_NLOS, fprintf('Warning: RMa NLOS dist %.1f > %.1f km limit\n', d_calc/1000, max_dist_RMa_NLOS/1000); end
                 pl1 = 32.4 + 20*log10(dist_3d_safe) + 20*log10(fc_GHz);
                 pl2 = 161.04 - 7.1*log10(20) + 7.5*log10(5) ...
                       - (24.37 - 3.7*(5/h_tx)^2)*log10(h_tx) ...
                       + (43.42 - 3.1*log10(h_tx))*(log10(dist_3d_safe) - 3) ...
                       + 20*log10(fc_GHz) - (3.2*(log10(11.75*h_rx))^2 - 4.97);
                 pl_dB = max(pl1, pl2);
             case '3GPP_RMa_LOS'
                 d_calc = max(min_dist_RMa, d_2d);
                 d_BP = max(min_dist_RMa, d_BP_prime);
                 h_eff = 5;
                 if d_calc <= d_BP
                      pl_dB = 20 * log10(40 * pi * dist_3d_safe * fc_GHz / 3) ...
                            + min(0.03 * h_eff^1.72, 10) * log10(dist_3d_safe) ...
                            - min(0.044 * h_eff^1.72, 14.77) ...
                            + 0.002 * log10(h_eff) * dist_3d_safe;
                 else
                      if d_calc > max_dist_RMa_LOS, fprintf('Warning: RMa LOS dist %.1f > %.1f km limit\n', d_calc/1000, max_dist_RMa_LOS/1000); end
                      d_BP_safe = max(1e-9, d_BP);
                      pl_at_bp = 20 * log10(40 * pi * d_BP_safe * fc_GHz / 3) ...
                               + min(0.03 * h_eff^1.72, 10) * log10(d_BP_safe) ...
                               - min(0.044 * h_eff^1.72, 14.77) ...
                               + 0.002 * log10(h_eff) * d_BP_safe;
                      pl_dB = pl_at_bp + 40 * log10(dist_3d_safe / d_BP_safe);
                 end
            case 'FSPL'
                 lambda = c_light / params.center_frequency;
                 pl_lin = ( (4 * pi * dist_3d_safe) / lambda )^2;
                 if pl_lin < 1.0, pl_lin = 1.0; end
                 pl_dB = 10*log10(pl_lin);
            otherwise
                warning('Unknown manual PL model: "%s". Using FSPL.', params.manualPL_model);
                 lambda = c_light / params.center_frequency;
                 pl_lin = ( (4 * pi * dist_3d_safe) / lambda )^2;
                 if pl_lin < 1.0, pl_lin = 1.0; end
                 pl_dB = 10*log10(pl_lin);
        end

        if isnan(pl_dB) || isinf(pl_dB) || pl_dB < 0
            if pl_dB < 0, warning('Manual PL = %.2f dB < 0. Clipping to 0 dB.', pl_dB); pl_dB = 0;
            else warning('Manual PL NaN/Inf for model %s. Returning Inf.', params.manualPL_model); pl_dB = Inf; end
        else
            fprintf('--> Manual PL calculated: %.2f dB\n', pl_dB);
        end

    catch ME_pl
        fprintf(2,'ERROR calculating manual PL for "%s": %s @ %s:%d\n', ...
            params.manualPL_model, ME_pl.message, ME_pl.stack(1).name, ME_pl.stack(1).line);
        pl_dB = Inf;
    end
end


% --- GUI ---
function togglePause(~, ~)
    global pauseFlag;
    pauseFlag = ~pauseFlag;
    if pauseFlag, disp('Script Paused.'); else disp('Script Resumed.'); end
end

function stopScript(~, ~)
    global stopFlag;
    disp('STOP button pressed. Script will terminate.');
    stopFlag = true;
end

function closeFigureCallback(src, ~)
    global stopFlag;
    disp('Figure close requested. Stopping script.');
    stopFlag = true;
    try if ishandle(src), delete(src); end; catch; end
end

function modeChangedCallback(popupHandle, ~)
    handles = guidata(popupHandle);
    selectedIndex = get(popupHandle, 'Value');
    modeOptions = get(popupHandle, 'String');
    selectedMode = modeOptions{selectedIndex};
    isManualMode = strcmp(selectedMode, 'Manual PL + QDR SSF');
    new_state = 'off';
    if isManualMode, new_state = 'on'; end

    set(handles.manualPLText, 'Enable', new_state);
    set(handles.manualPLPopup, 'Enable', new_state);
    set(handles.txHeightText, 'Enable', new_state);
    set(handles.txHeightEdit, 'Enable', new_state);
    set(handles.rxHeightText, 'Enable', new_state);
    set(handles.rxHeightEdit, 'Enable', new_state);
    disp(['Modeling mode changed. Manual PL controls: ', new_state]);
end

function applyParametersCallback(src, ~, modePopup, distEdit, scenarioPopup, clustersEdit, manualPLPopup, txHeightEdit, rxHeightEdit)
    global channelParams updateChannelFlag stopFlag;
    if stopFlag, disp('Update request ignored: Stop flag is set.'); return; end
    disp(' ');
    disp('>>> Apply & Update Channel button pressed <<<');

    try
        modeOptions = get(modePopup, 'String');
        selectedModeIndex = get(modePopup, 'Value');
        new_modelingMode = modeOptions{selectedModeIndex};
        if strcmp(new_modelingMode, 'Manual PL + QDR SSF'), new_modelingMode = 'Manual'; else new_modelingMode = 'FullQDR'; end

        new_dist_str = get(distEdit, 'String');
        new_dist_2d = str2double(new_dist_str);

        scenarioItems = get(scenarioPopup, 'String');
        selectedIndex = get(scenarioPopup, 'Value');
        new_scenario = scenarioItems{selectedIndex};

        new_clusters_str = get(clustersEdit, 'String');
        new_clusters = round(str2double(new_clusters_str));

        manualPLItems = get(manualPLPopup, 'String');
        selectedPLIndex = get(manualPLPopup, 'Value');
        new_manualPL_model = manualPLItems{selectedPLIndex};

        new_h_tx_str = get(txHeightEdit, 'String');
        new_h_tx = str2double(new_h_tx_str);

        new_h_rx_str = get(rxHeightEdit, 'String');
        new_h_rx = str2double(new_h_rx_str);

        validInput = true;
        if isnan(new_dist_2d) || new_dist_2d < 0
            errordlg('Invalid 2D distance.', 'Input Error', 'modal'); validInput = false;
        end
         if isnan(new_clusters) || new_clusters <= 0 || floor(new_clusters) ~= new_clusters
            errordlg('Invalid #Clusters.', 'Input Error', 'modal'); validInput = false;
        end
         if strcmp(new_modelingMode,'Manual')
             if isnan(new_h_tx) || new_h_tx <= 0
                errordlg('Invalid Tx Height.', 'Input Error', 'modal'); validInput = false;
             end
             if isnan(new_h_rx) || new_h_rx <= 0
                errordlg('Invalid Rx Height.', 'Input Error', 'modal'); validInput = false;
             end
         end
         if ~validInput
             disp('Input validation failed.');
             return;
         end

        tx_pos_current = channelParams.tx_pos;
        rx_pos_current = channelParams.rx_pos;
        direction_vector_2d = rx_pos_current(1:2) - tx_pos_current(1:2);
        current_dist_2d = norm(direction_vector_2d);
        if abs(current_dist_2d) < 1e-9
             warning('Current 2D distance near zero. Assuming direction along X-axis.');
             direction_vector_2d = [1; 0];
        else
             direction_vector_2d = direction_vector_2d / current_dist_2d;
        end
        new_rx_pos_xy = tx_pos_current(1:2) + direction_vector_2d * new_dist_2d;
        new_rx_pos = [new_rx_pos_xy; new_h_rx];
        new_tx_pos = [tx_pos_current(1:2); new_h_tx];

        mode_changed = ~strcmp(channelParams.modelingMode, new_modelingMode);
        distance_changed = abs(current_dist_2d - new_dist_2d) > 1e-6 * max([1, current_dist_2d, new_dist_2d]);
        scenario_changed = ~strcmp(channelParams.scenario, new_scenario);
        clusters_changed = abs(channelParams.numClusters - new_clusters) > 0;
        pl_model_changed = ~strcmp(channelParams.manualPL_model, new_manualPL_model);
        tx_h_changed = abs(channelParams.h_tx - new_h_tx) > 1e-6;
        rx_h_changed = abs(channelParams.h_rx - new_h_rx) > 1e-6;

        needs_update = mode_changed;
        if strcmp(new_modelingMode,'Manual')
            needs_update = needs_update || distance_changed || scenario_changed || clusters_changed || pl_model_changed || tx_h_changed || rx_h_changed;
        elseif strcmp(new_modelingMode,'FullQDR')
            needs_update = needs_update || distance_changed || scenario_changed || tx_h_changed || rx_h_changed || clusters_changed;
        end

        if needs_update
            disp('Parameters changed. Applying new values:');
            channelParams.modelingMode = new_modelingMode;
            channelParams.rx_pos = new_rx_pos;
            channelParams.tx_pos = new_tx_pos;
            channelParams.scenario = new_scenario;
            channelParams.numClusters = new_clusters;
            channelParams.manualPL_model = new_manualPL_model;
            channelParams.h_tx = new_h_tx;
            channelParams.h_rx = new_h_rx;

            fprintf('  Mode: %s\n', channelParams.modelingMode);
            fprintf('  Distance (2D): %.1f m\n', new_dist_2d);
            fprintf('  QDR Scenario: %s\n', channelParams.scenario);
            fprintf('  QDR #Clusters: %d\n', channelParams.numClusters);
            if strcmp(channelParams.modelingMode, 'Manual')
                fprintf('  Manual PL Model: %s\n', channelParams.manualPL_model);
                fprintf('  Tx Height (Manual): %.1f m\n', channelParams.h_tx);
                fprintf('  Rx Height (Manual): %.1f m\n', channelParams.h_rx);
            end
            fprintf('  New Tx pos (used by QDR): [%.1f, %.1f, %.1f]\n', channelParams.tx_pos(1), channelParams.tx_pos(2), channelParams.tx_pos(3));
            fprintf('  New Rx pos (used by QDR): [%.1f, %.1f, %.1f]\n', channelParams.rx_pos(1), channelParams.rx_pos(2), channelParams.rx_pos(3));

            updateChannelFlag = true;
            disp('Channel model update requested.');

            set(distEdit, 'String', sprintf('%.1f', new_dist_2d));
            set(clustersEdit,'String', num2str(new_clusters));
            set(txHeightEdit,'String', sprintf('%.1f', new_h_tx));
            set(rxHeightEdit,'String', sprintf('%.1f', new_h_rx));
            modeChangedCallback(modePopup);
        else
            disp('Input parameters effectively unchanged. No update requested.');
            set(distEdit, 'String', sprintf('%.1f', new_dist_2d));
            set(clustersEdit,'String', num2str(new_clusters));
            set(txHeightEdit,'String', sprintf('%.1f', new_h_tx));
            set(rxHeightEdit,'String', sprintf('%.1f', new_h_rx));
        end
        disp('>>> Apply & Update Channel finished <<<');

    catch ME_apply
        fprintf(2,'Error applying parameters: %s\n', ME_apply.message);
        errordlg(sprintf('Error applying parameters: %s', ME_apply.message), 'Parameter Error', 'modal');
        try
             dist_2d_old = norm(channelParams.tx_pos(1:2) - channelParams.rx_pos(1:2));
             set(distEdit, 'String', sprintf('%.1f', dist_2d_old));
             set(clustersEdit,'String', num2str(channelParams.numClusters));
             set(txHeightEdit,'String', sprintf('%.1f', channelParams.h_tx));
             set(rxHeightEdit,'String', sprintf('%.1f', channelParams.h_rx));
             if strcmp(channelParams.modelingMode,'Manual'), set(modePopup,'Value',1); else set(modePopup,'Value',2); end
             idx = find(strcmp(get(scenarioPopup, 'String'), channelParams.scenario), 1); if ~isempty(idx), set(scenarioPopup, 'Value', idx); end
             idx = find(strcmp(get(manualPLPopup, 'String'), channelParams.manualPL_model), 1); if ~isempty(idx), set(manualPLPopup, 'Value', idx); end
             modeChangedCallback(modePopup);
             warning('GUI values restored after error.');
         catch ME_restore
              fprintf(2, 'Could not restore GUI values after error: %s\n', ME_restore.message);
         end
    end
end

% --- Преобразование Данных ---
function complex_array = bytesToComplex(byte_array)
    complex_array = [];
    if isempty(byte_array) || mod(length(byte_array), 8) ~= 0
        % fprintf(2,'bytesToComplex: Input byte array empty or incorrect length (%d bytes).\n', length(byte_array));
        return;
    end
    try
        floatArray = typecast(uint8(byte_array(:)), 'single');
        complex_array = complex(floatArray(1:2:end), floatArray(2:2:end));
        mask = isnan(complex_array) | isinf(complex_array);
        if any(mask)
            % warning('bytesToComplex: NaN/Inf detected. Replacing with 0.');
            complex_array(mask) = 0;
        end
    catch ME_bytes
        fprintf(2,'Error in bytesToComplex: %s\n', ME_bytes.message);
        complex_array = [];
    end
end

function byte_array = complexToBytes(complex_array)
    byte_array = uint8([]);
    if isempty(complex_array), return; end
    try
        single_array = single(complex_array(:));
        real_part = real(single_array);
        imag_part = imag(single_array);
        mask_r = isnan(real_part) | isinf(real_part);
        mask_i = isnan(imag_part) | isinf(imag_part);
        if any(mask_r), % warning('complexToBytes: NaN/Inf in real part. Replacing with 0.');
            real_part(mask_r)=0; end
        if any(mask_i), % warning('complexToBytes: NaN/Inf in imag part. Replacing with 0.');
            imag_part(mask_i)=0; end
        floatArray = [real_part.'; imag_part.'];
        byte_array = typecast(floatArray(:), 'uint8');
    catch ME_complex
        fprintf(2,'Error in complexToBytes: %s\n', ME_complex.message);
        byte_array = uint8([]);
    end
end

% --- Очистка при выключении ---
function cleanup(socket, context, figHandle)
    disp(' ');
    disp('<<< Cleaning up ZMQ and GUI resources... >>>');
    global stopFlag;
    stopFlag = true;

    if ~isempty(socket)
        disp('Closing ZMQ socket...');
        try socket.close(); disp('ZMQ socket closed.');
        catch ME_sock_close, fprintf(2,'Warn: Error closing ZMQ socket: %s\n', ME_sock_close.message); end
    end
    if ~isempty(context) && ismethod(context,'term')
        disp('Terminating ZMQ context...');
        try context.term(); disp('ZMQ context terminated.');
        catch ME_ctx_term, fprintf(2,'Warn: Error terminating ZMQ context: %s\n', ME_ctx_term.message); end
    end
    if ~isempty(figHandle) && ishandle(figHandle)
        disp('Closing GUI figure...');
        try delete(figHandle); disp('GUI figure closed.');
        catch ME_fig_close, fprintf(2,'Warn: Error closing figure: %s\n', ME_fig_close.message); end
    end
    clear global pauseFlag stopFlag updateChannelFlag channelParams;
    disp('Global variables cleared.');
    disp('<<< Cleanup complete. MATLAB script finished. >>>');
end