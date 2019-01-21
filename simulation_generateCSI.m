function CSI = simulation_generateCSI(WLAN_paras)

% ��������֮����� ��λ��m
antenna_space = (WLAN_paras.speed_light/WLAN_paras.frequency) * WLAN_paras.antenna_space_ofwaveLen;

% ����洢CSIֵ�ľ���
CSI = complex(zeros(WLAN_paras.num_antenna,WLAN_paras.num_subcarrier,WLAN_paras.num_samples));

complex_gain = create_complexGain(WLAN_paras);

%����CSI
for s = 1:WLAN_paras.num_samples
    for k = 1:size(WLAN_paras.path_info_input,1)
        exp_AoA = exp((-1i * 2 * pi * antenna_space * cos(deg2rad(WLAN_paras.path_info_input(k,1))) * WLAN_paras.frequency / WLAN_paras.speed_light));
        exp_ToF = exp( -1i * 2 * pi * WLAN_paras.frequency_space * WLAN_paras.path_info_input(k,2) / WLAN_paras.speed_light);
        
        for t = 1:WLAN_paras.num_antenna
            tmp_AoA = exp_AoA.^(t-1);
            CSI(t,:,s) = CSI(t,:,s) + exp_ToF.^(0:(WLAN_paras.num_subcarrier-1)) * tmp_AoA * complex_gain(k,s);
        end
    end
end

if WLAN_paras.has_noise == 1
    CSI = awgn(CSI,WLAN_paras.SNR,'measured');
end

end

%% �����ض���complex_gain
function complex_gain = create_complexGain(Model_paras)

% ���ɳ�ʼ��Э�������
signalCovMat = 1*eye(size(Model_paras.path_info_input,1));

% Ϊ�˼�� �ҽ�����������ͬ��Դ֮������ϵ������ͬһ��ֵ
for t = 1:size(Model_paras.path_info_input,1)
    for k = 1:size(Model_paras.path_info_input,1)
        if t == k
            continue;
        else
            signalCovMat(t,k) = Model_paras.correlation_coefficient;
        end
    end
end

complex_gain = mvnrnd(zeros(size(Model_paras.path_info_input,1),1),signalCovMat,Model_paras.num_samples).';

% Ϊ�˼������ �����޸ģ������޸���񣬶���Ӱ�����ǵ�ʵ��Ч����
complex_gain = 2 * exp(1i * complex_gain);

end

