function [path_info_output] = music_SpotFi(CSI,WLAN_paras,signal_space)
%% smoothing
row = floor(WLAN_paras.num_subcarrier/2) * 2;% smoothed_CSI������
column=(WLAN_paras.num_subcarrier - row/2 + 1) * (WLAN_paras.num_antenna - 1);% smoothed_CSI������

if size(CSI,3) == 1
    %����smoothed_CSI�Ĵ洢����
smoothed_CSI = zeros(row,column,'like',CSI);

    %����smoothed_CSI
    for k = 1:(WLAN_paras.num_antenna-1)
        for t = 1:(WLAN_paras.num_subcarrier - row/2 + 1)
            smoothed_CSI(:,t + (k-1)*(WLAN_paras.num_subcarrier - row/2 + 1)) = [CSI(k,t:(t+row/2-1)),CSI(k+1,t:(t+row/2-1))].';
        end
    end
    correlation_matrix = smoothed_CSI * smoothed_CSI';
else
    %  ����洢��ؾ���Ŀռ�
    correlation_matrix = complex(zeros(row,row));
    
    %����smoothed_CSI
    for k = 1:(WLAN_paras.num_antenna-1)
        for t = 1:(WLAN_paras.num_subcarrier - row/2 + 1)
            tmp_CSI = [squeeze(CSI(k,t:(t+row/2-1),:));squeeze(CSI(k+1,t:(t+row/2-1),:))];
            correlation_matrix = correlation_matrix + tmp_CSI * tmp_CSI';
        end
    end
end


%% �㷨����
% ����smoothed_CSI������ֵ�����Ӧ����������
[E,D] = eig(correlation_matrix);

 % �ҵ�noise_space��Ӧ����������
[~,indx] = sort(diag(D),'descend');
eigenvects = E(:,indx);
noise_eigenvects = eigenvects(:,(signal_space+1):end);

antenna_space = (WLAN_paras.speed_light/WLAN_paras.frequency) * WLAN_paras.antenna_space_ofwaveLen; % ��������֮����� ��λ��m

%ȷ��������
X = 0:WLAN_paras.precison(1):180;
Y = WLAN_paras.precison(2):WLAN_paras.precison(2):40;
pseudo_spectrum = complex(zeros(length(X),length(Y)));

%����
for t = 1:length(X)
    for k = 1:length(Y)
        angleE = exp(-1i * 2 * pi * antenna_space * cos(deg2rad(X(t))) * WLAN_paras.frequency / WLAN_paras.speed_light);
        timeE = exp(-1i * 2 * pi * WLAN_paras.frequency_space * Y(k) / WLAN_paras.speed_light);
        steering_vector = complex(zeros(row,1));
        for n=0:1
            for m=1:row/2
                steering_vector((n*row/2)+m,1) = angleE.^n * timeE.^(m-1);
            end
        end
        pseudo_spectrum(t,k) = 1 / sum(abs(noise_eigenvects' * steering_vector).^2,1);
    end
end

pseudo_spectrum = 20*log10(pseudo_spectrum);
%% ������άͼ��
figure;
mesh(Y,X,pseudo_spectrum);
xlabel('X(TOF/m)');
ylabel('Y(AOA/��)');
zlabel('pseudo-spectrum(dB)');
shading interp;

%% 
%����洢��õ�AOA TOF�ľ���
path_info_output = zeros(signal_space,2);
max_N_value = zeros(1,signal_space);

format longE
%Ѱ��ǰsignal_space������ֵ��
for m = 2:length(X)-1
    for n = 2:length(Y)-1
        scope = [length(X),length(Y)];
        mark = 1;

        %�жϵ�ǰ���Ƿ�Ϊ����ֵ��
        for range_x = - WLAN_paras.window_size:1:WLAN_paras.window_size
            for range_y = - WLAN_paras.window_size:1:WLAN_paras.window_size
                temp_x = m + range_x;
                if temp_x < 1 || temp_x > scope(1)
                    temp_x = m;
                end
                temp_y = n + range_y;
                if temp_y < 1 ||temp_y > scope(2)
                    temp_y = n;
                end
                if pseudo_spectrum(m,n) < pseudo_spectrum(temp_x,temp_y)
                    mark = 0;
                    break;
                end
            end
            if mark == 0
                break;
            end
        end
       
        %���Ϊ����ֵ�㣬��洢����
        if mark == 1
            min_index = minI(max_N_value);
            if max_N_value(min_index) < pseudo_spectrum(m,n)
                max_N_value(min_index) =  pseudo_spectrum(m,n);
                path_info_output(min_index,:) = [X(m) Y(n)];
            end
        end
    end
end

end

%% ���������������СԪ�ص��±�

function index = minI(input)
    index  = 1;
    for k = 2:length(input)
        if input(k) < input(index)
            index = k;
        end
    end
end

