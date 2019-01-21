function [path_info_output,samples] = MUSIC_Origin(CSI,Model_paras,signal_space)
%% ��ԭʼ��MUSIC��û�н���ƽ��
% 
% [row,column] = size(CSI);
% correlation_matrix = CSI * CSI';
%% spatial smoothing
if size(CSI,2) == 1
    % ���ڵ���sample
    row = floor((length(CSI)+1) / 2);
    column = length(CSI)+1 - row;
    
    smoothed_CSI = complex(zeros(row,column));
    
    for c = 1:column
        smoothed_CSI(:,c) = CSI(c:(row+c-1));
    end
    CSI = smoothed_CSI;
    
    correlation_matrix = CSI * CSI';
else
    % ���ڶ��sample
    row = floor((size(CSI,1)+1) / 2);
    column = size(CSI,1) + 1 - row;
    correlation_matrix = complex(zeros(row,row));
    
    for c = 1:column
        correlation_matrix = correlation_matrix + CSI(c:(c + row - 1),:) * CSI(c:(c + row - 1),:)';
    end
end
%% �㷨����

% ����CSI������ֵ�����Ӧ����������

[E,D] = eig(correlation_matrix);

 % �ҵ�noise_space��Ӧ����������
[~,indx] = sort(diag(D),'descend');
eigenvects = E(:,indx);
noise_eigenvects = eigenvects(:,(signal_space+1):end);


% ������������֮��ľ���
antenna_space = (Model_paras.speed_light/Model_paras.frequency) * Model_paras.antenna_space_ofwaveLen; % ��������֮����� ��λ��m

%ȷ��������
X = (0:Model_paras.precison(1):180);

samples = complex(zeros(1,length(X)));

%����
for t = 1:length(X)
    Steering_Vector = complex(zeros(row,1));
    for m = 0:1:row-1
        Steering_Vector(m+1) = exp(-1i * 2*pi * antenna_space * cos(X(t) * pi / 180) * Model_paras.frequency / Model_paras.speed_light)^m;
    end
    samples(t) = 1/sum(abs(noise_eigenvects' * Steering_Vector).^2,1);

end

pseudo_spectrum = 20 * log10(samples);
%% ���ɶ�άͼ��

figure
plot(X,pseudo_spectrum,'m--');
title('MUSIC(no_smoothing)');
xlabel('angle��degree��');  %x����������
ylabel('pseudo-spectrum(dB)'); %y����������


%% 
%����洢��õ�AOA TOF�ľ���
path_info_output = zeros(1,signal_space);
max_N_value = zeros(1,signal_space);

%Ѱ��ǰsignal_space������ֵ��
for m = 1:length(samples)
    scopeX = length(X);
    for range = - Model_paras.window_size:1:Model_paras.window_size
        tmpX = m + range;
        if tmpX > scopeX || tmpX < 1
            tmpX = m;
        end
        mark = 1;
        if samples(m) < samples(tmpX)
            mark = 0;
            break;
        end
    end
    if mark == 1
        tmp_index = minI(max_N_value);
        if max_N_value(tmp_index) < samples(m)
            path_info_output(tmp_index) = X(m);
            max_N_value(tmp_index) = samples(m);
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

