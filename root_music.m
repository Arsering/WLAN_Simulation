function [angle_info_output] = root_music(CSI,Model_paras,signal_space)
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
    %�Զ��sample
    row = floor((size(CSI,1)+1)/2);
    column = size(CSI,1) + 1 - row;
    correlation_matrix = complex(zeros(row,row));
    
    for c = 1:column
        correlation_matrix = correlation_matrix + CSI(c:(c + row - 1),:) * CSI(c:(c + row - 1),:)';
    end
end

%% 

% ��������ֵ�����Ӧ����������
[E,D] = eig(correlation_matrix);

 % �ҵ�noise_space��Ӧ����������
[~,indx] = sort(diag(D),'descend');
eigenvects = E(:,indx);
noise_eigenvects = eigenvects(:,(signal_space+1):end);

C = noise_eigenvects * noise_eigenvects';


CL = zeros(1,size(C,1)-1);
for L = 1:size(C,1)-1
    for N = 1:size(C,1)-L
        CL(L) = CL(L) + C(N,N+L);
    end
end
pseudo_spectrum_inverse = [fliplr(CL),sum(diag(C)),conj(CL)];
Z = roots(pseudo_spectrum_inverse);
[absult_Z,I] = sort(abs(Z));
Z = Z(I);
[~,I] = sort(abs(absult_Z(1:floor(end /2)) - 1));
Z = Z(I);

% ��������֮����� ��λ��m
antenna_space = (Model_paras.speed_light/Model_paras.frequency) * Model_paras.antenna_space_ofwaveLen;

angle_info_output = zeros(1,signal_space);
for k = 1:signal_space
    angle_info_output(k) = rad2deg(acos( - angle(Z(k)) * Model_paras.speed_light / (2*pi * antenna_space * Model_paras.frequency)));
end

end
%% ����������������Ԫ�ص��±�

function index = maxI(input)
    index  = 1;
    for k = 2:length(input)
        if input(k) > input(index)
            index = k;
        end
    end
end

