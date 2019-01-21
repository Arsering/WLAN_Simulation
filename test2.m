%% ���������Եõ��ڲ�ͬ�������ֵ֮�� �����ض���·������ AOA��TOF�����ֲ�
%% �趨����
num_path = 1; % ·������
num_AP = 3; % ���ж��ٸ�AP
AOA_or_TOF = 2;% ����1��ʾ��AOA�Ĵ����ʣ�����2��ʾ��TOF�Ĵ�����
num_experiment = 500; % ÿ��SNR��ʵ�����
has_noise = 1; % Ϊ1��������� �����������
SNR_vector = [50,30,20,10,0]; % ���Ե�SNR��Χ
error = zeros(length(SNR_vector),num_experiment); % �洢�������
mark = {'r-+' 'g--o' 'b*:' 'c-.s' 'm-p' 'b--h'}; % ��ͬ����ͼ���в�ͬ����ɫ�����͡�����

%% �õ����ֵ
for m = 1:length(SNR_vector)
    for t = 1:num_experiment
        AP_index = mod(t,num_AP) + 1; % AP������
        [path_info_input,path_info_output] = simulation_environment(num_path,SNR_vector(m),AP_index,has_noise);

        %�����·����Ӧ���±�
        minTOF_index = 1;
        for k = 2:size(path_info_output,1)
            if path_info_output(k,2) < path_info_output(minTOF_index,2)
                minTOF_index = k;
            end
        end

        error(m,t) = abs(path_info_output(minTOF_index,AOA_or_TOF) - path_info_input(1,AOA_or_TOF)); % 
        
    end
end


%% ����ͼ��
    
%���AOA������������
format longE
max_error = max(max(error));
min_error = 0;

x_vector = linspace(min_error,max_error,15);
y_vector = zeros(length(SNR_vector),length(x_vector));

%��error�ķֲ�
for m = 1:length(SNR_vector)
    for t = 2:length(x_vector)
        for k = 1:length(error(m,:))
            if error(m,k) <= x_vector(t)
                y_vector(m,t) =  y_vector(m,t) + 1;
            end
        end
    end
end


y_vector = y_vector ./ num_experiment;
figure;
hold on;
for m = 1:length(SNR_vector)
    sign = [' SNR ',num2str(SNR_vector(m)),'dB'];
    plot(x_vector,y_vector(m,:),mark{m},'DisplayName',sign,'LineWidth',0.75);
end
legend();
if AOA_or_TOF == 1
    xlab = 'AOA';
else
    xlab = 'TOF';
end
xlabel([xlab,' error']);  %x����������
ylabel('CDF'); %y����������
hold off;

