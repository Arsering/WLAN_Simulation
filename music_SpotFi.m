function [path_info_output] = music_SpotFi(CSI,WLAN_paras,signal_space)
%% smoothing
row = floor(WLAN_paras.num_subcarrier/2) * 2;% smoothed_CSI的行数
column=(WLAN_paras.num_subcarrier - row/2 + 1) * (WLAN_paras.num_antenna - 1);% smoothed_CSI的列数

if size(CSI,3) == 1
    %定义smoothed_CSI的存储矩阵
smoothed_CSI = zeros(row,column,'like',CSI);

    %构造smoothed_CSI
    for k = 1:(WLAN_paras.num_antenna-1)
        for t = 1:(WLAN_paras.num_subcarrier - row/2 + 1)
            smoothed_CSI(:,t + (k-1)*(WLAN_paras.num_subcarrier - row/2 + 1)) = [CSI(k,t:(t+row/2-1)),CSI(k+1,t:(t+row/2-1))].';
        end
    end
    correlation_matrix = smoothed_CSI * smoothed_CSI';
else
    %  定义存储相关矩阵的空间
    correlation_matrix = complex(zeros(row,row));
    
    %构造smoothed_CSI
    for k = 1:(WLAN_paras.num_antenna-1)
        for t = 1:(WLAN_paras.num_subcarrier - row/2 + 1)
            tmp_CSI = [squeeze(CSI(k,t:(t+row/2-1),:));squeeze(CSI(k+1,t:(t+row/2-1),:))];
            correlation_matrix = correlation_matrix + tmp_CSI * tmp_CSI';
        end
    end
end


%% 算法主体
% 计算smoothed_CSI的特征值及其对应的特征向量
[E,D] = eig(correlation_matrix);

 % 找到noise_space对应的特征向量
[~,indx] = sort(diag(D),'descend');
eigenvects = E(:,indx);
noise_eigenvects = eigenvects(:,(signal_space+1):end);

antenna_space = (WLAN_paras.speed_light/WLAN_paras.frequency) * WLAN_paras.antenna_space_ofwaveLen; % 相邻天线之间距离 单位：m

%确定采样点
X = 0:WLAN_paras.precison(1):180;
Y = WLAN_paras.precison(2):WLAN_paras.precison(2):40;
pseudo_spectrum = complex(zeros(length(X),length(Y)));

%采样
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
%% 生成三维图像
figure;
mesh(Y,X,pseudo_spectrum);
xlabel('X(TOF/m)');
ylabel('Y(AOA/°)');
zlabel('pseudo-spectrum(dB)');
shading interp;

%% 
%定义存储求得的AOA TOF的矩阵
path_info_output = zeros(signal_space,2);
max_N_value = zeros(1,signal_space);

format longE
%寻找前signal_space个极大值点
for m = 2:length(X)-1
    for n = 2:length(Y)-1
        scope = [length(X),length(Y)];
        mark = 1;

        %判断当前点是否为极大值点
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
       
        %如果为极大值点，则存储起来
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

%% 求得输入数组中最小元素的下标

function index = minI(input)
    index  = 1;
    for k = 2:length(input)
        if input(k) < input(index)
            index = k;
        end
    end
end

