path_info_input = [60,15;...
    90,32;...
    130,21];
SNR = 10;
has_noise = 1;
num_samples = 1;
correlation_coefficient = 0.99;
format longE
path_info_output = simulation_environment(path_info_input,SNR,has_noise,num_samples,correlation_coefficient);
