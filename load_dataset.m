function data_dict = load_dataset(file_path)
    % 加载 .mat 文件
    dataset = load(file_path);
    
    % 将数据存储到结构体中
    data_dict = struct();
    data_dict.theta_G_t = double(dataset.theta_G_t);
    data_dict.phi_G_t = double(dataset.phi_G_t);
    data_dict.theta_L_t = double(dataset.theta_L_t);
    data_dict.phi_L_t = double(dataset.phi_L_t);
    data_dict.d_single_G = double(dataset.d_single_G);
    data_dict.d_single_L = double(dataset.d_single_L);
    data_dict.theta_GA_t = double(dataset.theta_GA_t);
    data_dict.phi_GA_t = double(dataset.phi_GA_t);
    data_dict.theta_LA_t = double(dataset.theta_LA_t);
    data_dict.phi_LA_t = double(dataset.phi_LA_t);
    data_dict.d_GA = double(dataset.d_GA);
    data_dict.d_LA = double(dataset.d_LA);
    data_dict.theta_JA_t = double(dataset.theta_JA_t);
    data_dict.phi_JA_t = double(dataset.phi_JA_t);
    data_dict.d_JA = double(dataset.d_JA);
    data_dict.theta_GA_r = double(dataset.theta_GA_r);
    data_dict.phi_GA_r = double(dataset.phi_GA_r);
    data_dict.theta_LA_r = double(dataset.theta_LA_r);
    data_dict.phi_LA_r = double(dataset.phi_LA_r);
    data_dict.theta_JA_r = double(dataset.theta_JA_r);
    data_dict.phi_JA_r = double(dataset.phi_JA_r);
end