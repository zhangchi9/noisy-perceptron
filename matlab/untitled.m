
noise_added = 30;

% file_obj = dir('/scratch/zhang.chi9/perceptron/data/different_network_size_load_at_capacity/perceptron*200*');
% get_retrieval_prob_length(file_obj,noise_added,'perceptron_200_retrieval.mat')
% get_dynamics_different_size(file_obj,'perceptron_200_dynamics.mat')
% get_structure_different_size(file_obj,'perceptron_200_structures.mat')
% 
% file_obj = dir('/scratch/zhang.chi9/perceptron/data/different_network_size_load_at_capacity/fmincon*200*');
% get_retrieval_prob_length(file_obj,noise_added,'fmincon_200_retrieval.mat')
% get_dynamics_different_size(file_obj,'fmincon_200_dynamics.mat')
% get_structure_different_size(file_obj,'fmincon_200_structures.mat')

% file_obj = dir('/scratch/zhang.chi9/perceptron/data/different_network_size_load_at_capacity/perceptron*400*');
% get_retrieval_prob_length(file_obj,noise_added,'perceptron_400_retrieval.mat')
% get_dynamics_different_size(file_obj,'perceptron_400_dynamics.mat')
% get_structure_different_size(file_obj,'perceptron_400_structures.mat')

% file_obj = dir('/scratch/zhang.chi9/perceptron/data/network_load_at_numerical_capacity/rin_-6a_40*');
% %get_retrieval_prob_length(file_obj,noise_added,'fmincon_400_retrieval.mat')
% get_dynamics_different_size(file_obj,'fmincon_400_dynamics.mat')
% get_structure_different_size(file_obj,'fmincon_400_structures.mat')

file_obj = dir('/scratch/zhang.chi9/perceptron/data/different_network_size_load_at_capacity/perceptron*800*');
get_retrieval_prob_length(file_obj,noise_added,'perceptron_800_retrieval.mat')
get_dynamics_different_size(file_obj,'perceptron_800_dynamics.mat')
get_structure_different_size(file_obj,'perceptron_800_structures.mat')

file_obj = dir('/scratch/zhang.chi9/perceptron/data/different_network_size_load_at_capacity/fmincon*800*');
get_retrieval_prob_length(file_obj,noise_added,'fmincon_800_retrieval.mat')
get_dynamics_different_size(file_obj,'fmincon_800_dynamics.mat')
get_structure_different_size(file_obj,'fmincon_800_structures.mat')
