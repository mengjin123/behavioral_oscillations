function [data_sig,test_threshold_H, test_threshold] = permutation_threshold_calc_1D(data,perm_data,alpha1,freq_correction_tag);

%%%perm_data: permute_size*x, data:x

test = perm_data;
H_threshold = floor(size(perm_data,1)*(1-alpha1));
L_threshold = floor(size(perm_data,1)*alpha1);

for i=1:size(test,2)
    temp=sort(squeeze(test(:,i)));
    test_threshold_H(i)=temp(H_threshold);
    test_threshold_L(i)=temp(L_threshold);
    test_threshold = test_threshold_H;
end

if freq_correction_tag == 1

    test_threshold_H = max(test_threshold_H); 
    test_threshold_L = min(test_threshold_L);%%%%%finding the maxium value in each frequency
end

data_sig = (data >= test_threshold_H);
