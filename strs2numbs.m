function numbers = strs2numbs(s)
%*------------------------------------------------*
% Developmental Plasticity of Sensory System Lab
% University of Science and Technology of China
% Author : Xu An, Aug. 2011
% anxu@mail.ustc.edu.cn
%*------------------------------------------------*

s = [s ' '];
n = 1;
itm = '';
bit10 = 1;
for i = 1:length(s)
    if s(i) ~= ' '
        itm(bit10) = s(i);
        bit10 = bit10+1;
    else
        if ~isempty(itm)
            numbers(n) = str2double(itm);
            if isnan(numbers(n))
                numbers = [];
                errordlg('Please check the format of your input values', 'ERROR');
                return;
            end
            n = n+1;
            itm = '';
            bit10 = 1;
        end
    end
end
numbers = sort(numbers);