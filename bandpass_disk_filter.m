function filtered_image=bandpass_disk_filter(image,r_small,r_big)
%*------------------------------------------------*
% Developmental Plasticity of Sensory System Lab
% University of Science and Technology of China
% Author : Xu An, Aug. 2011
% anxu@mail.ustc.edu.cn
%*------------------------------------------------*
if nargin==1
    low_pass=fspecial('disk',9);
    high_pass=fspecial('disk',33);
    filtered_image=filter2(low_pass,image,'same')-filter2(high_pass,image,'same');
elseif nargin == 2
    low_pass=fspecial('disk',r_small);
    filtered_image=filter2(low_pass,image,'same');
else
    low_pass=fspecial('disk',r_small);
    high_pass=fspecial('disk',r_big);
    filtered_image=filter2(low_pass,image,'same')-filter2(high_pass,image,'same');
end