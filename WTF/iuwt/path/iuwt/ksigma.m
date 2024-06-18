%% 截断高斯核函数，截断的程度取决于参数bias  
function ksigma = ksigma(ksize, sigma, center,bias)

% ksize  % 高斯大小   矩阵大小
% sigma   % 方差大小 
% center  = round(ksize/2);          % 中心点
% bias  = ksize*10/10;              % 偏移中心点量

ksigma=fspecial('gaussian',ksize, sigma);   % 构建高斯函数
% 截断高斯核函数，截断的程度取决于参数bias  
[m, n] =size(ksigma);  
for i = 1:m  
    for j = 1:n  
        if(  (i<center-bias)||(i>center+bias)||(j<center-bias)||(j>center+bias)  )  
            ksigma(i,j) = 0;  
        end  
    end  
end 

% figure; mesh(ksigma);
end


