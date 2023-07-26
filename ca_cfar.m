% %% CFAR implementation
%Slide Window through the complete Range Doppler Map
function [CUT,det_rangeindex,det_veloindex]=ca_cfar(RDM)
    % *%TODO* :
    %Select the number of Training Cells in both the dimensions.
    pfa=0.1;
    Tr = 4;
    Tc = 4;
    % *%TODO* :
    %Select the number of Guard Cells in both dimensions around the Cell under 
    %test (CUT) for accurate estimation
    Gc = 2;
    Gr = 2;
    length = 2*(Tr+Gr)+1;
    width =  2*(Tc+Gc)+1;
    referCellnum=length*width-(2*Gr+1)*(2*Gc+1);
    [row,col] = size(RDM);
    CUT = zeros(row,col);
    alpha_ca = referCellnum*(pfa^(-1/referCellnum) - 1);
    det_rangeindex=[];
    det_veloindex=[];
    %%注意循环开始的位置，意味着边缘处的门限不考虑
    %%属于二维遍历，下边大循环找出中心点  Cut
    for i = Tr+Gr+1:(row-Tr-Gr)
        for j = Tc+Gc+1:(col-Tc-Gc)
            noise_level = zeros(1,1);
            %%下边循环属于将矩形选出来  Guard cell
            for k = (i-Tr-Gr) : (i+Tr+Gr)
                for h = (j-Tc-Gc) : (j+Gc+Tc)
                    if(abs(k-i)>Gr||abs(h-j)>Gc) %%选出训练单元
                        noise_level = noise_level + db2pow(RDM(k,h));
                    end
                end
            end
            threshold(i,j) = pow2db( noise_level/referCellnum*alpha_ca); % ca-cfar
            if(RDM(i,j)>threshold)
                det_rangeindex=[det_rangeindex,i];
                det_veloindex=[det_veloindex,j];
                CUT(i,j) = 1;
            else
                CUT(i,j) = 0;
            end
        end
    end
end