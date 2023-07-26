function [CUT,det_rangeindex,det_veloindex]=os_cfar(RDM)
    %Select the number of Training Cells in both the dimensions.
    pfa=0.1;%10^(-6);
    Tr = 4;
    Tc = 4;
    % *%TODO* :
    %Select the number of Guard Cells in both dimensions around the Cell under 
    %test (CUT) for accurate estimation
    Gc = 2;
    Gr = 2;
    [row,col] = size(RDM);
    % RDM window padding
    RDM_padding=zeros(row+2*(Tr+Gr),col+2*(Gc+Tc));
    RDM_padding(Tr+Gr+1:Tr+Gr+row,Tc+Gc+1:Tc+Gc+col)=RDM;
    CUT = zeros(row,col);
    det_rangeindex=[];
    det_veloindex=[];
    rate=0.75; % hling H．Radar CFAR thresholding in clutter and multi— ple target situationsl J 3/4*M
    % 计算alpha
    d1_scope=[5,15];
    d2_precision=0.5;
    func=@form_PFA_os;
    %%属于二维遍历，下边大循环找出中心点  Cut
    for i = Tr+Gr+1:row+Tr+Gr
        for j = Tc+Gc+1:col+Tc+Gc
            noise_level = [];
            referCellnum=length(find(RDM_padding(i-Tr-Gr:i+Tr+Gr,j-Tc-Gc:j+Tc+Gc)))-length(find(RDM_padding(i-Gr:i+Gr,j-Gc:j+Gc)));  
            parameter=[referCellnum,rate,pfa];
            [ alpha_os ] = oscfar_alpha_binary( d1_scope,d2_precision,func,parameter);
            %%下边循环属于将矩形选出来  Guard cell
            for k = (i-Tr-Gr) : (i+Tr+Gr)
                for h = (j-Tc-Gc) : (j+Gc+Tc)
                    if((abs(k-i)>Gr||abs(h-j)>Gc) && RDM_padding(k,h)~=0) %%选出训练单元
                        noise_level = [noise_level,db2pow(RDM_padding(k,h))];  % sort
                    end
                end
            end
            noise_aftersort=sort(noise_level);
            threshold(i,j) = pow2db( noise_aftersort(ceil(referCellnum.*rate))*alpha_os); % os-cfar
            if(RDM_padding(i,j)>threshold(i,j))
                det_rangeindex=[det_rangeindex,i-Tr-Gr];
                det_veloindex=[det_veloindex,j-Tc-Gc];
                CUT(i-Tr-Gr,j-Tc-Gc) = 1;
            else
                CUT(i-Tr-Gr,j-Tc-Gc) = 0;
            end
        end
    end
end
