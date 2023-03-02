function [HrBt,HbrR,HbrBt] = gen_nearFieldChan0516(MtSet,MrSet,MrisSet,d1,d2,d3,f)
% d1: the distance between the center of receive and transmit antennas array
% d2: the distance between the plane of RIS and antennas
% d3: the distance between the center of ris and rx array in x-axis direction
% RIS and antenna space interval is lambda/2
c = 3e8;
lma = c/f;
dAnte = lma/2;
dRisEle = dAnte;
k = 2*pi/lma;
%% array setting
% receive antenna array
rxCenter = zeros(3,1);
rxArray = cal_coord_plane(rxCenter,MrSet,dAnte);
% transmit antenna array
txCenter = rxCenter+[d1;0;0];
txArray = cal_coord_plane(txCenter,MtSet,dAnte);
% ris antenna array
% d3 
risCenter = rxCenter+[d3;0;d2];
risArray = cal_coord_plane(risCenter,MrisSet,dRisEle);
% rx and tx (combined) to ris



%% test plot
% figure;
% myTestPlot(rxArray,'ro','markersize',10,'lineWidth',2);
% myTestPlot(txArray,'bx','markersize',10,'lineWidth',2);
% myTestPlot(risArray,'k^','markersize',10,'lineWidth',2);
% grid on
% xlabel('X (m)')
% ylabel('Y (m)')
% zlabel('Z (m)')
% axis equal
%% channel
% tx to rx
HbrBt = cal_chan(rxArray,txArray);
% tx to ris
HrBt = cal_chan(risArray,txArray);
% ris to rx
HbrR = cal_chan(rxArray,risArray);
% antenna array (both tx and rx) to ris

    function y = cal_coord_plane(center,planeSize,d)
        Mx = planeSize(1);
        My = planeSize(2);
        width = (Mx-1)*d;
        length = (My-1)*d;
        startPoint = center-[length/2;width/2;0];
        y = zeros(Mx,My,3);
        for ii=1:Mx
            for jj=1:My
                y(ii,jj,:) = startPoint+(ii-1)*[0;d;0]+(jj-1)*[d;0;0];
            end
        end 
    end       
        
    function y = cal_chan(array1,array2)
        M11 = size(array1,1);
        M12 = size(array1,2);
        M21 = size(array2,1);
        M22 = size(array2,2);
        Mr = M11*M12;
        Mt = M21*M22;
        array1Re = reshape(array1,Mr,3);
        array2Re = reshape(array2,Mt,3);
        y = zeros(Mr,Mt);
        for ii = 1 : Mr
            for jj = 1 : Mt
                y(ii,jj) = cal_coeff(array1Re(ii,:),array2Re(jj,:));
            end
        end
        
    end
    function y = cal_coeff(xAxis,yAxis)   
        d = sqrt((xAxis(1)-yAxis(1))^2+(xAxis(2)-yAxis(2))^2+(xAxis(3)-yAxis(3))^2);
        y = sqrt(1/4*(1/(k*d)^2-1/(k*d)^4+1/(k*d)^6))*exp(-1j*k*d); % Schantz. [2005]
    end
    function myTestPlot(array,plt1,plt2,width,plt3,plt4)
        [M1,M2,~] = size(array);
        for ii = 1 : M1
            for jj = 1 : M2
                tmp = array(ii,jj,:);
                plot3(tmp(1),tmp(2),tmp(3),plt1,plt2,width,plt3,plt4);
                hold on
            end
        end
        
        
    end
end