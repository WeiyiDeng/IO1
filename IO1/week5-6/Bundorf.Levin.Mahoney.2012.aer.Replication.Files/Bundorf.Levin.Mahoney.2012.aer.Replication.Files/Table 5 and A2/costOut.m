clear all; clc;

load costOut;

costs(1).string = 'Network Insurer Markup';
costs(1).beta = thetaHat(1);
costs(1).se = se(1);

costs(2).string = 'Integrated Insurer Markup';
costs(2).beta = thetaHat(2);
costs(2).se = se(2);

costs(3).string = 'NHMO';
costs(3).beta = thetaHat(3);
costs(3).se = se(3);

costs(7).string = 'NHMO X (Risk Score - 1)';
costs(7).beta = thetaHat(7);
costs(7).se = se(7);

costs(11).string = 'NHMO X Coinsurance';
costs(11).beta = thetaHat(11);
costs(11).se = se(11);

costs(4).string = 'NPPO';
costs(4).beta = thetaHat(4);
costs(4).se = se(4);

costs(8).string = 'NPPO X  (Risk Score - 1)';
costs(8).beta = thetaHat(8);
costs(8).se = se(8);

costs(12).string = 'NPPO X Coinsurance';
costs(12).beta = thetaHat(12);
costs(12).se = se(12);

costs(5).string = 'IHMO$';
costs(5).beta = thetaHat(5);
costs(5).se = se(5);

costs(13).string = 'IHMO X Coinsurance';
costs(13).beta = thetaHat(13);
costs(13).se = se(13);

costs(9).string = 'IHMO X (Risk Score - 1)';
costs(9).beta = thetaHat(9);
costs(9).se = se(9);

costs(6).string = 'IPOS';
costs(6).beta = thetaHat(6);
costs(6).se = se(6);

costs(10).string = 'IPOS X (Risk Score - 1)';
costs(10).beta = thetaHat(10);
costs(10).se = se(10);

costs(14).string = 'IPOS X Coinsurance';
costs(14).beta = thetaHat(14);
costs(14).se = se(14);

fid = fopen('costOut.txt','w');

fprintf(fid,'\n\n%s\n\n','Costs');

for i = 1:length(costs)

    fprintf(fid,'%s\t',costs(i).string);
    fprintf(fid,'%6.4f\n',costs(i).beta);
    fprintf(fid,'\t\t%6.4f\n',costs(i).se);

end

fclose(fid);


