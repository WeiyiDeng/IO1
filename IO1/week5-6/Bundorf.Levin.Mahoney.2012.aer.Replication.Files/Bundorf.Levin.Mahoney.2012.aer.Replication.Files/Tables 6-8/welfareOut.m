clear all; clc;

load welfareStats;
load structureWelfareData;

N = length(choiceData);

SIM = 11;

for sim = 1:SIM

    out(sim).noffered = zeros(4,1);  %#ok<AGROW>
    out(sim).marketshare = zeros(4,1); %#ok<AGROW>
    out(sim).risk = zeros(4,1); %#ok<AGROW>
    out(sim).price = zeros(4,1); %#ok<AGROW>
    out(sim).n = zeros(4,1); %#ok<AGROW>
    out(sim).gs = 0; %#ok<AGROW>
    out(sim).cost = 0;  %#ok<AGROW>
    out(sim).ss = 0; %#ok<AGROW>

    for n = 1:N

        if welfare(n).plan(1,sim) ~= 99

            plan = welfare(n).plan(1,sim);

            out(sim).n(plan) = out(sim).n(plan) + welfare(n).fam; %#ok<AGROW>

            out(sim).noffered = out(sim).noffered + welfare(n).offered*welfare(n).fam; %#ok<AGROW>

            out(sim).risk(plan) = out(sim).risk(plan) + welfare(n).pro*welfare(n).fam; %#ok<AGROW>

            out(sim).price = out(sim).price + welfare(n).price(:,sim)*welfare(n).fam; %#ok<AGROW>

            out(sim).gs = out(sim).gs + welfare(n).gs(1,sim)*welfare(n).fam; %#ok<AGROW>

            out(sim).cost = out(sim).cost + welfare(n).cost(1,sim)*welfare(n).fam; %#ok<AGROW>

            out(sim).ss = out(sim).ss + welfare(n).ss(1,sim)*welfare(n).fam; %#ok<AGROW>

        end

    end

end

for sim = 1:SIM

    out(sim).marketshare = out(sim).n/sum(out(sim).n); %#ok<AGROW>

    out(sim).risk = out(sim).risk./out(sim).n;%#ok<AGROW>

    out(sim).price = out(sim).price./out(sim).noffered;%#ok<AGROW>
    out(sim).price = (out(sim).price - out(sim).price(3)*ones(4,1)); %#ok<AGROW>

    out(sim).gs = out(sim).gs/sum(out(sim).n);%#ok<AGROW>

    out(sim).cost = out(sim).cost/sum(out(sim).n);%#ok<AGROW>

    out(sim).ss = out(sim).ss/sum(out(sim).n);%#ok<AGROW>

end

gsBase = out(1).gs;
ssBase = out(1).ss;

for sim = 1:SIM
    
    out(sim).gs = out(sim).gs - gsBase; %#ok<AGROW>
    out(sim).ss = out(sim).ss - ssBase; %#ok<AGROW>
    
end

fid = fopen('welfare.txt','w');

out(1).label = 'Observed';
out(2).label = 'Efficient';
out(3).label = 'Uniform';
out(4).label = 'Uniform by Firm';
out(5).label = 'Uniform by Tier';
out(6).label = 'Uniform by Firm by Tier';
out(7).label = 'AHMO';
out(8).label = 'APPO';
out(9).label = 'BHMO';
out(10).label = 'BPOS';
out(11).label = 'Optimal Single Plan by Firm';

for sim = 1:SIM

    fprintf(fid, '%s\n', out(sim).label);

    fprintf(fid,'%6.3f %6.3f %6.3f %6.3f %6.1f %6.1f %6.1f\n ', [out(sim).marketshare; out(sim).gs; out(sim).cost; out(sim).ss]);

    fprintf(fid,'%6.3f %6.3f %6.3f %6.3f\n ',out(sim).risk);

    fprintf(fid,'%6.1f %6.1f %6.1f %6.1f\n ',out(sim).price);

end

fclose(fid);

