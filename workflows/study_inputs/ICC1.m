function [ICC1,ICC3,p,ratioR] = ICC1(x,label)
    data = x;
    groupv = label;
    [p,tbl] = anova1(data,groupv,"off");
    ngroup = length(unique(groupv));
    agroups = length(data)/ngroup;
    BETW_v = cell2mat(tbl(2,4));%BETW_n = cell2mat(tbl(2,3));
    WITH_v = cell2mat(tbl(3,4));%WITH_n = cell2mat(tbl(3,3));
    ICC1 = (BETW_v - WITH_v)/(BETW_v + WITH_v*(agroups-1));
%     ICC2(vvv) =  agroups*ICC1(vvv)/(1 + (agroups-1)*ICC1(vvv));
    ICC3 =  (BETW_v - WITH_v)/(BETW_v);
    ratioR = cell2mat(tbl(2,2))/(cell2mat(tbl(2,2)) + cell2mat(tbl(3,2)));
end
