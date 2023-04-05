%% Return number of rows that have two/one equations

function [two_equations_num, one_equations_num] = mac_count_equations(map)

map(find(map(:)<0))=0;
map(find(map(:)>0))=1;

A_INDEX=1;
B_INDEX=2;
X_INDEX=3;

two_equations_num = length(find(map(:,A_INDEX)==1 & map(:,B_INDEX)==1 & map(:,X_INDEX)==0)) + ...
                    length(find(map(:,A_INDEX)==1 & map(:,B_INDEX)==0 & map(:,X_INDEX)==1)) + ...
                    length(find(map(:,A_INDEX)==0 & map(:,B_INDEX)==1 & map(:,X_INDEX)==1)) + ...
                    length(find(map(:,A_INDEX)==1 & map(:,B_INDEX)==1 & map(:,X_INDEX)==1));
one_equations_num = length(find(map(:,A_INDEX)==1 & map(:,B_INDEX)==0 & map(:,X_INDEX)==0)) + ...
                    length(find(map(:,A_INDEX)==0 & map(:,B_INDEX)==1 & map(:,X_INDEX)==0)) + ...
                    length(find(map(:,A_INDEX)==0 & map(:,B_INDEX)==0 & map(:,X_INDEX)==1));