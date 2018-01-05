function res = Dummies(x);

L = length(x);
Sorted_x = sort(x);

Values_Ind = [1; 1+find(Sorted_x(1:L-1)~=Sorted_x(2:L))];
Values = Sorted_x(Values_Ind);

D = length(Values);
res = zeros(L,D);

for i = 1:D,
    res(:,i) = (x==Values(i));
end;






