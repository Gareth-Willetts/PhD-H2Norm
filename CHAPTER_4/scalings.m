m = 11;
odd = []; 
even = [];
odd_max = [];
even_max = [];

for i = 3:m+2
    if mod(i,2) == 1
        odd = [odd (i-3)/2];
    else
        even = [even (i/2)-1];
    end
end

for i = m+3:m+7
    if mod(i,2) == 1
        odd_max = [odd_max floor((m-1)/2)];
    else
        even_max = [even_max ceil((m-1)/2)-1];
    end
end

%%
clc
display("Odd");
2*odd

display("Odd max");
2*odd_max
%%
clc
display("Even");
2*even - 1

display("Even max");
2*even_max + 1

%% i = m+4, m+6, m+8...
% Equations (4.60) and (4.64) - even m
clc
m = 10;
i = m+4;
ceil((m-1)/2) - 1
ceil(m/2) - 1

% Equations (4.59) and (4.63) - odd m
m = 11;
i = m+4;
floor((m-1)/2)
floor(m/2)

%% i = m+3
% Equations (4.59) and (4.61)
clc
m = 10;
i = m+3;
2*floor((m-1)/2)
2*((i-3)/2)

%% i = m+3
% Equations (4.60) and (4.62)
clc
m = 9;
i = m+3;
2*(ceil((m-1)/2) - 1) + 1
2*(i/2 - 1) - 1