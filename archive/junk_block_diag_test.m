% block diagonal test

clear all
clc


G = rand(10,9);
for i = 1:size(G,1)
    for j = 1:size(G,2)
        if i <= 3
            G(i,4:end) = 0;
        elseif i <= 7
            G(i,1:3) = 0;
            G(i,7:end) = 0;
        else
            G(i,1:6) = 0;
        end
    end
end
G(1:3,:) = 0;

Cd = diag(rand(size(G,1),1));

GtCG = G.'*Cd*G
%inv(GtCG)

Cm = diag(diag(GtCG));
%Cm = GtCG;

F = diag(4*ones(size(Cm,1),1));

for i = 1:size(F,1)
    if i + 3 <= size(F,1)
        F(i,i+3) = -1;
    end
    if i - 3 > 0
        F(i,i-3) = -1;
    end
end

Cmm = F.'*Cm*F
L = GtCG + Cmm
inv(GtCG + Cmm)

cnt = 0;
for i = 1:size(L,1)
    for j = 1:size(L,2)
        if L(i,j) ~= 0
            cnt = cnt + 1;
        end
    end
end

cnt
cnt/(size(L,1)*size(L,2))

