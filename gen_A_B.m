%% GEN A with two eigenvalues > 1
while 1
A=2*rand(3)-1;
    if sum(abs(eig(A))>1)>=2
        break
    end
end
B=2*rand(3,1)-1;
disp(eig(A))
