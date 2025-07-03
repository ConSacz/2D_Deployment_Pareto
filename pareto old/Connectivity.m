function check = Connectivity(pop,Rc)
numberNodes = numel(pop)/2;
adj_matrix = zeros(numberNodes,numberNodes);
for i=1:numberNodes
    for j=1:numberNodes
        if (((pop(2*i-1)-pop(2*j-1))^2+(pop(2*i)-pop(2*j))^2)<=(Rc)^2)
            adj_matrix(i,j) = 1;
        end
    end
end
for i=1:numberNodes
    adj_matrix(i,i) = 0;
end
G= graph(adj_matrix);
v = dfsearch(G,1);
if (size(v,1)==numberNodes)
    check = 1;
else
    check = 0;
end