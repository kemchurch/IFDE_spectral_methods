function [a_new,epsilon] = desingularize_transform(a,N,u1)
a = reshape(a,[N+1,u1]);
g = a(1,1) + 2*sum(a(2:end,1));
a_new = reshape(a/g,[],1);
epsilon = g;
end