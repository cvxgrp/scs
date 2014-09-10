function z = proj_soc(tt)
if isempty(tt)
    z=[];
    return;
elseif length(tt)==1
    z = max(tt,0);
    return;
end
v1=tt(1);v2=tt(2:end);
if norm(v2)<=-v1
    v2=zeros(length(v2),1);v1=0;
elseif norm(v2)> abs(v1)
    v2=0.5*(1+v1/norm(v2))*v2;
    v1=norm(v2);
end
z=[v1;v2];
end