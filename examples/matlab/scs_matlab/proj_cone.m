function z = proj_cone(z,c)
    z = z + proj_dual_cone(-z,c);
end