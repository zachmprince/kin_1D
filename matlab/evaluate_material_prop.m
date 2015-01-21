function values=evaluate_material_prop(mat_prop,time,x)

values = mat_prop.ft(time) * mat_prop.fx(x);

return
end
