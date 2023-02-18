function angle = angleBetween(v1,v2)
    angle = abs(acos(dot(v1,v2)/(norm(v1)*norm(v2))));
end

