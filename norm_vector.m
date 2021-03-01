function norm_vector = norm_vector(vector)
norm_vector= (vector - min(vector))/(max(vector)-min(vector));
end
