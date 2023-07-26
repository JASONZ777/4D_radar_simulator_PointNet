function mat_to_pts( cloudpoints, filename)
npoints= size( cloudpoints.x,2);
fid = fopen(filename, 'w');
% zero-padding or random sample to 50 samples
if npoints<60
    cloudpoints.x=[cloudpoints.x,zeros(1,60-npoints)];
    cloudpoints.y=[cloudpoints.y,zeros(1,60-npoints)];
    cloudpoints.z=[cloudpoints.z,zeros(1,60-npoints)];
    cloudpoints.velocity=[cloudpoints.velocity,zeros(1,60-npoints)];
else
    id=randsample(1:npoints,60);
    cloudpoints.x=cloudpoints.x(id);
    cloudpoints.y=cloudpoints.y(id);
    cloudpoints.z=cloudpoints.z(id);
    cloudpoints.velocity=cloudpoints.velocity(id);
end
for j = 1 : 60
    fprintf(fid,'%f %f %f %f\n', cloudpoints.x(j), cloudpoints.y(j),cloudpoints.z(j),cloudpoints.velocity(j));                 
end
fclose(fid);
end