function [quality_indexes,min_angles]=getMeshQuality(faces,nodes)
        quality_indexes=zeros(size(faces,1)*3,1);
        min_angles=zeros(size(faces,1),1);
        for count_face=1:size(faces,1)
                nodel=faces(count_face,:);
                face_nodes=nodes(nodel,:);
                vec1=face_nodes(2,:)-face_nodes(1,:);
                vec1_norm=vec1/norm(vec1);
                vec2=face_nodes(3,:)-face_nodes(2,:);
                vec2_norm=vec2/norm(vec2);
                vec3=face_nodes(3,:)-face_nodes(1,:);
                vec3_norm=vec3/norm(vec3);
                angle1=rad2deg(acos(dot(vec1_norm,vec3_norm)));
                angle2=rad2deg(acos(dot(vec2_norm,-vec1_norm)));
                angle3=180-angle1-angle2;
                quality_indexes(count_face*3-2)=angle1;
                quality_indexes(count_face*3-1)=angle2;
                quality_indexes(count_face*3)=angle3;
                min_angles(count_face)=min([angle1,angle2,angle3]);
        end

end