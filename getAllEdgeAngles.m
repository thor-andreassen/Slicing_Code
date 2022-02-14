function edge_angles=getAllEdgeAngles(faces,vertices)
        temp_mesh=triangulation(faces,vertices);
        adjacent_faces=neighbors(temp_mesh);
        face_normals=zeros(size(faces,1),3);
        for count_face=1:size(faces,1)
                nodel=faces(count_face,:);
                node_vertices=vertices(nodel,:);
                face_normals(count_face,:)=findTriFaceNormalArea(node_vertices);
        end

        adjacent_angles=zeros(size(faces,1),3);
        for count_face=1:size(faces,1)
                current_face_normal=face_normals(count_face,:);
                for count_edge=1:3
                        adjacent_face=adjacent_faces(count_face,count_edge);
                        adjacent_face_normal=face_normals(adjacent_face,:);
                        current_angle=real(acosd(dot(current_face_normal,adjacent_face_normal)));
                        adjacent_angles(count_face,count_edge)=current_angle;
                end
        end
        edge_angles=reshape(adjacent_angles,[],1);
end