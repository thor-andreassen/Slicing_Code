function vertex_normal=findVertexNormalFromMesh(faces,nodes,node_num,faces_with_node)
        
        [faces_with_node,~]=find(faces==node_num);
        total_area=0;
        vertex_normal=zeros(1,3);
        for count_face=1:length(faces_with_node)
                nodel=faces(faces_with_node(count_face),:);
                face_pts=nodes(nodel,:);
                [face_normal,face_area]=findTriFaceNormalArea(face_pts);
                total_area=total_area+face_area;
                vertex_normal=vertex_normal+total_area*face_normal;
        end
        vertex_normal=vertex_normal/total_area;
        vertex_normal=vertex_normal/norm(vertex_normal);

end