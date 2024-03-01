function vertex_normal=findVertexNormalFromMesh(faces,nodes,node_num)
        %% This function determines the Approximate Normal at a given Vertex
        % Created by Thor Andreasen
        % 1/10/22
        % The function determines the approximate vertex normal as the
        % average of the weighted some of normals of the attached faces
        % based on the area of those faces.
        %
        %
        % Inputs:
                % faces: (n x 3) matrix of connections of nodes defined by
                        % corner node indices for a triangulated mesh
                % nodes: (mx3) matrix of node coordinates (x, y, z) with each row
                        % corresponding to the assumed node number and
                        % corresponding to the nodes used in the "faces" input
                % node_num: a single number correpsonding to the node to
                        % calculate the normal of. this number should
                        % corresponding to the row of the node in the nodes
                        % input
        % Outputs:
                % vertex_normal: the x,y,z normalized normal direction for
                % the vertex
        [faces_with_node,~]=find(faces==node_num);
        total_area=0;
        vertex_normal=zeros(1,3);
        for count_face=1:length(faces_with_node)
                nodel=faces(faces_with_node(count_face),:);
                face_pts=nodes(nodel,:);
                [face_normal]=findTriFaceNormalArea(face_pts);
                vertex_normal=vertex_normal+face_normal;
        end
        vertex_normal=vertex_normal/norm(vertex_normal);

end
