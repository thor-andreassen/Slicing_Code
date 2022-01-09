function [edges]=findEdgesMesh(faces)
        edges=zeros(size(faces,1)*3,2);
        for count_edge=1:size(faces,1)
                edges(count_edge*3-2,:)=faces(count_edge,[1,2]);
                edges(count_edge*3-1,:)=faces(count_edge,[2,3]);
                edges(count_edge*3,:)=faces(count_edge,[1,3]);
        end
end