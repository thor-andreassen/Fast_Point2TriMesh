function [distances,project_pts,outside]=fastPoint2TriMesh(inputs,pts,use_parallel)
    %% Main function to Determine nearest points between an arbitrary point and triangulated surface
    % This code takes in a triangulated surface of faces and nodes/vertices
    % such as those present in STLs. It then takes a set of arbitrary
    % points in space and determines the nearest point on the triangulated
    % surface to this point, as well as the distance away. This code is
    % similar to an existing algoirthm by Daniel Frisch called
    % "point2TriMesh", but is written more efficiently and provides
    % significant speed increases.
    
    % For example, to determine the projection
    % of 10,000 points to the included example surface, the time in
    % Point2TriMesh is approximately 100 seconds. Using this implementation
    % the exact same points can be obtained in approximately 0.3 seconds
    % without parallel computing and 0.18 seconds with parallel computing.
    % Additionally, the proposed algorithm allows values to be stored
    % allowing for subsequent computations to take less time.
    
    % The algorithm works by either taking in or calculating the circle
    % incenter of every triangle, as well as the normal direction to each
    % face or calculating it. Additionally, it either takes in or
    % calculates an initial KDTree for the original face incenters. The
    % code then determines the nearest average face to every one of the
    % inputted test/query points using a KNN Search and the previously
    % created KDTree. These points are used to determine the initial
    % distance and projection. The code then checks every projection point,
    % and first ensures that it is inside the given triangle. If the
    % calcualted projection point is outside of the triangle, the code
    % determines the minimum distance to each of the face's 3 edges to the
    % arbitrary point and then uses the minimum distance and correpsonding
    % point as the projection. This is repeated for every point to get the
    % distance, and projection point for every query point. The face
    % normals are used to determine a positive or negative sign depending
    % on whether the point is nearest the positive normal direction
    % (outside) or the negative normal direction (inside)
    
    % INPUTS:
    % inputs.faces = (N x 3) the faces of the original triangulated surface (Required)
    % inputs.nodes = (M x 3) the nodes or vertices of the original triangulated surface (Required)
    % inputs.face_mean_nodes = (N x 3) the location of the incenter of each
    %                     triangle. This can either be calculated using "getTriInCenter.m" ahead
    %                     of time, or calculated herein.
    % inputs.face_normals = (N x 3) the unit normal direction to each face.
    %                     This can either be calculated using "getFaceCenterAndNormals.m" ahead
    %                     of time, or calculated herein.
    % inputs.tree_model = (model) the KDTree trained to each triangle
    %                     incenter using "KDTreeSearcher.m" This can either be calculated using "KDTreeSearcher.m"
    %                     ahead of time, or calculated herein.
    % pts = (Q x 3) Arbitrary points that we are trying to project onto the
    %                     surface (Required)
    % use_parallel = (0 or 1) A binary that determines whether to use (1) the parallel computing
    %                      or not use (0) the parallel computing. NOTE this requires the parallel computing toolbox(Required)
    
    % OUTPUTS:
    % distances= (Q x 1) The signed distance for each arbitary point (pts)
    %                     to the nearest point on the triangulated surface. Positive means the
    %                     point is off of the positive face normal, where as negative means the
    %                     point is nearest the opposite direction to the
    %                     face normal.
    % project_pts= (Q x 3) The location of the projected point for each arbitary point (pts)
    %                     to the nearest point on the triangulated surface.
    % outside= (Q x 1) A boolean array representing if the given point is on the positive direction
    %                     (1) or the negative direction (0) of the corresponding surface
    
    
    % Written by Thor Andreassen
    % University of Denver
    % 5/9/2023
    
    
    %% CODE:
    % determine proper inputs or calculate values
    faces=inputs.faces;
    nodes=inputs.nodes;
    get_norm=0;
    if isfield(inputs,'face_mean_nodes')
        face_mean_nodes=inputs.face_mean_nodes;
    else
        get_norm=1;
    end
    if isfield(inputs,'face_normals');
        face_normals=inputs.face_normals;
    else
        get_norm=1;
    end
    
    get_tree=0;
    if isfield(inputs,'tree_model')
        tree_model=inputs.tree_model;
    else
        get_tree=1;
    end
    
    if get_norm==1
        [face_mean_nodes,face_normals]=getFaceCenterAndNormals(faces,nodes);
    end
    
    if get_tree==1
        tree_model=KDTreeSearcher(face_mean_nodes);
    end
    
    % determine the location of the nearest face to each querry point
    near_id=knnsearch(tree_model,pts);
    
    % determine the distance between the point and the querry point
    distances=dot(pts-face_mean_nodes(near_id,:),face_normals(near_id,:),2);
    
    % determine the signed direction of the original point location
    signs=sign(distances);
    
    % determine projection point
    direction_vector=face_normals(near_id,:);
    project_pts=pts-(distances.*direction_vector);
    
    if use_parallel==1
        parfor count_pt=1:size(pts,1)
            % check if the point is within the triangle using the
            % barycenter
            dist=[0,0,0];
            q=zeros(3,3);
            project_pt_cur=project_pts(count_pt,:);
            node_ids=faces(near_id(count_pt),:);
            cor1=nodes(node_ids(1),:);
            cor2=nodes(node_ids(2),:);
            cor3=nodes(node_ids(3),:);
            cor4=project_pt_cur;
            check=isInsideTriangle( cor1,cor2,cor3,cor4);
            if check==0
                % the projected point is not within the triangle so the
                % nearest point will lie on the nearest line segment. Check
                % nearest point to each line segment/edge of face and use
                % minimum distance as the correct projection point
                [dist(1),q(1,:)] = projectPointToLineSegment(cor1,cor2,cor4);
                [dist(2),q(2,:)] = projectPointToLineSegment(cor2,cor3,cor4);
                [dist(3),q(3,:)] = projectPointToLineSegment(cor1,cor3,cor4);
                [~,id]=min(dist);
                project_pts(count_pt,:)=q(id,:);
            end
        end

    else
        for count_pt=1:size(pts,1)
            project_pt_cur=project_pts(count_pt,:);
            node_ids=faces(near_id(count_pt),:);
            cor1=nodes(node_ids(1),:);
            cor2=nodes(node_ids(2),:);
            cor3=nodes(node_ids(3),:);
            cor4=project_pt_cur;
            check=isInsideTriangle( cor1,cor2,cor3,cor4);
            if check==0
                [dist(1),q(1,:)] = projectPointToLineSegment(cor1,cor2,cor4);
                [dist(2),q(2,:)] = projectPointToLineSegment(cor2,cor3,cor4);
                [dist(3),q(3,:)] = projectPointToLineSegment(cor1,cor3,cor4);
                [~,id]=min(dist);
                project_pts(count_pt,:)=q(id,:);
            end
        end
    end
    % recalculate the distances based on the final projection points and
    % the correct sign value
    distances=vecnorm(pts-project_pts,2,2).*signs;
    
    % return the sign of each point
    outside=signs>=0;
end





%% helper functions
function [ores] = isInsideTriangle( A,B,C,P)
    v0 = C - A;
    v1 = B - A;
    v2 = P - A;

    dot00 = dot(v0, v0);
    dot01 = dot(v0, v1);
    dot02 = dot(v0, v2);
    dot11 = dot(v1, v1);
    dot12 = dot(v1, v2);

    invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    ores =(u >= 0) && (v >= 0) && (u + v < 1);
end


function [dist,q] = projectPointToLineSegment(A,B,p)
    AB = (B-A);
    AB_squared = dot(AB,AB);
    if(AB_squared == 0)
        q = A;
    else
        Ap = (p-A);
        t = dot(Ap,AB)/AB_squared;
    end
    if (t < 0.0)
        q = A;
    elseif (t > 1.0)
        q = B;
    else
        q = A + t * AB;
    end
    dist=norm(p-q);
end
