%% MyCrust
%
%Simple surface recostruction program based on d algorithm
%Given a set of 3D points returns a triangulated tight surface.
%
%The more points there are the best the surface will be fitted,
%although you will have to wait more. For very large models an 
%help memory errors may occurs.
%It is important even the point distribution, generally uniformly 
% distributed points with denser zones in high curvature features
% give the best results.
% 
% Remember  crust algorithom needs a cloud representing a volume 
% so open surface may give inaccurate results. 
%   
%   
%Here is a simple example:
%
%load Dino.mat%load input points from mat file
%
%[t]=MyCrust(p);
%
% figure(1)
%         hold on
%         title('Output Triangulation','fontsize',14)
%         axis equal
%         trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b')
%
%Input:
%              p is a Nx3 array containing the 3D set of points
%Output:
%              t are points id contained in triangles Nx3 array too.

%


function [t]=MyCrust(b)
%%   Main
starttime=clock;


%add points to the given ones, this is usefull
%to create outside tetraedrom
tic
b=AddShield(b);
fprintf('Addedded Shield: %4.4f s\n',toc)


tic
tetr=delaunayn(b);%creating tedraedron
tetr=int32(tetr);%use integer to save memory
fprintf('Delaunay Triangulation Time: %4.4f s\n',toc)


%connectivity data
%find triangles to tetraedrom and tetraedrom to triangles connectivity data
tic
[t2tetr,tetr2t]=Connectivity(tetr);
fprintf('Connectivity Time: %4.4f s\n',toc)
tic
[cc,r]=CC();%Circumcenters of tetraedroms
fprintf('Circumcenters Time: %4.4f s\n',toc)
clear n



tic
t=Walking();%Flagging tetraedroms as inside or outside
fprintf('Walking Time: %4.4f s\n',toc)


time=etime(clock,starttime);
fprintf('Total Time: %4.4f s\n',time)





%% Circumcenters(Nested)
    function [cc,r]=CC()
         %finds circumcenters fro a set of tetraedrom
        
         %points of tetraedrom 
        p1=(b(tetr(:,1),:));
        p2=(b(tetr(:,2),:));
        p3=(b(tetr(:,3),:));
        p4=(b(tetr(:,4),:));

        %vectors of tetraedrom edges
        v21=b(tetr(:,1),:)-b(tetr(:,2),:);
        v31=b(tetr(:,3),:)-b(tetr(:,1),:);
        v41=b(tetr(:,4),:)-b(tetr(:,1),:);

        %preallocation
        cc=zeros(size(tetr,1),3);






         %Solve the system using cramer method
        d1=sum(v41.*(p1+p4)*.5,2);
        d2=sum(v21.*(p1+p2)*.5,2);
        d3=sum(v31.*(p1+p3)*.5,2);

        det23=(v21(:,2).*v31(:,3))-(v21(:,3).*v31(:,2));
        det13=(v21(:,3).*v31(:,1))-(v21(:,1).*v31(:,3));
        det12=(v21(:,1).*v31(:,2))-(v21(:,2).*v31(:,1));

        Det=v41(:,1).*det23+v41(:,2).*det13+v41(:,3).*det12;



        detx=d1.*det23+...
            v41(:,2).*(-(d2.*v31(:,3))+(v21(:,3).*d3))+...
            v41(:,3).*((d2.*v31(:,2))-(v21(:,2).*d3));



        dety=v41(:,1).*((d2.*v31(:,3))-(v21(:,3).*d3))+...
            d1.*det13+...
            v41(:,3).*((d3.*v21(:,1))-(v31(:,1).*d2));



        detz=v41(:,1).*((v21(:,2).*d3)-(d2.*v31(:,2)))...
            +v41(:,2).*(d2.*v31(:,1)-v21(:,1).*d3)...
            +d1.*(det12);

        %Circumcenters
        cc(:,1)=detx./Det;
        cc(:,2)=dety./Det;
        cc(:,3)=detz./Det;


         %Circumradius
        r=((sum((p2-cc).^2,2))).^.5;%quadrato del raggio
    end



%% Connectivity(Nested)

    function [t2tetr,tetr2t]=Connectivity(tetr)

       

        numt = size(tetr,1);
        vect = 1:numt;
        t = [tetr(:,[1,2,3]); tetr(:,[2,3,4]); tetr(:,[1,3,4]);tetr(:,[1,2,4])];%triangles not unique
        [t,j,j] = unique(sort(t,2),'rows');%triangles
        t2tetr = [j(vect), j(vect+numt), j(vect+2*numt),j(vect+3*numt)];%each tetraedrom has 3 triangles


        % triang-to-tetr connectivity
      
        nume = size(t,1);
        tetr2t  = zeros(nume,2,'int32');
        count= ones(nume,1,'int8');
        for k = 1:numt

            for j=1:4
                ce = t2tetr(k,j);
                tetr2t(ce,count(ce)) = k;
                count(ce)=count(ce)+1;
            end

        end


    end      % connectivity()



%% Walking(Nested)
    function t=Walking()
        np=size(b,1)-540;%540 = number of shield points put at the end of array
        numtetr=size(tetr,1);
        nt=size(tetr2t,1);


        deleted=true(numtetr,1);%deleted tetraedroms
        checked=false(numtetr,1);%checked tetraedros 
        onfront=false(nt,1);%tetraedroms that need to be checked
        countchecked=0;
        %First flag as outsde tetraedroms with Shield points
        for i=1:numtetr
            for j=1:4
                if tetr(i,j)>np;
                    deleted(i)=true;
                    checked(i)=true;
                    onfront(t2tetr(i,:))=true;
                    countchecked=countchecked+1;
                    break
                end
            end
        end

        %tollerances to mark as in or out
        toll=zeros(nt,1)+.95; 
        level=0;
        
        alpha=zeros(nt,1);%intersection factor
        %it is computed from radius of the tetraedroms circumscribed sphere
        % and the distance between their center
        for i=1:nt
            if tetr2t(i,2)>0 %jump boundary tetraedrom
                    distcc=sum((cc(tetr2t(i,1),:)-cc(tetr2t(i,2),:)).^2,2);%distance from circumcenters
                    %intersection factor
                    alpha(i)=(-distcc+r(tetr2t(i,1))^2+r(tetr2t(i,2))^2)/(2*r(tetr2t(i,1))*r(tetr2t(i,2)));
            end
        end
        clear cc
        
%         Now we scan all tetraedroms. When one is scanned put on front is
%         neighobur. This means that now even the neighobour can be
%         checked. At the begining only tetraedroms with shield points are
%         on front, because we are sure the are out.
%         Tetraedrom with high intersction factor  will be marked as equal
%         else different.
%         When I say high i mean under a set tollerance that becames lower as the algorithm
%         progresses. This Aims to avoid errors propagation when a tetraedrom is wrong marked.
%         
        
        while countchecked<numtetr && level<50
             level=level+1;%level of scan reached
             
            for id=1:nt%loop trough triangles

                if onfront(id)

                    tetr1=tetr2t(id,1);tetr2=tetr2t(id,2);%tetraedroms linked to triangle under analysis
                    if  tetr2==0 %do not check boundary triangles
                        onfront(id)=false;
                        continue
                       
                    elseif (checked(tetr1) && checked(tetr2)) %tetraedroms are already checked
                         onfront(id)=false;
                        continue
                       
                    end

                    if alpha(id)>=toll(id) %flag as equal
                        if checked(tetr1)%find the checked one between the two
                            deleted(tetr2)=deleted(tetr1) ;%flag as equal
                            checked(tetr2)=true;%check
                            countchecked=countchecked+1;
                            onfront(t2tetr(tetr2,:))=true;%put on front all tetreadrom triangles
                        else
                            deleted(tetr1)=deleted(tetr2) ;%flag as equal
                            checked(tetr1)=true;%check
                            countchecked=countchecked+1;
                            onfront(t2tetr(tetr1,:))=true;%put on front all tetreadrom triangles
                        end
                         onfront(id)=false;%remove from front
                   
                    
                    elseif alpha(id)<-toll(id)%flag as different
                        if checked(tetr1)%find the checked one between the two
                            deleted(tetr2)=~(deleted(tetr1)) ;%flag as different
                            checked(tetr2)=true;%check
                            countchecked=countchecked+1;
                            onfront(t2tetr(tetr2,:))=true;%put on front all tetreadrom triangles
                        else
                            deleted(tetr1)=~(deleted(tetr2)) ;%flag as different
                            checked(tetr1)=true;%check
                            countchecked=countchecked+1;
                            onfront(t2tetr(tetr1,:))=true;%put on front all tetreadrom triangles
                        end
                         onfront(id)=false;%remove from front

                         
                    else
                        toll(id)=toll(id)-.05;%tolleraces were too high next time will be lower

                    end



                end
            end
            
            if level==31 %brute continuation(this may appens when the triangulation is corrupt)
                 warning('Brute continuation necessary')
                onfront(t2tetr(~(checked),:))=true;%force onfront collocation
            end
        end
        
        %delete tetraedroms marked as outside 
        tetr(deleted,:)=[];
        
        
        
        %take boundary triangles from tetraedroms
        %this is the raw surface and needs improvements to be used in CAD
        %systems. Maybe in my next revision I will add surface post treatments.
        %Anyway for grafical purpose this should be good.

        %extract boundary triangles 
        t=BoundTriangles(tetr);
        
        
        %Output Data
        numchecked=countchecked/numtetr;
        if level==50
            warning([num2str(level),' th level was reached\n'])
        else
         fprintf('%4.4f th level was reached\n',level)
        end
        fprintf('%4.4f %% of Tetraedrom were checked\n',numchecked*100)



    end




end

%% AddAddShield
function pnew=AddShield(p)

%find the bounding box
maxx=max(p(:,1));
maxy=max(p(:,2));
maxz=max(p(:,3));
minx=min(p(:,1));
miny=min(p(:,2));
minz=min(p(:,3));

%give offset to the bounding box
step=max(abs([maxx,minx,maxy,miny,maxz,minz]));

maxx=maxx+step;
maxy=maxy+step;
maxz=maxz+step;
minx=minx-step;
miny=miny-step;
minz=minz-step;
step=step/100;
N=10;



%creating a grid lying on the bounding box
vx=linspace(minx,maxx,N);
vy=linspace(miny,maxy,N);
vz=linspace(minz,maxz,N);




[x,y]=meshgrid(vx,vy);
facez1=[x(:),y(:),ones(N*N,1)*maxz];
facez2=[x(:),y(:),ones(N*N,1)*minz];
[x,y]=meshgrid(vy,vz-step);
facex1=[ones(N*N,1)*maxx,x(:),y(:)];
facex2=[ones(N*N,1)*minx,x(:),y(:)];
[x,y]=meshgrid(vx-step,vz);
facey1=[x(:),ones(N*N,1)*maxy,y(:)];
facey2=[x(:),ones(N*N,1)*miny,y(:)];

%add points to the p array
pnew=[p;
 facex1;
 facex2;
 facey1;
 facey2;
 facez1;
 facez2];

% figure(2)
% plot3(pnew(:,1),pnew(:,2),pnew(:,3),'.g')


end



%% BoundTriangles
function t=BoundTriangles(tetr)

numt = size(tetr,1);
vect = 1:numt;
t = [tetr(:,[1,2,3]); tetr(:,[2,3,4]); tetr(:,[1,3,4]);tetr(:,[1,2,4])];
[t,j,j] = unique(sort(t,2),'rows');
t2tetr = [j(vect), j(vect+numt), j(vect+2*numt),j(vect+3*numt)];


% triang-to-tetr connectivity
% Each row has two entries corresponding to the triangle numbers
% associated with each triangle. Boundary triangles have only one tetraedrom
nume = size(t,1);
tetr2t  = zeros(nume,2,'int32');
count= ones(nume,1,'int32');
for k = 1:numt

    for j=1:4
        ce = t2tetr(k,j);
        tetr2t(ce,count(ce)) = k;
        count(ce)=count(ce)+1;
    end

end
tbound=false(numt,1);
tbound=tetr2t(:,2)==0;
t=t(tbound,:);
end


