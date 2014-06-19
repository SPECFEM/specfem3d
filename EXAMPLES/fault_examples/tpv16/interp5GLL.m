data = load('tpv17_input_file.txt');

   ix = data(:,1)+1;
   iy = data(:,2)+1;
   xx = data(:,3);
   yy = data(:,4);
   zz5 = data(:,5);
   zz6 = data(:,6); % initial shear stress
   zz7 = data(:,7);
   zz8 = data(:,8);
   zz9 = data(:,9);
   zz10 = data(:,10);
   zz11 = data(:,11);
   zz12 = data(:,12);
   zz13 = data(:,13);
   zz14 = data(:,14);

   for i=1:length(ix)
       X(ix(i),iy(i)) = xx(i);
       Y(ix(i),iy(i)) = yy(i);
       Z5(ix(i),iy(i)) = zz5(i);
       Z6(ix(i),iy(i)) = zz6(i);
       Z7(ix(i),iy(i)) = zz7(i);
       Z8(ix(i),iy(i)) = zz8(i);
       Z9(ix(i),iy(i)) = zz9(i);
       Z10(ix(i),iy(i)) = zz10(i);
       Z11(ix(i),iy(i)) = zz11(i);
       Z12(ix(i),iy(i)) = zz12(i);
       Z13(ix(i),iy(i)) = zz13(i);
       Z14(ix(i),iy(i)) = zz14(i);
   end

   figure;
   subplot(2,1,1);
   pcolor(X,Y,Z6); shading flat,colorbar,axis image;
   %orient tall, wysiwyg;

   %**** Set here the parameters of the square box domain and mesh : ****
   LX = 48e3;
   LY = 19.5e3;
   NELX = 160;
   NELY = 65;
   P = 4;

   dxe = LX/NELX;
   dye = LY/NELY;
   NEL = NELX*NELY;
   NGLL = P+1; % number of GLL nodes per element

   %[iglob,x,y] = MeshBox(LX,LY,NELX,NELY,NGLL);

   XGLL = GetGLL(NGLL);   % local cordinate of GLL quadrature points

   iglob = zeros(NGLL,NGLL,NELX*NELY);	% local to global index mapping
   nglob = (NELX*(NGLL-1)+1)*(NELY*(NGLL-1)+1);	% number of global nodes

   x     = zeros(nglob,1);		% coordinates of GLL nodes
   y     = zeros(nglob,1);

   e=0;
   last_iglob = 0;
   igL = reshape([1:NGLL*(NGLL-1)],NGLL-1,NGLL);
   igB = reshape([1:NGLL*(NGLL-1)],NGLL,NGLL-1);
   igLB = reshape([1:(NGLL-1)*(NGLL-1)],NGLL-1,NGLL-1);
   xgll = repmat( 0.5*(1+XGLL) , 1,NGLL);
   ygll = dye*xgll';
   xgll = dxe*xgll;

   for ey=1:NELY,
       for ex=1:NELX,
           e = e+1;
           % Take care of redundant nodes at element edges :
           if e==1  % first element: bottom-left
               ig = reshape([1:NGLL*NGLL],NGLL,NGLL);
           else
               if ey==1 	%  elements on first (bottom) row
                   ig(1,:) = iglob(NGLL,:,e-1); 		% left edge
                   ig(2:NGLL,:) = last_iglob + igL; 		% the rest
               elseif ex==1 % elements on first (left) column
                   ig(:,1) = iglob(:,NGLL,e-NELX); 		% bottom edge
                   ig(:,2:NGLL) = last_iglob + igB; 		% the rest
               else 	% other elements
                   ig(1,:) = iglob(NGLL,:,e-1); 		% left edge
                   ig(:,1) = iglob(:,NGLL,e-NELX); 		% bottom edge
                   ig(2:NGLL,2:NGLL) = last_iglob + igLB;
               end
           end
           iglob(:,:,e) = ig;
           last_iglob = ig(NGLL,NGLL);

           % Global coordinates of the computational (GLL) nodes
           x(ig) = dxe*(ex-1)+xgll;
           y(ig) = dye*(ey-1)+ygll;

       end
   end

   [ys I]=sort(y);
   xs = x(I);
   for i=1:length(x)
       XSEM(ix(i),iy(i)) = xs(i);
       YSEM(ix(i),iy(i)) = ys(i);
   end

   ZSEM5 = griddata(X,Y,Z5,XSEM,YSEM,'linear');
   ZSEM6 = griddata(X,Y,Z6,XSEM,YSEM,'linear');
   ZSEM7 = griddata(X,Y,Z7,XSEM,YSEM,'linear');
   ZSEM8 = griddata(X,Y,Z8,XSEM,YSEM,'linear');
   ZSEM9 = griddata(X,Y,Z9,XSEM,YSEM,'linear');
   ZSEM10 = griddata(X,Y,Z10,XSEM,YSEM,'linear');
   ZSEM11 = griddata(X,Y,Z11,XSEM,YSEM,'linear');
   ZSEM12 = griddata(X,Y,Z12,XSEM,YSEM,'linear');
   ZSEM13 = griddata(X,Y,Z13,XSEM,YSEM,'linear');
   ZSEM14 = griddata(X,Y,Z14,XSEM,YSEM,'linear');

   subplot(2,1,2);
   pcolor(XSEM,YSEM,ZSEM6); shading flat,colorbar,axis image;

   file1 = sprintf('tpv17_input_file_interpGLL5.txt');
   fid = fopen(file1,'w');
   for i=1:length(ix)
       fprintf(fid,'%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e\n',...
           ix(i)-1,iy(i)-1,XSEM(ix(i),iy(i)),YSEM(ix(i),iy(i)),ZSEM5(ix(i),iy(i)),ZSEM6(ix(i),iy(i)),...
           ZSEM7(ix(i),iy(i)),ZSEM8(ix(i),iy(i)),ZSEM9(ix(i),iy(i)),ZSEM10(ix(i),iy(i)),ZSEM11(ix(i),iy(i)),...
           ZSEM12(ix(i),iy(i)),ZSEM13(ix(i),iy(i)),ZSEM14(ix(i),iy(i)));  
   end
   fclose(fid);

   dataSEM = load('tpv17_input_file_interpGLL5.txt');
