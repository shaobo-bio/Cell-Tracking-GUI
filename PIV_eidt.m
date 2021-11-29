function pushbutton_start_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ittWidth = 32;
ittHeight = 32;
ovlapHor = 32;
ovlapVer = 32;

s2ntype = 2;
s2nl = 1;
sclt = 1;
outl = 100;

preprocess = 0;
if preprocess
    prepfun = str2func(handles.preprocess);
else
    prepfun = inline('x');
end

if isfield(handles,'rect') && ~isempty(handles.rect)
    cropvec = handles.rect;
else
    cropvec = [0 0 0 0];
end
jump = 1;
if jump == -1
    jump = 0;
end

% Main loop
set(handles.figure1,'pointer','watch')


image1 = fullfile(handles.path,handles.files{1});
image2 = fullfile(handles.path,handles.files{2});
[a,b,a1,b1] = read_pair_of_images_rect(image1,image2,cropvec,ittWidth,ittHeight,ovlapHor,ovlapVer);
if isempty(a) || isempty(b)
    errordlg('Something wrong with your images')
end

switch handles.filesType
    case{'sequence'}
        
        for fileind = 1:handles.amount-jump	% main loop, for whole file list
            %      while (get(handles.pushbutton_start,'UserData') ==1)
            
            %             set(handles.edit_num,'string',sprintf('%d/%d',fileind,handles.amount-jump));
            set(handles.edit_num,'string',sprintf('%d',fileind));
            
            image1 = fullfile(handles.path,handles.files{fileind});
            image2 = fullfile(handles.path,handles.files{fileind+jump});
            
            [a,b,a1,b1,origin] = read_pair_of_images_rect(image1,image2,cropvec,ittWidth,ittHeight,ovlapHor,ovlapVer);
            
            a1 = prepfun(a1);
            b1 = prepfun(b1);
            
            
            [verSize,horSize]= size(a1);
            
            % Prepare the results storage;
            numcols = floor((horSize-ittWidth)/ovlapHor+1);
            numrows = floor((verSize-ittHeight)/ovlapVer+1);
            res = zeros(numcols*numrows,5);
            resind = 0;
            
            a2 = zeros(ittHeight,ittWidth);
            b2 = zeros(ittHeight,ittWidth);
            NfftWidth = 2*ittWidth;
            NfftHeight = 2*ittHeight;
            
            %%%%%% Start the loop for each interrogation block %%%%%%%
            axes(handles.axes1);
            % imshow(imadjust(a),[]);
            imshow(prepfun(a),[]);
            hold on
            
            for m = 1:ovlapVer:verSize - ittHeight + 1 % vertically
                for k = 1:ovlapHor:horSize-ittWidth+1 % horizontally
                    if (get(hObject,'UserData') == 1)
                        a2 = a1(m:m+ittHeight-1,k:k+ittWidth-1);
                        b2 = b1(m:m+ittHeight-1,k:k+ittWidth-1);
                        
                        %                         a2 = prepfun(a2);
                        %                         b2 = prepfun(b2);
                        
                        c = cross_correlate_rect(a2,b2,NfftHeight,NfftWidth);
                        % c = cross_correlate_rect(a2,b2,Nfftx,Nffty);
                        if ~any(c(:)), % completely "black"
                            u = 0;
                            v = 0;
                            y = origin(1) + m + ittHeight/2 - 1;  % y set interests points
                            x = origin(2) + k + ittWidth/2 -  1;  % x
                            continue
                        end
                        
                        [peak1,peak2,pixi,pixj] = find_displacement_rect(c,s2ntype);
                        
                        [peakVer,peakHor,s2n] = sub_pixel_velocity_rect(c,pixi,pixj,peak1,peak2,s2nl,sclt,ittWidth,ittHeight);
                        
                        % Scale the pixel displacement to the velocity
                        u = (ittWidth-peakHor)*sclt;
                        v = (ittHeight-peakVer)*sclt;
                        y = origin(2) + m + ittHeight/2-1;
                        x = origin(1) + k + ittWidth/2-1;
                        
                        resind = resind + 1;
                        res(resind,:) = [x y u v s2n];
                        % quiver(x+cropvec(1),y+cropvec(2),u,v,'y');
                        
                        
                        if u ~= 0 || v ~= 0
                            %                             quiver(x,y,u,v,5,'y','Linewidth',1);
                            %                             drawnow;
                            % plotarrow(x,y,u,v,'y',2);
                            % drawnow
                        end
                        
                        
                    end
                end
            end
            
            
            
            % NO_FILT_RES will be stored in '.._noflt.txt' file at the end of program
            no_filt_res = res;
            
            % Reshape U and V matrices in two-dimensional grid and produce
            % velocity vector in U + i*V form (real and imaginary parts):
            
            u = reshape(res(:,3), numrows,numcols);
            v = reshape(res(:,4), numrows,numcols);
            vector = u + sqrt(-1)*v;
            
            % Remove outlayers - GLOBAL FILTERING
            vector(abs(vector)>mean(abs(vector(find(vector))))*outl) = 0;
            u = real(vector);
            v = imag(vector);
            
            % Adaptive Local Median filtering
            
            kernel = [-1 -1 -1; -1 8 -1; -1 -1 -1];
            tmpv = abs(conv2(v,kernel,'same'));
            tmpu = abs(conv2(u,kernel,'same'));
            
            % WE HAVE TO DECIDE WHICH LIMIT TO USE:
            % 1. Mean + 3*STD for each one separately OR
            % 2. For velocity vector length (and angle)
            % 3. OR OTHER.
            
            lmtv = mean(tmpv(find(tmpv))) + 3*std(tmpv(find(tmpv)));
            lmtu = mean(tmpu(find(tmpu))) + 3*std(tmpu(find(tmpu)));
            u_out = find(tmpu>lmtu);
            v_out = find(tmpv>lmtv);
            
            % Let's throw the outlayers out:
            u(u_out) = 0; u(v_out) = 0;
            v(v_out) = 0; v(u_out) = 0;
            vector = u + sqrt(-1)*v;
            
            res(:,3) = reshape(real(vector),numrows*numcols,1);
            res(:,4) = reshape(imag(vector),numrows*numcols,1);
            
            
            vector = fill_holes(vector,numrows,numcols);
            res(:,3) = reshape(real(vector),numrows*numcols,1);
            res(:,4) = reshape(imag(vector),numrows*numcols,1);
            
            %res(:,3) the u for x
            %res(:,4) the v for y
            
            % res the result!!!
 