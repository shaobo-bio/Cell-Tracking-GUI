function varargout = Trakcing(varargin)
% TRAKCING MATLAB code for Trakcing.fig
%      TRAKCING, by itself, creates a new TRAKCING or raises the existing
%      singleton*.
%
%      H = TRAKCING returns the handle to a new TRAKCING or the handle to
%      the existing singleton*.
%
%      TRAKCING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRAKCING.M with the given input arguments.
%
%      TRAKCING('Property','Value',...) creates a new TRAKCING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Trakcing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Trakcing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Trakcing

% Last Modified by GUIDE v2.5 05-May-2014 14:24:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Trakcing_OpeningFcn, ...
    'gui_OutputFcn',  @Trakcing_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Trakcing is made visible.
function Trakcing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Trakcing (see VARARGIN)

% Choose default command line output for Trakcing

warning off all;
handles.TrackingWay=1;
handles.Thred=0;
handles.c=[];
handles.tmp=[];
handles.Res=[];
handles.YD=[];
handles.Current=0;
handles.p=[];
handles.OriginalC=[];
set(gcf,'menubar','figure');
set(gcf,'ToolBar','figure');

set(handles.slider4,'Value',0.95);
set(handles.edit3,'String','');
set(handles.edit5,'String','');
set(handles.edit1,'String','');
set(handles.edit2,'String','');
set(handles.edit4,'String','');

im=imread('Opening.png');
imshow(im);
axes(handles.axes1);

imshow(im,[]);
handles.Method=1;
handles.Thred2=0.95;
handles.output = hObject;
handles.AnDone=0;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Trakcing wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Trakcing_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Add_Point.
function Add_Point_Callback(hObject, eventdata, handles)
% hObject    handle to Add_Point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
but=1;
i=1;
L=length(handles.tmp);
set(handles.edit4,'String',num2str(L));

while but==1
    i=1+i;
    [x y but]=ginput(1);
    if but==1
        set(handles.edit1,'String',num2str(x));
        set(handles.edit2,'String',num2str(y));
        
        x=round(x);
        y=round(y);
        Pi=plot(x,y,'r+','MarkerSize',4);
        
        handles.tmp=[handles.tmp;x,y];
        Li=length(handles.P);
        handles.P(Li+1)=Pi;
        L=length(handles.tmp);
        set(handles.edit4,'String',num2str(L));
    end
end



handles.output = hObject;

%Update handles structure
guidata(hObject, handles);




function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Quit.
function Quit_Callback(hObject, eventdata, handles)
% hObject    handle to Quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);

% --- Executes on button press in New_Project.
function New_Project_Callback(hObject, eventdata, handles)
% hObject    handle to New_Project (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.c=[];

handles.Flag=1;

handles.AnaTime=1;
Method=1;

[filename,pathname]=uigetfile({'*.tif';'*.h5'},'Choose images');

Filetype=finfo(filename);

if isequal(filename,0)
    disp('Cancled');
else
    text2=get(handles.edit6,'String');
    if isempty(get(handles.edit6,'String'))
        Start=1;
    else
        Start=str2double(text2);
    end
    
    str=[pathname filename];
    
    if strcmp(Filetype,'im')
        Type=1;     % for original images
        [OImg,Slice,TimeZ,Channel] = ULoadTiff_Z(str);
        
        if ndims(OImg)==4
            if Channel == 0
                handles.Flag=0;
                for k=1:TimeZ
                    ImgNew(:,:,k)=max(OImg(:,:,:,k),[],3);
                end
            else
                ImgNew=OImg(:,:,:,1);
                
            end
        else
            ImgNew=OImg;
            OImg=[];
        end
        
        PredictImg=[];
        im=ImgNew(:,:,Start);
        %im = adapthisteq(im);
        
        [centers,~,~] = imfindcircles(im,[2 4],'Method','TwoStage','Sensitivity',0.95);
    elseif strcmp(Filetype,'h5')
        Type=2;
        Method=3;% for h5 file which has processed the original images
        Info = h5info(str);
        Mask = permute(h5read(str,'/volume/data'),[4 3 2 1]);
        ImgNew=Mask(:,:,:,1);
        OImg=[];
        im=ImgNew(:,:,Start);
        Mask = permute(h5read(str,'/volume/prediction'),[4 3 2 1]);
        PredictImg=Mask(:,:,:,2);
        centers=h5findcenter(PredictImg(:,:,Start),0); % h5findcenter is a function used for finding centers of processed images
    end
    
    axes(handles.axes1);
    imshow(im,[]);
end


centers=round(centers);
hold on;

P=[];
for i=1:length(centers)
    
    P(i)=plot(centers(i,1),centers(i,2),'g+','MarkerSize',4);
    
end
%
% hSP = imscrollpanel(hFig,hIm);
%         set(hSP,'Units','normalized',...
%                 'Position',[0 .1 1 .9])

L=length(centers);
set(handles.edit4,'String',num2str(L));
set(handles.edit3,'String',0);
set(handles.edit5,'String',1);

%% greater contract images
% K=mat2gray(ImgNew);
% K=K.*7;
% K(find(K>1))=1;
% ImgNew=K;

handles.Method=Method;
handles.PredictImg=PredictImg;
handles.Type=Type;
handles.P=P;
handles.Size=length(ImgNew(1,1,:));
handles.c=[];
handles.Img=ImgNew;
handles.OrigImg=OImg;
handles.tmp=centers;
handles.ImgIn=Start;
handles.Current=1;

handles.Start=Start;

% Choose default command line output for Tracking_Point
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);




function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Next_Img.
function Next_Img_Callback(hObject, eventdata, handles)
% hObject    handle to Next_Img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Check=handles.ImgIn+handles.TrackingWay*handles.Current;
if  Check== 0 || Check==handles.Size
    set(handles.text6,'String','The Last Image!!!!!');
else
    set(handles.text6,'String','');
    Current=handles.Current;
    handles.Current=handles.Current+1;
    A=handles.c;
    B=handles.tmp;
    l=max([length(A),length(B)]);
    [~,len]=size(handles.c);
    if len < 2*Current
        
        handles.c=[padarray(A,[l-length(A) 0],'post') padarray(B,[l-length(B) 0],'post')];
        [handles.tmp,handles.P]=ShowImg(hObject, eventdata, handles,0);
    elseif len==2*Current
        II=Current;
        A=padarray(A,[l-length(A) 0],'post');
        B=padarray(B,[l-length(B) 0],'post');
        A(:,2*II-1:2*II)=B;
        handles.c=A;
        [handles.tmp,handles.P]=ShowImg(hObject, eventdata, handles,0);
        
    else
        II=Current;
        A=padarray(A,[l-length(A) 0],'post');
        B=padarray(B,[l-length(B) 0],'post');
        A(:,2*II-1:2*II)=B;
        handles.c=A;
        [handles.tmp,handles.P]=ShowImg(hObject, eventdata, handles,1);
    end
end
% when click next save the previous information into handles.c

% Choose default command line output for Tracking_Point
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in Previous_Img.
function Previous_Img_Callback(hObject, eventdata, handles)
% hObject    handle to Previous_Img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.Current==1
    set(handles.text6,'String','The First Image!!!!!');
else
    set(handles.text6,'String','');
    Current=handles.Current;
    handles.Current=handles.Current-1;
    A=handles.c;
    B=handles.tmp;
    l=max([length(A),length(B)]);
    [~,len]=size(handles.c);
    if len < 2*Current
        handles.c=[padarray(A,[l-length(A) 0],'post') padarray(B,[l-length(B) 0],'post')];
        [handles.tmp,handles.P]=ShowImg(hObject, eventdata, handles,1);
    elseif len==2*Current
        II=Current;
        A=padarray(A,[l-length(A) 0],'post');
        B=padarray(B,[l-length(B) 0],'post');
        A(:,2*II-1:2*II)=B;
        handles.c=A;
        [handles.tmp,handles.P]=ShowImg(hObject, eventdata, handles,1);
    else
        II=Current;
        A=padarray(A,[l-length(A) 0],'post');
        B=padarray(B,[l-length(B) 0],'post');
        A(:,2*II-1:2*II)=B;
        handles.c=A;
        [handles.tmp,handles.P]=ShowImg(hObject, eventdata, handles,1);
    end
    
    
end
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Thred=get(hObject,'Value');
handles.Thred=Thred;
[handles.tmp,handles.P]=ShowImg(hObject, eventdata, handles,0);

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Analysis.
function Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Current=handles.Current;
A=handles.c; B=handles.tmp;
l=max([length(A),length(B)]);
[~,len]=size(handles.c);
if len < 2*Current
    handles.c=[padarray(A,[l-length(A) 0],'post') padarray(B,[l-length(B) 0],'post')];
else
    II=Current;
    A=padarray(A,[l-length(A) 0],'post');
    B=padarray(B,[l-length(B) 0],'post');
    A(:,2*II-1:2*II)=B;
    handles.c=A;
end

% save the last image into handles.c

% load centers.mat;
% handles.c=c;

% Res saved the centers position!!
% YD saved the Y distance of cells in a pair!!
handles.c(all(handles.c==0,2),:)=[];
Res=[];
Res=handles.c(:,1:2);
Res=Res(find(Res(:,1)>0),:);
Res=sortrows(Res,1);
Res=sortrows(Res,2);


Ww=length(handles.c(1,:))/2;

hwait=waitbar(0,'Please wait...>>>>>>>>');

for i=1:Ww-1
 
    Tmp=handles.c(:,2*i+1:2*i+2);
    Tmp=Tmp(find(Tmp(:,1)>0),:);
    Res=FindTrack(Res,Tmp,i,handles);
    if Ww-i<3
        waitbar(i/(Ww-1),hwait,'Almost Finished...');
    else
        PerStr=fix(100*i/(Ww-1));
        str=['Please wait...',num2str(PerStr),'%'];
        waitbar(i/(Ww-1),hwait,str);
    end
end

close(hwait);

% define the pairs of cells

x=[];y=[];
[L,Ww]=size(Res);

W=2;

Define=GUI_AnaInterface;

if strcmp(Define,'Yes')
    if handles.Current ~= 1
        handles.Current=1;
        [handles.tmp,handles.P]=ShowImg(hObject, eventdata, handles,1);
    end
    Sure='No';
    while strcmp(Sure,'No')
        delete(findobj(handles.axes1,'Color','m'));
        x=[];y=[];
        i=1;but=1;
        while but==1
            
            [a,b,but]=ginput(1);
            if but==1
                x(i)=b;y(i)=a;
                if i==1
                    plot(a,b,'m*');
                else
                    plot(y,x,'m--');
                end
                
                i=i+1;
            end
        end
        
        Sure=GUI_SureInterface;
        
    end
else
    
    % Res=sortrows(Res,W);
    for i=1:L/2
        x=[x;(Res(2*i,W)+Res(2*i-1,W))/2];
        y=[y;(Res(2*i,W-1)+Res(2*i-1,W-1))/2];
    end
end

P=polyfit(x,y,1);

Z=Res(:,W)*P(1)-Res(:,W-1)+P(2);
R1=Res(find(Z<0),:);
R2=Res(find(Z>0),:);

R1=sortrows(R1,2);
R2=sortrows(R2,2);
P=polyfit(x,y,1);
R1(:,Ww+1)=R1(:,1)+R1(:,2)/P(1);

R2(:,Ww+1)=R2(:,1)+R2(:,2)/P(1);

W1=size(R1,1);
W2=size(R2,1);

% [M,I]=max([W1,W2]);
Res=[];
List=[];
% for i=1:W1
% %     if I==1
%         Res(2*i-1,:)=R1(i,:);
%         [~,Index]=min(abs(R1(i,Ww+1)-R2(:,Ww+1)));
%         Res(2*i,:)=R2(Index,:);
% %     else
% %         Res(2*i-1,:)=R2(i,:);
% %         [~,Index]=min(abs(R2(i,Ww+1)-R1(:,Ww+1)));
% %         Res(2*i,:)=R1(Index,:);
% %     end
%     List=[List;Index];
% end
% List=unique(List);
% R2(List,:)=[];
% if ~isempty(R2)
%     for j=1:length(R2(:,1))
%         [~,Index]=min(abs(R2(j,Ww+1)-R1(:,Ww+1)));
%         Res(2*(W1+j)-1,:)=R1(Index,:);
%         Res(2*(W1+j),:)=R2(j,:);
%     end
% end
for i=1:W1
    RTmp(i,1)=i;
    [~,Index]=min(abs(R1(i,Ww+1)-R2(:,Ww+1)));
    RTmp(i,2)=Index;
    List=[List;Index];
end
List=unique(List);
for i=1:W2
    R2T(i,1)=R2(i,Ww+1);
    R2T(i,2)=i;
end
R2T(List,:)=[];
if ~isempty(R2T)
    for j=1:length(R2T(:,1))
        IndexJ=R2T(j,2);
        RTmp(W1+j,2)=IndexJ;
        [~,Index]=min(abs(R1(:,Ww+1)-R2T(j,1)));
        RTmp(W1+j,1)=Index;
    end
end

RTmp=sortrows(RTmp,1);

for i=1:length(RTmp(:,1))
    Res(2*i-1,:)=R1(RTmp(i,1),:);
    Res(2*i,:)=R2(RTmp(i,2),:);
end


Res(:,Ww+1)=[];
% L=length(Res);
% R1=Res(1:2:L,:);
% R2=Res(2:2:L,:);
% R1=sortrows(R1,2);
% R2=sortrows(R2,2);
% 
% for i=1:L/2
%     Res(2*i-1,:)=R1(i,:);
%     Res(2*i,:)=R2(i,:);
% end
    
% Reorder the position to define the matched cells


[YD,PFit]=FindD2C(Res);
% YD=FindYD(Res);

handles.AnaTime=0;
handles.Result=Res;

% save('YD.mat','YD');
% save('Res.mat','Res');
% save('handles.mat','handles');

handles.OriginalC=handles.c;
handles.Res=Res;

handles.c=unique(handles.c,'rows');
handles.YD=YD;
[handles.tmp,handles.P]=ShowImg(hObject, eventdata, handles,1);
set(handles.text21,'String','green: correct points.') 
set(handles.text23,'String','red: not used points.')
set(handles.text22,'String','yellow: added point.');


GUI_Analysis(YD,Res,handles.Img,handles.ImgIn,handles.OrigImg,handles.Flag,PFit,handles.TrackingWay,handles.Type,handles.OriginalC);


handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


function Result=FindTrack(Res,Tmp,Index,handles)
Im=2*Index;
TmpS=Res(:,Im-1:Im);
L=length(TmpS);
RT=[];
OIn=2*(Index-1);
Way=handles.TrackingWay;
Index=(Index-1)*Way;

ProbePos=[];
img1=handles.Img(:,:,Index+handles.ImgIn);
img2=handles.Img(:,:,Index+handles.ImgIn+Way);

for i=1:length(TmpS)
    TX=TmpS(i,1);
    TY=TmpS(i,2);
    
    ProbeX=linspace(TX-2,TX+2,3);
    ProbeY=linspace(TY-2,TY+2,3);
    [X,Y]=meshgrid(ProbeX,ProbeY);
    ProbePosT = [reshape(X,[],1),reshape(Y,[],1)];
    
    ProbePos=[ProbePos;ProbePosT];
end

[PredictPos,PredictVel,Credibility]=pivTrack(img1,img2,ProbePos,30);
TmpPre=[0,0];

for i=L:-1:1
    
    ProbeSquare=PredictVel(9*i-8:9*i,:);

    
    VPS=median(ProbeSquare);
    V=std(ProbeSquare);
    Wrong=find(abs(ProbeSquare(:,1)-VPS(1))>2*abs(V(1))); 
    Wrong=[Wrong;find(abs(ProbeSquare(:,2)-VPS(2))>2*abs(V(2)))];
    Wrong=unique(Wrong);
%     Wrong
    ProbeSquare(Wrong,:)=[];
    
    VPS=mean(ProbeSquare);
    if i < 20
        Thred=20;
    else
        Thred=25;
    end
    if abs(VPS(1))>Thred
        VPS=[0 0];
    end
    
%     if mod(i,2)==0 && VPS(1) > 0
%         VPS=[0 0];
%     end
        
    LX=TmpS(i,1)+VPS(1);
    LY=TmpS(i,2)+VPS(2);
    
    
    if i < L-1
        TmpPre=RT(L-i-1,:);
    end
    
    %     In=Findmin(Tmp,LX,LY,TmpPre);
    
    %     LX=PredictPos(i,1);
    %     LY=PredictPos(i,2);
    
    if isempty(Tmp)
        RT=[round([LX,LY]);RT];
    else
        [D In]=min((Tmp(:,1)-LX).^2+(Tmp(:,2)-LY).^2);
        
        if D < 1600
            RT=[Tmp(In,:);RT];
            Tmp(In,:)=[];
            
        else
            RT=[round([LX,LY]);RT];
            
        end
    end
end




% for i=0:(L/2-1)
%     %%%%%%% By perivious speed
%     Ind1=L/2-i;
%     Ind2=L/2+1+i;
% %     LX=TmpS(Ind1,1);%+SpeedX;
% %     LY=TmpS(Ind1,2);%+SpeedY;
%
%     LX=PredictPos(Ind1,1);
%     LY=PredictPos(Ind1,2);
%
%     In=Findmin(Tmp,LX,LY,Tmp1);
% %     [D In]=min((Tmp(:,1)-LX).^2+2*(Tmp(:,2)-LY).^2+(Tmp(:,2)-Tmp1(2)).^2);
%     Tmp1=Tmp(In,:);
%
%     RT=[Tmp(In,:);RT];
%     Tmp(In,:)=[];
%
% %     LX=TmpS(Ind2,1);%+SpeedX;
% %     LY=TmpS(Ind2,2);%+SpeedY;
%
%     LX=PredictPos(Ind2,1);
%     LY=PredictPos(Ind2,2);
%     In=Findmin(Tmp,LX,LY,Tmp2);
%
% %     [D In]=min((Tmp(:,1)-LX).^2+3*(Tmp(:,2)-LY).^2+(Tmp(:,2)-Tmp2(2)).^2);
%     Tmp2=Tmp(In,:);
%     RT=[RT;Tmp(In,:)];
%     Tmp(In,:)=[];
% end


Result=[Res,RT];

function In=Findmin(Tmp,LX,LY,TmpP)
% Find the tracking Point
[D In]=min((Tmp(:,1)-LX).^2+2*(Tmp(:,2)-LY).^2+(Tmp(:,2)-TmpP(2)).^2);



function YD=FindYD(Res)
%caculate the Y distance between the two matched cell

W=length(Res(1,:))/2;
L=length(Res(:,1))/2;
YDi=[];
for i=1:L
    T=[];
    for j=1:W
        X=Res(2*i,2*j)-Res(2*i-1,2*j);
        T=[T;X];
    end
    YDi=[YDi,T];
end
YD=YDi;

%%%%%%%%%%%%
function [D2Central,PFit]=FindD2C(Res)
% caculate the distance between the paired cells according to the central
% line

W=length(Res(1,:))/2;
L=length(Res(:,1))/2;
D2Central=[];
PFit=[];
for j=1:W
    x=[];
    y=[];
    
    for i=1:L
        x=[x;(Res(2*i,2*j)+Res(2*i-1,2*j))/2];
        y=[y;(Res(2*i,2*j-1)+Res(2*i-1,2*j-1))/2];
    end
    
    P=polyfit(x,y,1);
    
%  for i=L:-1:1
    for i=1:L
        a=Res(2*i,2*j)/P(1)+Res(2*i,2*j-1);
        b=Res(2*i-1,2*j)/P(1)+Res(2*i-1,2*j-1);
        
        D=(a-b)*P(1)/sqrt(P(1)^2+1);
%         D=abs(Res(2*i,2*j)-Res(2*i-1,2*j));
        D2Central(j,i)=D;
    end
    PFit=[PFit;P];
    
end




function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=get(hObject,'String');
index=str2double(index);
[~,len]=size(handles.c);
len=len/2;

if len < index
    set(handles.text6,'String','Not existing!!!!!');
else
    set(handles.text6,'String','');
    
    handles.Current=index;
    
    [handles.tmp,handles.P]=ShowImg(hObject, eventdata, handles,1);
    
end
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Start=get(hObject,'String');
Start=str2double(Start);
handles.ImgIn=Start;
[handles.tmp,handles.P]=ShowImg(hObject, eventdata, handles,0);
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in Select_Rect.
function Select_Rect_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Rect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
handles.rect = getrect;

XL=handles.rect(1);XM=XL+handles.rect(3);
YL=handles.rect(2);YM=YL+handles.rect(4);
xv=[XL;XL;XM;XM;XL];
yv=[YL;YM;YM;YL;YL];
handles.polygon=[xv,yv];
% rect
rectangle('position',handles.rect,'Curvature',[0 0],'Edgecolor','m','Linestyle',':');
guidata(handles.figure1,handles);

% --- Executes on button press in Reset_Rect.
function Reset_Rect_Callback(hObject, eventdata, handles)
% hObject    handle to Reset_Rect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'rect')
    handles.rect = [];
    delete(findobj(handles.axes1,'type','rectangle'))
end
guidata(handles.figure1,handles);


% --- Executes on button press in Selcct_Polygon.
function Selcct_Polygon_Callback(hObject, eventdata, handles)
% hObject    handle to Selcct_Polygon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

but=1;
i=0;
L=length(handles.tmp);
set(handles.edit4,'String',num2str(L));

while but==1
    
    [x y but]=ginput(1);
    if but==1
        i=1+i;
        polygon(i,:)=[x,y];
        if i==1
            handles.ploy(i)=plot(x,y,'m*');
        else
            handles.ploy(i)=plot(polygon(i-1:i,1),polygon(i-1:i,2),'m--');
        end
    end
end

handles.polygon=[polygon;polygon(1,:)];

handles.ploy(i+1)=plot(handles.polygon(i:i+1,1),handles.polygon(i:i+1,2),'m--');
guidata(handles.figure1,handles);


% --- Executes on button press in Free_Hand_Draw.
function Free_Hand_Draw_Callback(hObject, eventdata, handles)
% hObject    handle to Free_Hand_Draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[myobj,xs,ys,zs] = freehanddraw(gca,'color','m','linewidth',2);
handles.polygon=[xs,ys];
handles.ploy=myobj;
guidata(handles.figure1,handles);

% --- Executes on button press in Reset_Polygon.
function Reset_Polygon_Callback(hObject, eventdata, handles)
% hObject    handle to Reset_Polygon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.polygon=[];

if isfield(handles,'ploy')
    delete(findobj(handles.axes1,'Color','m'))
end

if isfield(handles,'rect')
    handles.rect = [];
    delete(findobj(handles.axes1,'type','rectangle'))
end

guidata(handles.figure1,handles);

% --- Executes on button press in Batch_Keep.
function Batch_Keep_Callback(hObject, eventdata, handles)
% hObject    handle to Batch_Keep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

xv=handles.polygon(:,1);
yv=handles.polygon(:,2);

A=handles.tmp;
% I=find(XL<A(:,1) & A(:,1)<XM & YL< A(:,2) & A(:,2)<YM);

in=inpolygon(A(:,1),A(:,2),xv,yv);

I=find(in==0);

set(handles.P(I),'visible','off');

handles.tmp(I,:)=[];
handles.P(I)=[];
L=length(handles.tmp);
set(handles.edit4,'String',num2str(L));


handles.output = hObject;
%Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Batch_Delete.
function Batch_Delete_Callback(hObject, eventdata, handles)
% hObject    handle to Batch_Delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

xv=handles.polygon(:,1);
yv=handles.polygon(:,2);

A=handles.tmp;
% I=find(XL<A(:,1) & A(:,1)<XM & YL< A(:,2) & A(:,2)<YM);

in=inpolygon(A(:,1),A(:,2),xv,yv);

I=find(in==1);

set(handles.P(I),'visible','off');

handles.tmp(I,:)=[];
handles.P(I)=[];
L=length(handles.tmp);
set(handles.edit4,'String',num2str(L));


handles.output = hObject;
%Update handles structure
guidata(hObject, handles);



% --- Executes on button press in Load_Mat.
function Load_Mat_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Mat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile({'*.mat';},'Choose existing project');
if isequal(filename,0)
    disp('Cancled');
else
    str=[pathname filename];
    load(char(str));
    
    
    handles.Img=Pro.Img;
    if isfield(Pro,'c')
        handles.c=Pro.c;
        handles.OrignalC=Pro.c;
    else
        handles.c=Pro.Res;
    end
    if isfield(Pro,'Res')
        handles.Res=Pro.Res;
    end
    handles.OrigImg=Pro.OrigImg;
    if isfield(Pro,'Type')
        handles.Type=Pro.Type;
    else
        handles.Type=1;
    end
    handles.Flag=1;%Pro.Flag;

    handles.Size=length(handles.Img(1,1,:));
    handles.ImgIn=Pro.ImgIn; % Start slice
    if isfield(Pro,'TrackingWay')
        handles.TrackingWay=Pro.TrackingWay;
    else 
        handles.TrackingWay=1;
    end
    axes(handles.axes1);
    imshow(handles.Img(:,:,handles.ImgIn),[]);
    hold on;
    
    L=length(find(handles.c(:,1)~=0));
    set(handles.edit4,'String',num2str(L));
    set(handles.edit3,'String',length(handles.c(1,:))/2);
    set(handles.edit5,'String',1);
    
    
    centers=handles.c(:,1:2);
    
    for i=1:length(centers)
        
        P(i)=plot(centers(i,1),centers(i,2),'g+','MarkerSize',4);
        
    end
    handles.P=P;
    
    handles.tmp=centers;
    
    handles.Current=1;
end
% Choose default command line output for Tracking_Point
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uiputfile({'*.mat';},'Saving current project');
str=[pathname filename];

% Pro=handles;

Pro.Img=handles.Img;
Pro.c=handles.c;
Pro.Size=handles.Size;
Pro.OrigImg=handles.OrigImg;
Pro.Type=handles.Type;
Pro.Flag=handles.Flag;
Pro.Thred2=handles.Thred2;
Pro.ImgIn=handles.ImgIn;
Pro.TrackingWay=handles.TrackingWay;
Pro.Res=handles.Res;

%Pro.Res=handles.Res;
%Pro.YD=handles.YD;

save(char(str),'Pro');

% --- Executes on button press in Delete_Point.
function Delete_Point_Callback(hObject, eventdata, handles)
% hObject    handle to Delete_Point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
but=1;
P=handles.P;

while but==1
    [x y but]=ginput(1);
    if but==1
        set(handles.edit1,'String',num2str(x));
        set(handles.edit2,'String',num2str(y));
        %[L,I]=min(abs(handles.c(:,1)-x+handles.c(:,2)-y));
        [L,I]=min((handles.tmp(:,1)-x).^2+(handles.tmp(:,2)-y).^2);
        if L<=16
            set(handles.P(I),'visible','off');
            Mark=I;
            handles.tmp(Mark,:)=[];
            handles.P(Mark)=[];
            L=length(handles.tmp);
            set(handles.edit4,'String',num2str(L));
        end
    end
end

handles.output = hObject;

%Update handles structure
guidata(hObject, handles);

function [C P]=ShowImg(hObject,evendata,handles,check)
%Show Image!!!!
%hold off;
%when check==0 predict the new image
%when check==1 show the exsiting centers

TrackingWay=handles.TrackingWay;
Method=handles.Method;
hold off;
Ind=handles.ImgIn+TrackingWay*(handles.Current-1);
Index=handles.Current;
Type=handles.Type; % type==1 for tif (original images) type==2 for h5 (processed img);
I=handles.Img(:,:,Ind);
if check==0
    
    
    switch Method
        case 1
            im=zeros(size(I));
            im(I>handles.Thred)=I(I>handles.Thred);
            [centers,~,~] = imfindcircles(im,[2 4],'Method','TwoStage','Sensitivity',handles.Thred2);
            centers=round(centers);
        case 2
            [I,centers]=FindcenterSobel(I,handles.Thred,handles.Thred2);
        case 3
            if Type==2
                centers=h5findcenter(handles.PredictImg(:,:,Ind),handles.Thred/2000,Type);
            else
                centers=FindcenterBW(I,handles.Thred);
            end
    end
    
    
elseif check==1   
    centers=handles.c(:,2*Index-1:2*Index);   
end

axes(handles.axes1);
imshow(I,[]);
centers(all(centers==0,2),:)=[];
centers=unique(centers,'rows');
hold on;

Col=['g','r','y'];
M=1;
if isempty(handles.Res)
    LOR=0;
else
    LOR=length(handles.Res(1,:));
end
if isempty(handles.Res)==0 && handles.AnDone==0 && handles.Current*2 <= LOR
    
    COP=handles.Res(:,2*Index-1:2*Index);
    COld=handles.OriginalC(:,2*Index-1:2*Index);
    
    for i=1:length(centers)
        M=2;
        if ismember(centers(i,:),COP,'rows')
            M=1;
        end
        if ismember(centers(i,:),COld,'rows')==0
            M=3;
        end
        P(i)=plot(centers(i,1),centers(i,2),strcat(Col(M),'+'),'MarkerSize',4);
    end
    set(handles.text21,'String','green: correct points.') 
    set(handles.text23,'String','red: not used points.')
    set(handles.text22,'String','yellow: added points.');
else

    for i=1:length(centers)
        P(i)=plot(centers(i,1),centers(i,2),strcat(Col(M),'+'),'MarkerSize',4);
    end

    set(handles.text21,'String','');
    set(handles.text23,'String','');
    set(handles.text22,'String','');
end


C=centers;


L=length(C);
set(handles.edit4,'String',num2str(L));
[~,LengthC]=size(handles.c);
set(handles.edit3,'String',num2str(LengthC/2));
set(handles.edit5,'String',num2str(handles.Current));

%handles.c=cat(2,handles.c,round(centers));

% Choose default command line output for Tracking_Point
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in SliceBrowser.
function SliceBrowser_Callback(hObject, eventdata, handles)
% hObject    handle to SliceBrowser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.OrigImg)
    SliceBrowser(handles.Img);
else
    SliceBrowser(handles.OrigImg);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Thred2=get(hObject,'Value');
handles.Thred2=Thred2;
[handles.tmp handles.P]=ShowImg(hObject, eventdata, handles,0);

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in Delete_Plane.
function Delete_Plane_Callback(hObject, eventdata, handles)
% hObject    handle to Delete_Plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

L=size(handles.c,2);
if handles.Current == L/2
    set(handles.text6,'String','Delete The Last Image!!!!!');
    
    Current=handles.Current;
    handles.c=handles.c(:,1:L-2);
    handles.Current=handles.Current-1;
    
    [handles.tmp,handles.P]=ShowImg(hObject, eventdata, handles,1);
    set(handles.text6,'String','');
end

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in Method_Choice.
function Method_Choice_Callback(hObject, eventdata, handles)
% hObject    handle to Method_Choice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Method_Choice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Method_Choice


% --- Executes during object creation, after setting all properties.
function Method_Choice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Method_Choice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str=get(hObject,'String');
val=get(hObject,'Value');

handles.Method=val;
% switch str{val};
%     case 'ImfindCircles'
%         handles.Method=1;
%     case 'Sobel'
%         handles.Method=2;
%     case 'Watershed'
%         handles.Method=3;
% end

[handles.tmp,handles.P]=ShowImg(hObject, eventdata, handles,0);
guidata(hObject,handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel5.
function uipanel5_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel5 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
val=get(hObject,'Tag');
if strcmp(val,'Tracking_Bak')
    handles.TrackingWay=-1;
else
    handles.TrackingWay=1;
end

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function uipanel6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes when selected object is changed in uipanel6.
function uipanel6_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel6 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
val=get(hObject,'Tag');


if isempty(handles.Res)
    set(handles.Original,'Value',1);
else
    if strcmp(val,'Tracked')
        handles.AnDone=1;
        handles.c=handles.Res;
    else
        handles.AnDone=0;
        handles.c=handles.OriginalC;
    end
    [handles.tmp,handles.P]=ShowImg(hObject, eventdata, handles,1);
    
end

guidata(hObject,handles);
