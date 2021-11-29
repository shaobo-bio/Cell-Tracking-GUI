function varargout = GUI_Analysis(varargin)
% GUI_ANALYSIS MATLAB code for GUI_Analysis.fig
%      GUI_ANALYSIS, by itself, creates a new GUI_ANALYSIS or raises the existing
%      singleton*.
%
%      H = GUI_ANALYSIS returns the handle to a new GUI_ANALYSIS or the handle to
%      the existing singleton*.
%
%      GUI_ANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_ANALYSIS.M with the given input arguments.
%
%      GUI_ANALYSIS('Property','Value',...) creates a new GUI_ANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Analysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Analysis

% Last Modified by GUIDE v2.5 30-Apr-2014 16:52:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_Analysis_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_Analysis_OutputFcn, ...
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


% --- Executes just before GUI_Analysis is made visible.
function GUI_Analysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Analysis (see VARARGIN)

% Choose default command line output for GUI_Analysis
% Res saved the centers position!!
% YD saved the Y distance of cells in a pair!!


set(gcf,'menubar','figure');
set(gcf,'ToolBar','figure');
handles.in=[];
if nargin > 3
    handles.Line=1;
    YD=varargin{1};
    Res=varargin{2};
    Img=varargin{3};
    Start=varargin{4};
    OrigImg=varargin{5};
    handles.Flag=varargin{6};
    PFit=varargin{7};
    handles.TrackingWay=varargin{8};
    handles.Type=varargin{9};
    handles.c=varargin{10};
   
    handles.D=[];
    handles.P=[];
    
    L=length(Res(1,:));
    LD=[];DD=[];DP=[];
    axes(handles.axes1);
    imshow(Img(:,:,Start),[]);
    
    for j=1:1:L/2
        x=[];
        y=[];
        
        for i=1:length(Res)/2
            x=[x;(Res(2*i,2*j)+Res(2*i-1,2*j))/2];
            y=[y;(Res(2*i,2*j-1)+Res(2*i-1,2*j-1))/2];
        end
        
        P=PFit(j,:);
        
        %     for i=1:length(Res)
        %         Col=['r','g','b','y'];
        %
        %         A=Res(i,1:j*2);
        %         A=reshape(A,2,j);
        %         A=A';
        %         M=mod(round(i/2),4);
        %         if handles.Line==1
        %             plot(A(:,1),A(:,2),'-g','LineWidth',1);
        %         end
        %         plot(A(j,1),A(j,2),strcat(Col(M+1),'*'),'MarkerSize',3)
        %     end
        x=x';
        y=y';
        D=((P(1)*x)-y+P(2))/sqrt(P(1)^2+1);
        
        LD=[LD;sum(D.^2)];
        DD=[DD;D]; % DD is the distance of the centers of the pair cells to the fit central line
        
        DoP=abs(P(1)*Res(:,2*j)-Res(:,2*j-1)+P(2))/sqrt(P(1)^2+1);
        DP(j,:)=DoP;
        
    end
    
    % axes(handles.axes2);
    % plot(YD);
    %
    % axes(handles.axes3);
    % plot(LD,'r');
    [~,W]=size(Res);
    for i=2:W/2
        Speed(:,i-1)=sqrt((Res(:,i*2)-Res(:,i*2-2)).^2+(Res(:,2*i-1)-Res(:,2*i-3)).^2);
    end
    
    for i=1:length(DP(1,:))/2
        DP2(:,i)=DP(:,2*i)+DP(:,2*i-1);
    end
    handles.DP=DP2;
    handles.Speed=Speed;
    handles.P=PFit;
    handles.Res3D=[0,0,0];
    handles.L=L/2;
    handles.Start=Start;
    handles.D=DD;
    handles.YD=YD;
    handles.LD=LD;
    handles.Res=Res;
    handles.Img=Img;
    handles.OrigImg=OrigImg;
end
handles.output = hObject;
% LD saved the distance of the central point of cell pairs to the fitted
% central line!
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Analysis_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Show Next Image  Detail!!!!!!!!!!


if handles.Index< length(handles.Res(1,:))/2
    set(handles.text5,'String','');
    handles.Index=handles.Index+1;
    handles.Pi=Show_Detail(hObject, eventdata, handles,handles.Index);
else
    set(handles.text5,'String','This is the last Image!!','ForegroundColor','red');
end


handles.output = hObject;

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Start Show Image Detail!!!!!!!!!!

handles.Pi=Show_Detail(hObject, eventdata, handles,1);

handles.Index=1;
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

function Pi=Show_Detail(hObject, eventdata, handles,Index)
% Show Image one by one


YD=handles.YD;
Res=handles.Res;
Img=handles.Img;
D=handles.D;
DP=handles.DP;
Way=handles.TrackingWay;
InJ=handles.Start+Way*(Index-1);
in=handles.in;

axes(handles.axes1);
imshow(Img(:,:,InJ),[]);
hold on;

P=handles.P(Index,:);
x1=linspace(0,500);
y1=polyval(P,x1);
%plot(y,x,'--w');
[L,W]=size(Res);
PT(L+1)=plot(y1,x1,'--m');

% for i=1:L
% for i=L:-1:1
%     Col=['r','g','b','y'];
%     
%     A=Res(i,1:Index*2);
%     A=reshape(A,2,Index);
%     A=A';
%     M=mod(round(i/2),4);
%     PT(i)=plot(A(:,1),A(:,2),'-g','LineWidth',1);
%     
%     plot(A(Index,1),A(Index,2),strcat(Col(M+1),'*'),'MarkerSize',3);
%     
% end

Count=1;
for i=1:L
    Col=['r','g','b','y'];
    
    A=Res(i,1:Index*2);
    A=reshape(A,2,Index);
    A=A';
    if i > 2 && mod(i,2)
        if ~isequal(Res(i,1:2),Res(i-2,1:2))
            if ~isequal(Res(i+1,1:2),Res(i-1,1:2))
                Count=Count+1;
            end
        end
    end
    
    M=mod(Count,4);
    PT(i)=plot(A(:,1),A(:,2),'-g','LineWidth',1);
    
    plot(A(Index,1),A(Index,2),strcat(Col(M+1),'*'),'MarkerSize',3); 
    
    
end

if ~isempty(in)
        
        Point=Res(:,Index*2-1:Index*2);
        plot(Point(in,1),Point(in,2),'co','MarkerSize',9);
end

hold off;

set(handles.edit2,'String',num2str(Index));



axes(handles.axes2);
plot(YD(Index,:),'rh');
Ymax=max(YD(:));Ymin=min(YD(:));
axis([0 length(YD(1,:)) Ymin Ymax]);
hold on;
if ~isempty(in)
    ind=unique(round(find(in==1)/2));
    Point=YD(Index,ind);
    plot(ind,Point,'ko');
end
hold off;


axes(handles.axes3);
plot(D(Index,:),'gh');
Ymax=max(D(:));Ymin=min(D(:));
axis([0 length(D(1,:)) Ymin Ymax]);
hold on;
if ~isempty(in)
    Point=D(Index,ind);
    plot(ind,Point,'ko');
end
hold off;

axes(handles.axes4);

Shape=Res(:,Index*2-1:Index*2);
XX=Res(:,1:2:W);
Xmin=min(XX(:))-5;
Xmax=max(XX(:))+5;
Shape1=Shape(1:2:L-1,:);
Shape2=Shape(2:2:L,:);
Shape1=sortrows(Shape1,2);
Shape2=sortrows(Shape2,2);

plot(Shape1(:,1),512-Shape1(:,2),'-g','LineWidth',2);
set(gca,'color',[0,0,0]);
axis([Xmin Xmax 0 512]);
% axis off;
hold on;
plot(Shape2(:,1),512-Shape2(:,2),'-g','LineWidth',2);
if ~isempty(in);
    Point=Shape(find(in==1),:);
    plot(Point(:,1),512-Point(:,2),'yo');
end
hold off;

axes(handles.axes5);
cla reset;
SpeedMax=max(handles.Speed(:));
if Index >= 2
%     Pre=Res(:,Index*2-3:Index*2-2);
%     Speed=sqrt((Shape(:,1)-Pre(:,1)).^2+(Shape(:,2)-Pre(:,2)).^2);
    Speed=handles.Speed(:,Index-1);
    
    plot(Speed,'bh');
    axis([0 L 0 SpeedMax]);
    if ~isempty(in)
        hold on;
        InP=find(in==1);
        Point=Speed(InP);
        plot(InP,Point,'ro');
        hold off;
    end
end




axes(handles.axes6);
cla;
L=length(DP(1,:));
Y=max(DP(Index,:))+3;
plot(DP(Index,:),'b');
axis([0 L 0 Y]);
if ~isempty(in)
    hold on;   
    Point=handles.DP(Index,:);
    plot(ind,Point(ind),'ro');
    hold off;
end
Pi=PT;


%%

% figure;
% subplot(2,2,[1,3]);
% 
% imshow(Img(:,:,InJ),[]);
% hold on;
% 
% P=handles.P(Index,:);
% x1=linspace(0,500);
% y1=polyval(P,x1);
% %plot(y,x,'--w');
% [L,W]=size(Res);
% PT(L+1)=plot(y1,x1,'--m');
% 
% % for i=1:L
% % for i=L:-1:1
% %     Col=['r','g','b','y'];
% %     
% %     A=Res(i,1:Index*2);
% %     A=reshape(A,2,Index);
% %     A=A';
% %     M=mod(round(i/2),4);
% %     PT(i)=plot(A(:,1),A(:,2),'-g','LineWidth',1);
% %     
% %     plot(A(Index,1),A(Index,2),strcat(Col(M+1),'*'),'MarkerSize',3);
% %     
% % end
% 
% Count=1;
% for i=1:L
%     Col=['r','g','b','y'];
%     
%     A=Res(i,1:Index*2);
%     A=reshape(A,2,Index);
%     A=A';
%     if i > 2 && mod(i,2)
%         if ~isequal(Res(i,1:2),Res(i-2,1:2))
%             if ~isequal(Res(i+1,1:2),Res(i-1,1:2))
%                 Count=Count+1;
%             end
%         end
%     end
%     
%     M=mod(Count,4);
%     PT(i)=plot(A(:,1),A(:,2),'-g','LineWidth',1);
%     
%     plot(A(Index,1),A(Index,2),strcat(Col(M+1),'*'),'MarkerSize',3); 
%     
%     
% end
% 
% if ~isempty(in)
%         
%         Point=Res(:,Index*2-1:Index*2);
%         plot(Point(in,1),Point(in,2),'co','MarkerSize',9);
% end
% 
% hold off;
% 
% subplot(2,2,4)
% L=length(DP(1,:));
% Y=max(DP(Index,:))+3;
% plot(DP(Index,:),'b');
% axis([0 L 0 Y]);
% if ~isempty(in)
%     hold on;   
%     Point=handles.DP(Index,:);
%     plot(ind,Point(ind),'ro');
%     hold off;
% end
% title('The outline of the heart cells');
% 
% subplot(2,2,2)
% Shape=Res(:,Index*2-1:Index*2);
% XX=Res(:,1:2:W);
% Xmin=min(XX(:))-5;
% Xmax=max(XX(:))+5;
% Shape1=Shape(1:2:L-1,:);
% Shape2=Shape(2:2:L,:);
% Shape1=sortrows(Shape1,2);
% Shape2=sortrows(Shape2,2);
% 
% plot(Shape1(:,1),512-Shape1(:,2),'-g','LineWidth',2);
% set(gca,'color',[0,0,0]);
% axis([Xmin Xmax 0 512]);
% % axis off;
% hold on;
% plot(Shape2(:,1),512-Shape2(:,2),'-g','LineWidth',2);
% if ~isempty(in);
%     Point=Shape(find(in==1),:);
%     plot(Point(:,1),512-Point(:,2),'yo');
% end
% hold off;
% title('The Distance of paired cells');

    




% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Show details of total stastic information!!!!

axes(handles.axes2);
plot(handles.YD);

axes(handles.axes3);
plot(handles.LD,'r');

axes(handles.axes5);
cla;
if handles.Index >= 2
%     Pre=Res(:,Index*2-3:Index*2-2);
%     Speed=sqrt((Shape(:,1)-Pre(:,1)).^2+(Shape(:,2)-Pre(:,2)).^2);
    Speed=handles.Speed;
    plot(Speed');
   
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Show details of individual stastic information!!!!

axes(handles.axes2);
plot(handles.YD(handles.Index,:),'rh');
axis([0 length(handles.YD(1,:)) -20 20]);

axes(handles.axes3);
plot(handles.D(handles.Index,:),'gh');
axis([0 length(handles.D(1,:)) -10 10]);
axes(handles.axes5);
cla;
if handles.Index >= 2
%     Pre=Res(:,Index*2-3:Index*2-2);
%     Speed=sqrt((Shape(:,1)-Pre(:,1)).^2+(Shape(:,2)-Pre(:,2)).^2);
    Speed=handles.Speed(:,handles.Index-1);
    bar(Speed,'g');
    axis([0 length(handles.Speed) 0 15]);
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Show the previous image!!!

if handles.Index > 1
    set(handles.text5,'String','');
    handles.Index=handles.Index-1;
    handles.Pi=Show_Detail(hObject, eventdata, handles,handles.Index);
else
    set(handles.text5,'String','This is the first Image!!','ForegroundColor','red');
end

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=get(hObject,'String');
index=str2double(index);

if length(handles.Res(1,:))/2 < index
    set(handles.text5,'String','Not existing!!!!!');
else
    set(handles.text5,'String','');
    
    handles.Index=index;
    
    handles.Pi=Show_Detail(hObject, eventdata, handles,handles.Index);
end
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
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


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Save the analysis result!!!!

%[filename,pathname]=uiputfile({'';},'Saving current project');
%str=[pathname filename];
xlswrite('.\Result\Centers.xls',handles.Res);
xlswrite('.\Result\Centers3D.xls',handles.Res3D);
xlswrite('.\Result\Distance_of_Cell_pairs.xls',handles.YD);
xlswrite('.\Result\Distance2Central_Line.xls',handles.LD);



% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Lines Off!!!
set(handles.Pi,'visible','off');



% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Lines ON!!!
set(handles.Pi,'visible','on');


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 3D modelling

Way=handles.TrackingWay;
 handles.OrigImg=[];

if any(handles.Res3D) == 0

    if handles.Flag==1 && isempty(handles.OrigImg)
        
        [filename,pathname]=uigetfile({'*.tif'},'Choose orignal images');
        
        if isequal(filename,0)
            disp('Cancled');
        else
            
            str=[pathname filename];
            ImgO=ULoadTiff_Z(str);
            
        end
    else
        ImgO=handles.OrigImg;
        
    end
    
    [W,L,SliceZ,~]=size(ImgO);
    Z_Res=[];
    
    for i=1:handles.L
        Index=handles.Start+Way*(i-1);
        for j=1:length(handles.Res)
            
            %         IndS=SliceZ*(handles.Start+i-2)+1;
            %         IndE=SliceZ*(handles.Start+i-1);
%             IndS=SliceZ*(handles.L-(handles.Start+i-1))+1;
%             IndE=SliceZ*(handles.L-(handles.Start+i-2));
% 
%             IndS=SliceZ*(Index-1)+1;
%             IndE=SliceZ*Index;
            X=handles.Res(j,2*i-1);
            Y=handles.Res(j,2*i);
           
            if(X+2 < W && Y+2<L && X-2 >0 && Y-2>0)
                ProbeX=linspace(X-2,X+2,3);
                ProbeY=linspace(Y-2,Y+2,3);
                [D Z]=max(ImgO(ProbeY,ProbeX,:,Index),[],3);
            elseif Y>L
 
                [D Z]=max(ImgO(512,X,:,Index),[],3);
            else
                [D Z]=max(ImgO(Y,X,:,Index),[],3);
            end
            
            ZL(j)=mean(Z(:));
        end
        TmpRes(:,3*i-2:3*i-1)=handles.Res(:,2*i-1:2*i);
        TmpRes(:,3*i)=ZL';
        
    end
    
    handles.Res3D=TmpRes;
end

% axes(handles.axes1);
T3Dmodel_Z(handles.Res3D);

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in Save_Project.
function Save_Project_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Project (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uiputfile({'*.mat';},'Saving current project');
str=[pathname filename];

Pro.Res3D=handles.Res3D;
Pro.L=handles.L;

Pro.ImgIn=handles.Start;
Pro.D=handles.D;
Pro.YD=handles.YD;
Pro.LD=handles.LD;
Pro.Res=handles.Res;
Pro.Img=handles.Img;
% Pro.OrigImg=handles.OrigImg;
Pro.P=handles.P;
Pro.Line=handles.Line;
Pro.TrackingWay=handles.TrackingWay;
Pro.Speed=handles.Speed;
Pro.Flag=handles.Flag;
Pro.DP=handles.DP;
if isfield(handles,'c')
    Pro.c=handles.c;
end
if isfield(handles,'Type')
    Pro.Type=handles.Type;
else
    Pro.Type=1;
end

save(char(str),'Pro');


% --- Executes on button press in Load_Project.
function Load_Project_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Project (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile({'*.mat';},'Choose existing project');
if isequal(filename,0)
    disp('Cancled');
else
    
    str=[pathname filename];
    load(char(str));
end

handles.Res3D=Pro.Res3D;
handles.L=Pro.L;
if isfield(Pro,'Start')
    handles.Start=Pro.Start;
else
    handles.Start=Pro.ImgIn;
end
handles.D=Pro.D;
handles.YD=Pro.YD;
handles.LD=Pro.LD;
handles.Res=Pro.Res;
handles.Img=Pro.Img;
% handles.OrigImg=Pro.OrigImg;
handles.P=Pro.P;
handles.Line=Pro.Line;
handles.TrackingWay=Pro.TrackingWay;
handles.Flag=Pro.Flag;

% % Res=handles.Res;
% % L=length(handles.Res(1,:));
% %     LD=[];DD=[];DP=[];D2Central=[];
% %    for j=1:1:L/2
% %         x=[];
% %         y=[];
% %         
% %         for i=1:length(handles.Res(:,1))/2
% %             x=[x;(Res(2*i,2*j)+Res(2*i-1,2*j))/2];
% %             y=[y;(Res(2*i,2*j-1)+Res(2*i-1,2*j-1))/2];
% %         end
% %         
% %         P=handles.P(j,:);
% %         x=x';
% %         y=y';
% %         D=((P(1)*x)-y+P(2))/sqrt(P(1)^2+1);
% %         
% %         LD=[LD;sum(D.^2)];
% %         DD=[DD;D]; % DD is the distance of the centers of the pair cells to the fit central line
% %         
% %         DoP=abs(P(1)*Res(:,2*j)-Res(:,2*j-1)+P(2))/sqrt(P(1)^2+1);
% %         DP(j,:)=DoP;
% %         
% %         for i=1:length(handles.Res(:,1))/2
% %             a=Res(2*i,2*j)/P(1)+Res(2*i,2*j-1);
% %             b=Res(2*i-1,2*j)/P(1)+Res(2*i-1,2*j-1);
% %             
% %             D=(a-b)*P(1)/sqrt(P(1)^2+1);
% %             %         D=abs(Res(2*i,2*j)-Res(2*i-1,2*j));
% %             D2Central(j,i)=D;
% %         end
% %         
% %     end
% %     
% %     % axes(handles.axes2);
% %     % plot(YD);
% %     %
% %     % axes(handles.axes3);
% %     % plot(LD,'r');
% %     handles.YD=D2Central;
% %     [~,W]=size(Res);
% %     for i=2:W/2
% %         Speed(:,i-1)=sqrt((Res(:,i*2)-Res(:,i*2-2)).^2+(Res(:,2*i-1)-Res(:,2*i-3)).^2);
% %     end
% %     
% %     for i=1:length(DP(1,:))/2
% %         DP2(:,i)=DP(:,2*i)+DP(:,2*i-1);
% %     end
% %     handles.DP=DP2;
% %     handles.Speed=Speed;


handles.DP=Pro.DP;
if isfield(Pro,'Speed')
    handles.Speed=Pro.Speed;
else
    handles.Speed=[];
end

axes(handles.axes1);
imshow(handles.Img(:,:,handles.Start),[]);


handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Show_Process.
function Show_Process_Callback(hObject, eventdata, handles)
% hObject    handle to Show_Process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Res=handles.Res;
axes(handles.axes1);
Way=handles.TrackingWay;
[L W]=size(Res);
for j=1:1:handles.L
%     hold off;
%     set(handles.edit2,'String',num2str(j));
%     InJ=handles.Start+Way*(j-1);
%     imshow(handles.Img(:,:,InJ),[]);
%     hold on;
%     
%     x1=linspace(0,500);
%     y1=polyval(handles.P(j,:),x1);
%     %plot(y,x,'--w');
%     plot(y1,x1,'--m');
%     
% %     for i=length(Res):-1:1
% for i=1:length(Res)
%         Col=['r','g','b','y'];
%         
%         A=Res(i,1:j*2);
%         A=reshape(A,2,j);
%         A=A';
%         M=mod(round(i/2),4);
%         if handles.Line==1
%             plot(A(:,1),A(:,2),'-g','LineWidth',1);
%         end
%         plot(A(j,1),A(j,2),strcat(Col(M+1),'*'),'MarkerSize',3)
%     end
%     
%     
%     

P=Show_Detail(hObject, eventdata, handles,j);
pause(0.6);
end

handles.Index=j;
handles.Pi=P;

axes(handles.axes2);
plot(handles.YD);

axes(handles.axes3);
plot(handles.LD,'r');

axes(handles.axes5);
if isempty(handles.Speed)==1

for i=2:W/2
    Speed(:,i-1)=sqrt((Res(:,i*2)-Res(:,i*2-2)).^2+(Res(:,2*i-1)-Res(:,2*i-3)).^2);
end

handles.Speed=Speed;
else
    Speed=handles.Speed;
end
plot(Speed');
guidata(hObject, handles);


% --- Executes on button press in Select_Region.
function Select_Region_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'ploy')
    delete(findobj(handles.axes1,'Color','c'));
%     delete(findobj(handles.axes1,'Color','c'));
    delete(findobj(handles.axes4,'Color','y'));
    delete(findobj(handles.axes6,'Color','r'));
end

[myobj,xs,ys,zs] = freehanddraw(gca,'color','c','linewidth',2);
in=[];
if zs ==1 
    A=handles.Res(:,handles.Index*2-1:handles.Index*2);
    in=inpolygon(A(:,1),A(:,2),xs,ys);
    
    Speed=handles.Speed;
    Speed(find(in==0),:)=[];
    axes(handles.axes5);
    cla reset;
    plot(Speed');
    handles.ploy=myobj;
    guidata(handles.figure1,handles);
    
    index=unique(round(find(in==1)/2));
    YD=handles.YD(:,index);
    axes(handles.axes2);
    cla reset;
    plot(YD);
    
    D=zeros(length(handles.D),1);
    D(index)=abs(handles.D(handles.Index,index));
    
    axes(handles.axes3);
    cla reset;
    bar(D,'c');
    
    axes(handles.axes4);
    hold on;
    Point=A(find(in==1),:);
    plot(Point(:,1),512-Point(:,2),'yo');
    hold off;
    
    axes(handles.axes6);
    hold on;
    Point=handles.DP(handles.Index,:);
    plot(index,Point(index),'ro');
    hold off;
end
handles.in=in;
guidata(hObject, handles);
