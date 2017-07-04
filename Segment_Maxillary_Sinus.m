
%Computed Tomography-Based Volumetric Tool for Standardized Measurement of the Maxillary Sinus 

% Read images
clear all
close all

[filename, pathname] = uigetfile('*','Select CT exam (all slices)','MultiSelect','on');
num= length(filename);
bbbb=1;


step=input('Step of image reading: \')


for aaaa = 1:step:num
    xinfo=dicominfo([pathname,char(filename(aaaa))]);
    pxsp=cat(2,xinfo.PixelSpacing);
    x=dicomread([pathname,char(filename(aaaa))])+cat(2,xinfo.RescaleIntercept);
    k=x;
    k = im2bw(k,0.49);
    k = imfill(k,'holes');
    cc = bwconncomp(k);
    stats = regionprops(cc,'Area');
    A = [stats.Area];
    [~,biggest] = max(A);
    k(labelmatrix(cc)~=biggest) = 0;
    x(k~=1)=-2000;     
    masccranio(:,:,bbbb)=k;
    cranio(:,:,bbbb)=x;
    cranio_full(:,:,bbbb)=x;
    bbbb=bbbb+1;
    end

%% middle slice of the CT exam

continua=2;
ff=round((size(cranio,3)/3));

while (continua==2 | continua==3)
x=cranio(:,:,ff);
xlim=zeros(512);
for i = 1:size(x,1)
    for j = 1:size(x,2)
        if x(i,j)>=150; 
            xlim(i,j)=1;
        else xlim(i,j)=0;
        end
    end
end

%Morphological operators
BW2=bwmorph(xlim,'bridge');
BW3 =  bwmorph(BW2,'majority');

%Watershed
w=watershed(BW3);

%Rule-based system
w(450:end,:)=1;
figure,imshow(w,[])  
w2=zeros(512);
final=zeros(512);
w2(w~=0)=1;

a1=w2;
cc = bwconncomp(w2);
stats = regionprops(cc,'Area'); 
A = [stats.Area];
[~,biggest] = max(A);
a1(labelmatrix(cc)~=biggest) = 0;
w2(a1==1)=0; 

% Major areas 
for i=1:9
    seg=zeros(512);
    A(biggest)=0;
    [~,biggest] = max(A);
    seg(labelmatrix(cc)== biggest) = 1;
    segment(:,:,i)=seg;
    
end

%Correlation
segment_ok=[];
segment_inv=[];
a=1;
for i=1:9
    if (sum(sum(segment(:,:,i)))>0)
        segment_ok(:,:,a)=segment(:,:,i);
        segment_inv(:,:,a)= flipdim(segment(:,:,i) ,2);
        a=a+1;
    end
end

i=1;
for aa=1:size(segment_ok,3)-1
    coef(i)=corr2(segment_ok(:,:,1),segment_inv(:,:,aa+1));
    i=i+1;
end

if size(segment_ok,3)>2
    for aa=1:size(segment_ok,3)-2
        coef(i)=corr2(segment_ok(:,:,2),segment_inv(:,:,aa+2));
        i=i+1;
    end
end

if size(segment_ok,3)>3
    for aa=1:size(segment_ok,3)-3
        coef(i)=corr2(segment_ok(:,:,3),segment_inv(:,:,aa+3));
        i=i+1;
    end
end

if size(segment_ok,3)>4
    for aa=1:size(segment_ok,3)-4
        coef(i)=corr2(segment_ok(:,:,4),segment_inv(:,:,aa+4));
        i=i+1;
    end
end

if size(segment_ok,3)>5
    for aa=1:size(segment_ok,3)-5
        coef(i)=corr2(segment_ok(:,:,5),segment_inv(:,:,aa+5));
        i=i+1;
    end
end

if size(segment_ok,3)>6
    for aa=1:size(segment_ok,3)-6
        coef(i)=corr2(segment_ok(:,:,6),segment_inv(:,:,aa+6));
        i=i+1;
    end
end

if size(segment_ok,3)>7
    for aa=1:size(segment_ok,3)-7
        coef(i)=corr2(segment_ok(:,:,7),segment_inv(:,:,aa+7));
        i=i+1;
    end
end

if size(segment_ok,3)>8
    for aa=1:size(segment_ok,3)-8
        coef(i)=corr2(segment_ok(:,:,8),segment_inv(:,:,aa+8));
        i=i+1;
    end
end
[a,b]=max(coef);

%% Selecting only the MS areas 
    if size(segment_ok,3)==2
        sm1=segment_ok(:,:,1);
        sm2=segment_ok(:,:,2);
    end
    
    if size(segment_ok,3)==3
        if b==1
            sm1=segment_ok(:,:,1);
            sm2=segment_ok(:,:,2);
        else if b==2
                sm1=segment_ok(:,:,1);
                sm2=segment_ok(:,:,3);
            else if b==3
                    sm1=segment_ok(:,:,2);
                    sm2=segment_ok(:,:,3);
                end
            end
        end
    end
    
    if size(segment_ok,3)==4
        if b==1
            sm1=segment_ok(:,:,1);
            sm2=segment_ok(:,:,2);
        else if b==2
                sm1=segment_ok(:,:,1);
                sm2=segment_ok(:,:,3);
            else if b==3
                    sm1=segment_ok(:,:,1);
                    sm2=segment_ok(:,:,4);
                else if b==4
                        sm1=segment_ok(:,:,2);
                        sm2=segment_ok(:,:,3);
                    else if b==5
                            sm1=segment_ok(:,:,2);
                            sm2=segment_ok(:,:,4);
                        else if b==6
                                sm1=segment_ok(:,:,3);
                                sm2=segment_ok(:,:,4);
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    if size(segment_ok,3)==5
        if b==1
            sm1=segment_ok(:,:,1);
            sm2=segment_ok(:,:,2);
        else if b==2
                sm1=segment_ok(:,:,1);
                sm2=segment_ok(:,:,3);
            else if b==3
                    sm1=segment_ok(:,:,1);
                    sm2=segment_ok(:,:,4);
                else if b==4
                        sm1=segment_ok(:,:,1);
                        sm2=segment_ok(:,:,5);
                    else if b==5
                            sm1=segment_ok(:,:,2);
                            sm2=segment_ok(:,:,3);
                        else if b==6
                                sm1=segment_ok(:,:,2);
                                sm2=segment_ok(:,:,4);
                            else if b==7
                                    sm1=segment_ok(:,:,2);
                                    sm2=segment_ok(:,:,5);
                                else if b==8
                                        sm1=segment_ok(:,:,3);
                                        sm2=segment_ok(:,:,4);
                                    else if b==9
                                            sm1=segment_ok(:,:,3);
                                            sm2=segment_ok(:,:,5);
                                        else if b==10
                                                sm1=segment_ok(:,:,4);
                                                sm2=segment_ok(:,:,5);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    if size(segment_ok,3)==6
        if b==1
            sm1=segment_ok(:,:,1);
            sm2=segment_ok(:,:,2);
        else if b==2
                sm1=segment_ok(:,:,1);
                sm2=segment_ok(:,:,3);
            else if b==3
                    sm1=segment_ok(:,:,1);
                    sm2=segment_ok(:,:,4);
                else if b==4
                        sm1=segment_ok(:,:,1);
                        sm2=segment_ok(:,:,5);
                    else if b==5
                            sm1=segment_ok(:,:,1);
                            sm2=segment_ok(:,:,6);
                        else if b==6
                                sm1=segment_ok(:,:,2);
                                sm2=segment_ok(:,:,3);
                            else if b==7
                                    sm1=segment_ok(:,:,2);
                                    sm2=segment_ok(:,:,4);
                                else if b==8
                                        sm1=segment_ok(:,:,2);
                                        sm2=segment_ok(:,:,5);
                                    else if b==9
                                            sm1=segment_ok(:,:,2);
                                            sm2=segment_ok(:,:,6);
                                        else if b==10
                                                sm1=segment_ok(:,:,3);
                                                sm2=segment_ok(:,:,4);
                                            else if b==11
                                                    sm1=segment_ok(:,:,3);
                                                    sm2=segment_ok(:,:,5);
                                                else if b==12
                                                        sm1=segment_ok(:,:,3);
                                                        sm2=segment_ok(:,:,6);
                                                    else if b==13
                                                            sm1=segment_ok(:,:,4);
                                                            sm2=segment_ok(:,:,5);
                                                        else if b==14
                                                                sm1=segment_ok(:,:,4);
                                                                sm2=segment_ok(:,:,6);
                                                            else if b==15
                                                                    sm1=segment_ok(:,:,5);
                                                                    sm2=segment_ok(:,:,6);
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    if size(segment_ok,3)==7
        if b==1
            sm1=segment_ok(:,:,1);
            sm2=segment_ok(:,:,2);
        else if b==2
                sm1=segment_ok(:,:,1);
                sm2=segment_ok(:,:,3);
            else if b==3
                    sm1=segment_ok(:,:,1);
                    sm2=segment_ok(:,:,4);
                else if b==4
                        sm1=segment_ok(:,:,1);
                        sm2=segment_ok(:,:,5);
                    else if b==5
                            sm1=segment_ok(:,:,1);
                            sm2=segment_ok(:,:,6);
                        else if b==6
                                sm1=segment_ok(:,:,1);
                                sm2=segment_ok(:,:,7);
                            else if b==7
                                    sm1=segment_ok(:,:,2);
                                    sm2=segment_ok(:,:,3);
                                else if b==8
                                        sm1=segment_ok(:,:,2);
                                        sm2=segment_ok(:,:,4);
                                    else if b==9
                                            sm1=segment_ok(:,:,2);
                                            sm2=segment_ok(:,:,5);
                                        else if b==10
                                                sm1=segment_ok(:,:,2);
                                                sm2=segment_ok(:,:,6);
                                            else if b==11
                                                    sm1=segment_ok(:,:,2);
                                                    sm2=segment_ok(:,:,7);
                                                else if b==12
                                                        sm1=segment_ok(:,:,3);
                                                        sm2=segment_ok(:,:,4);
                                                    else if b==13
                                                            sm1=segment_ok(:,:,3);
                                                            sm2=segment_ok(:,:,5);
                                                        else if b==14
                                                                sm1=segment_ok(:,:,3);
                                                                sm2=segment_ok(:,:,6);
                                                            else if b==15
                                                                    sm1=segment_ok(:,:,3);
                                                                    sm2=segment_ok(:,:,7);
                                                                else if b==16
                                                                        sm1=segment_ok(:,:,4);
                                                                        sm2=segment_ok(:,:,5);
                                                                    else if b==17
                                                                            sm1=segment_ok(:,:,4);
                                                                            sm2=segment_ok(:,:,6);
                                                                        else if b==18
                                                                                sm1=segment_ok(:,:,4);
                                                                                sm2=segment_ok(:,:,7);
                                                                            else if b==19
                                                                                    sm1=segment_ok(:,:,5);
                                                                                    sm2=segment_ok(:,:,6);
                                                                                else if b==20
                                                                                        sm1=segment_ok(:,:,5);
                                                                                        sm2=segment_ok(:,:,7);
                                                                                    else if b==21
                                                                                            sm1=segment_ok(:,:,6);
                                                                                            sm2=segment_ok(:,:,7);
                                                                                        end
                                                                                    end
                                                                                end
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    if size(segment_ok,3)==8
        if b==1
            sm1=segment_ok(:,:,1);
            sm2=segment_ok(:,:,2);
        else if b==2
                sm1=segment_ok(:,:,1);
                sm2=segment_ok(:,:,3);
            else if b==3
                    sm1=segment_ok(:,:,1);
                    sm2=segment_ok(:,:,4);
                else if b==4
                        sm1=segment_ok(:,:,1);
                        sm2=segment_ok(:,:,5);
                    else if b==5
                            sm1=segment_ok(:,:,1);
                            sm2=segment_ok(:,:,6);
                        else if b==6
                                sm1=segment_ok(:,:,1);
                                sm2=segment_ok(:,:,7);
                            else if b==7
                                    sm1=segment_ok(:,:,1);
                                    sm2=segment_ok(:,:,8);
                                else if b==8
                                        sm1=segment_ok(:,:,2);
                                        sm2=segment_ok(:,:,3);
                                    else if b==9
                                            sm1=segment_ok(:,:,2);
                                            sm2=segment_ok(:,:,4);
                                        else if b==10
                                                sm1=segment_ok(:,:,2);
                                                sm2=segment_ok(:,:,5);
                                            else if b==11
                                                    sm1=segment_ok(:,:,2);
                                                    sm2=segment_ok(:,:,6);
                                                else if b==12
                                                        sm1=segment_ok(:,:,2);
                                                        sm2=segment_ok(:,:,7);
                                                    else if b==13
                                                            sm1=segment_ok(:,:,2);
                                                            sm2=segment_ok(:,:,8);
                                                        else if b==14
                                                                sm1=segment_ok(:,:,3);
                                                                sm2=segment_ok(:,:,4);
                                                            else if b==15
                                                                    sm1=segment_ok(:,:,3);
                                                                    sm2=segment_ok(:,:,5);
                                                                else if b==16
                                                                        sm1=segment_ok(:,:,3);
                                                                        sm2=segment_ok(:,:,6);
                                                                    else if b==17
                                                                            sm1=segment_ok(:,:,3);
                                                                            sm2=segment_ok(:,:,7);
                                                                        else if b==18
                                                                                sm1=segment_ok(:,:,3);
                                                                                sm2=segment_ok(:,:,8);
                                                                            else if b==19
                                                                                    sm1=segment_ok(:,:,4);
                                                                                    sm2=segment_ok(:,:,5);
                                                                                else if b==20
                                                                                        sm1=segment_ok(:,:,4);
                                                                                        sm2=segment_ok(:,:,6);
                                                                                    else if b==21
                                                                                            sm1=segment_ok(:,:,4);
                                                                                            sm2=segment_ok(:,:,7);
                                                                                        else if b==22
                                                                                                sm1=segment_ok(:,:,4);
                                                                                                sm2=segment_ok(:,:,8);
                                                                                            else if b==23
                                                                                                    sm1=segment_ok(:,:,5);
                                                                                                    sm2=segment_ok(:,:,6);
                                                                                                else if b==24
                                                                                                        sm1=segment_ok(:,:,5);
                                                                                                        sm2=segment_ok(:,:,7);
                                                                                                    else if b==25
                                                                                                            sm1=segment_ok(:,:,5);
                                                                                                            sm2=segment_ok(:,:,8);
                                                                                                        else if b==26
                                                                                                                sm1=segment_ok(:,:,6);
                                                                                                                sm2=segment_ok(:,:,7);
                                                                                                            else if b==27
                                                                                                                    sm1=segment_ok(:,:,6);
                                                                                                                    sm2=segment_ok(:,:,8);
                                                                                                                else if b==28
                                                                                                                        sm1=segment_ok(:,:,7);
                                                                                                                        sm2=segment_ok(:,:,8);
                                                                                                                    end
                                                                                                                end
                                                                                                            end
                                                                                                        end
                                                                                                    end
                                                                                                end
                                                                                            end
                                                                                        end
                                                                                    end
                                                                                end
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    if size(segment_ok,3)==9
        if b==1
            sm1=segment_ok(:,:,1);
            sm2=segment_ok(:,:,2);
        else if b==2
                sm1=segment_ok(:,:,1);
                sm2=segment_ok(:,:,3);
            else if b==3
                    sm1=segment_ok(:,:,1);
                    sm2=segment_ok(:,:,4);
                else if b==4
                        sm1=segment_ok(:,:,1);
                        sm2=segment_ok(:,:,5);
                    else if b==5
                            sm1=segment_ok(:,:,1);
                            sm2=segment_ok(:,:,6);
                        else if b==6
                                sm1=segment_ok(:,:,1);
                                sm2=segment_ok(:,:,7);
                            else if b==7
                                    sm1=segment_ok(:,:,1);
                                    sm2=segment_ok(:,:,8);
                                else if b==8
                                        sm1=segment_ok(:,:,1);
                                        sm2=segment_ok(:,:,9);
                                    else if b==9
                                            sm1=segment_ok(:,:,2);
                                            sm2=segment_ok(:,:,3);
                                        else if b==10
                                                sm1=segment_ok(:,:,2);
                                                sm2=segment_ok(:,:,4);
                                            else if b==11
                                                    sm1=segment_ok(:,:,2);
                                                    sm2=segment_ok(:,:,5);
                                                else if b==12
                                                        sm1=segment_ok(:,:,2);
                                                        sm2=segment_ok(:,:,6);
                                                    else if b==13
                                                            sm1=segment_ok(:,:,2);
                                                            sm2=segment_ok(:,:,7);
                                                        else if b==14
                                                                sm1=segment_ok(:,:,2);
                                                                sm2=segment_ok(:,:,8);
                                                            else if b==15
                                                                    sm1=segment_ok(:,:,2);
                                                                    sm2=segment_ok(:,:,9);
                                                                else if b==16
                                                                        sm1=segment_ok(:,:,3);
                                                                        sm2=segment_ok(:,:,4);
                                                                    else if b==17
                                                                            sm1=segment_ok(:,:,3);
                                                                            sm2=segment_ok(:,:,5);
                                                                        else if b==18
                                                                                sm1=segment_ok(:,:,3);
                                                                                sm2=segment_ok(:,:,6);
                                                                            else if b==19
                                                                                    sm1=segment_ok(:,:,3);
                                                                                    sm2=segment_ok(:,:,7);
                                                                                else if b==20
                                                                                        sm1=segment_ok(:,:,3);
                                                                                        sm2=segment_ok(:,:,8);
                                                                                    else if b==21
                                                                                            sm1=segment_ok(:,:,3);
                                                                                            sm2=segment_ok(:,:,9);
                                                                                        else if b==22
                                                                                                sm1=segment_ok(:,:,4);
                                                                                                sm2=segment_ok(:,:,5);
                                                                                            else if b==23
                                                                                                    sm1=segment_ok(:,:,4);
                                                                                                    sm2=segment_ok(:,:,6);
                                                                                                else if b==24
                                                                                                        sm1=segment_ok(:,:,4);
                                                                                                        sm2=segment_ok(:,:,7);
                                                                                                    else if b==25
                                                                                                            sm1=segment_ok(:,:,4);
                                                                                                            sm2=segment_ok(:,:,8);
                                                                                                        else if b==26
                                                                                                                sm1=segment_ok(:,:,4);
                                                                                                                sm2=segment_ok(:,:,9);
                                                                                                            else if b==27
                                                                                                                    sm1=segment_ok(:,:,5);
                                                                                                                    sm2=segment_ok(:,:,6);
                                                                                                                else if b==28
                                                                                                                        sm1=segment_ok(:,:,5);
                                                                                                                        sm2=segment_ok(:,:,7);
                                                                                                                    else if b==29
                                                                                                                            sm1=segment_ok(:,:,5);
                                                                                                                            sm2=segment_ok(:,:,8);
                                                                                                                        else if b==30
                                                                                                                                sm1=segment_ok(:,:,5);
                                                                                                                                sm2=segment_ok(:,:,9);
                                                                                                                            else if b==31
                                                                                                                                    sm1=segment_ok(:,:,6);
                                                                                                                                    sm2=segment_ok(:,:,7);
                                                                                                                                else if b==32
                                                                                                                                        sm1=segment_ok(:,:,6);
                                                                                                                                        sm2=segment_ok(:,:,8);
                                                                                                                                    else if b==33
                                                                                                                                            sm1=segment_ok(:,:,6);
                                                                                                                                            sm2=segment_ok(:,:,9);
                                                                                                                                        else if b==34
                                                                                                                                                sm1=segment_ok(:,:,7);
                                                                                                                                                sm2=segment_ok(:,:,8);
                                                                                                                                            else if b==35
                                                                                                                                                    sm1=segment_ok(:,:,7);
                                                                                                                                                    sm2=segment_ok(:,:,9);
                                                                                                                                                else if b==36
                                                                                                                                                        sm1=segment_ok(:,:,8);
                                                                                                                                                        sm2=segment_ok(:,:,9);
                                                                                                                                                    end
                                                                                                                                                end
                                                                                                                                            end
                                                                                                                                        end
                                                                                                                                    end
                                                                                                                                end
                                                                                                                            end
                                                                                                                        end
                                                                                                                    end
                                                                                                                end
                                                                                                            end
                                                                                                        end
                                                                                                    end
                                                                                                end
                                                                                            end
                                                                                        end
                                                                                    end
                                                                                end
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
x=cranio(:,:,ff);
[y1 x1]=find(sm1==1);
[y2 x2]=find(sm2==1);
stats=regionprops(sm1,'Centroid','MajorAxisLength','MinorAxisLength');
centsm1=[stats.Centroid];
stats=regionprops(sm2,'Centroid','MajorAxisLength','MinorAxisLength');
centsm2=[stats.Centroid]; 

figure,imshow(x,[]);
hold on
plot(centsm1(1),centsm1(2),'xr')
plot(centsm2(1),centsm2(2),'xg')

continua=input('For middle slice, was maxillary sinus correctly identified? (1) Yes (2) No (ff+) (3) No (ff-): \')

if continua==2
    ff=ff+1;
end
if continua==3
    ff=ff-1;
end
end

%%
segmentacao1=zeros(size(cranio));
segmentacao2=zeros(size(cranio));

% Upper slices (ff + 1) 
ok1=sm1;
ok2=sm2;
stats=regionprops(sm1,'Centroid','MajorAxisLength','MinorAxisLength');
centsm1=[stats.Centroid];
stats=regionprops(sm2,'Centroid','MajorAxisLength','MinorAxisLength');
centsm2=[stats.Centroid]; 

for aa=ff:size(cranio,3)
x=cranio(:,:,aa);
xlim=zeros(512);
for i = 1:size(x,1)
    for j = 1:size(x,2)
        if x(i,j)>=150; 
            xlim(i,j)=1;
        else xlim(i,j)=0;
        end
    end
end

%Morphological Operators
BW2=bwmorph(xlim,'bridge');
BW3 =  bwmorph(BW2,'majority');
BW3 =  bwmorph(BW3,'majority');
BW3 =  bwmorph(BW3,'majority');
SE = strel('line',3,0);
xlim2=imdilate(BW3,SE);
SE = strel('line',20,90);
xlim3=imdilate(xlim2,SE);
SE = strel('line',5,0);
xlim4=imdilate(xlim3,SE);
xlim5=bwmorph(xlim4,'bridge');

%Watershed
w=watershed(xlim5);

%Rule-based system
w(450:end,:)=1;
w2=zeros(512);
final=zeros(512);
w2(w~=0)=1;

a1=w2;
cc = bwconncomp(w2);
stats = regionprops(cc,'Area'); 
A = [stats.Area];
[~,biggest] = max(A);
a1(labelmatrix(cc)~=biggest) = 0;
w2(a1==1)=0; 

%Major areas
for i=1:9
    seg=zeros(512);
    A(biggest)=0;
    [~,biggest] = max(A);
    seg(labelmatrix(cc)== biggest) = 1;
    segment(:,:,i)=seg;
     
end
            
verificador1=0;
verificador2=0;
for i=1:9
    seg=segment(:,:,i);
  if  seg(round(centsm1(2)),round(centsm1(1))) == 1
      ok1=seg;
      verificador1=1;
  else  if seg(round(centsm2(2)),round(centsm2(1)))==1  
          ok2=seg;
          verificador2=1;
      end
  end
end

if (verificador1==0 || verificador2==0)
    break
end
se=strel('disk',3);
ok1=imerode(ok1,se);
ok2=imerode(ok2,se);
ok1=imfill(ok1,'holes');
ok2=imfill(ok2,'holes');

stats=regionprops(ok1,'Centroid','MajorAxisLength','MinorAxisLength');
centsm1=[stats.Centroid];
stats=regionprops(ok2,'Centroid','MajorAxisLength','MinorAxisLength');
centsm2=[stats.Centroid]; 

x=cranio(:,:,aa);
segmentacao1(:,:,aa)=ok1;
segmentacao2(:,:,aa)=ok2;

end

% Lower slices (ff - 1)
stats=regionprops(sm1,'Centroid','MajorAxisLength','MinorAxisLength');
centsm1=[stats.Centroid];
stats=regionprops(sm2,'Centroid','MajorAxisLength','MinorAxisLength');
centsm2=[stats.Centroid]; 
ok1=sm1;
ok2=sm2;

for aa=ff:-1:1
x=cranio(:,:,aa);
xlim=zeros(512);
for i = 1:size(x,1)
    for j = 1:size(x,2)
        if x(i,j)>=150; 
            xlim(i,j)=1;
        else xlim(i,j)=0;
        end
    end
end

%Morphological Operators
BW2=bwmorph(xlim,'bridge');
BW3 =  bwmorph(BW2,'majority');
BW3 =  bwmorph(BW3,'majority');
BW3 =  bwmorph(BW3,'majority');
SE = strel('line',10,90);
xlim2=imdilate(BW3,SE);
SE = strel('line',5,0);
xlim3=imdilate(xlim2,SE);
xlim4=bwmorph(xlim3,'bridge');

%Watershed
w=watershed(xlim4);

%Rule-based system
w(450:end,:)=1;
w2=zeros(512);
final=zeros(512);
w2(w~=0)=1;

a1=w2;
cc = bwconncomp(w2);
stats = regionprops(cc,'Area'); 
A = [stats.Area];
[~,biggest] = max(A);
a1(labelmatrix(cc)~=biggest) = 0;
w2(a1==1)=0; 

%Major areas
for i=1:9
    seg=zeros(512);
    A(biggest)=0;
    [~,biggest] = max(A);
    seg(labelmatrix(cc)== biggest) = 1;
    segment(:,:,i)=seg;
     
end
            
verificador1=0;
verificador2=0;
for i=1:9
    seg=segment(:,:,i);
  if  seg(round(centsm1(2)),round(centsm1(1))) == 1
      ok1=seg;
      verificador1=1;
  else  if seg(round(centsm2(2)),round(centsm2(1)))==1  
          ok2=seg;
          verificador2=1;
      end
  end
end


if (verificador1==0 || verificador2==0)
    break
end
se=strel('disk',3);
ok1=imerode(ok1,se);
ok2=imerode(ok2,se);

ok1=imfill(ok1,'holes');
ok2=imfill(ok2,'holes');

stats=regionprops(ok1,'Centroid','MajorAxisLength','MinorAxisLength');
centsm1=[stats.Centroid];
stats=regionprops(ok2,'Centroid','MajorAxisLength','MinorAxisLength');
centsm2=[stats.Centroid]; 

x=cranio(:,:,aa);
segmentacao1(:,:,aa)=ok1;
segmentacao2(:,:,aa)=ok2;

end

seg3d=segmentacao1+segmentacao2;

%% Threshold (air-free and involvement)
ar_full=[];
osso_full=[];
muco_full=[];

for i=1:size(cranio,3)
    
    im=cranio(:,:,i);
    seg_seio=seg3d(:,:,i);
    
    se=strel('disk',4);
    seg_seio=imerode(seg_seio,se);
    
    im(seg_seio~=1)=-2000;
    seio_real(:,:,i)=im;
    muco=zeros(512);
    ar=zeros(512);
    osso=zeros(512);
    for j = 1:512
        for k = 1:512
            if im(j,k)>1000;
                osso(j,k)=1;
            else if  im(j,k)>-200 &&  im(j,k)<1000;
                    muco(j,k)=1;
                else if im(j,k)>-1500 &&  im(j,k)<-200;
                        ar(j,k)=1;
                    end
                end
            end
        end
    end
    
ar_full(:,:,i)=ar;
osso_full(:,:,i)=osso;
muco_full(:,:,i)=muco;
end
 
%Volumetry

%Air-free volume 
cranio_sm=zeros(size(masccranio));
cranio_sm(:,:,inicio:fim)=seg3d;

ar_full_d=(ar_full);
ar_full_d(:,1:256,:)=0;
ar_full_e=(ar_full);
ar_full_e(:,257:512,:)=0;

vol_ar_mm3_d=sum(ar_full_d(:)).*(pxsp(1)*pxsp(2)*step*.3);
vol_ar_mm3_e=sum(ar_full_e(:)).*(pxsp(1)*pxsp(2)*step*.3);

%Involvement volume
muco_full_cranio=zeros(size(masccranio));
muco_full_cranio(:,:,inicio:fim)=muco_full;

vol_comp_mm3=sum(muco_full(:)).*(pxsp(1)*pxsp(2)*step*.3);

muco_full_d=(muco_full_cranio);
muco_full_d(:,1:256,:)=0;
muco_full_e=(muco_full_cranio);
muco_full_e(:,257:512,:)=0;

vol_muco_mm3_d=sum(muco_full_d(:)).*(pxsp(1)*pxsp(2)*step*.3);
vol_muco_mm3_e=sum(muco_full_e(:)).*(pxsp(1)*pxsp(2)*step*.3);

% Total volume 
vol_tot_mm3_e=vol_muco_mm3_e+vol_ar_mm3_e;
vol_tot_mm3_d=vol_muco_mm3_d+vol_ar_mm3_d;

%3D visualization
figure,isosurface(smooth3(smooth3(masccranio)));
hold on
alpha .2
isosurface(smooth3(smooth3(cranio_sm)));
hold on
alpha .4
isosurface(smooth3(smooth3(muco_full_cranio)));
grid on

results=[vol_tot_mm3_d; vol_tot_mm3_e; vol_muco_mm3_d;vol_muco_mm3_d];

disp(['Left total volume (mm³): ', num2str(vol_tot_mm3_e)])
disp(['Right total volume (mm³): ', num2str(vol_tot_mm3_d)])
disp(['Left involvement volume(mm³): ', num2str(vol_muco_mm3_e)])
disp(['Right involvement volume(mm³): ', num2str(vol_muco_mm3_d)])
