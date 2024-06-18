%% 计算 UVW，并设定频率观测通道：
% pro uvgen,U,V
% Nant=75
% ANT_POS=FLTARR(Nant,3)
% BSL=FLTARR(Nant,Nant,3)
% DXYZ=FLTARR(Nant,Nant,3)
% V=FLTARR(Nant,Nant)
% U=FLTARR(Nant,Nant)
% LAT=42.+12.710/60
% LAT=!PI LAT/182
% LONT=115.+15.030/60.
% LOC_TIME=LONT/15.
% OPENR,LUN1,'ANT_POS.TXT',/GET_LUN;%当地坐标 xyz
% FOR I=0,Nant-1 DO BEGIN
% READF,LUN1,X,Y,Z
% ANT_POS[I, ]=[Y,X,Z]
% ENDFOR
% FREE_LUN,LUN1
% FOR I=0,Nant-1 DO BEGIN;
% FOR J=0,Nant-1 DO BEGIN
% BSL[I,J, ]=ANT_POS[I, ]-ANT_POS[J, ]
% DXYZ[I,J,0]=-SIN(LAT) BSL[I,J,0]+COS(LAT) BSL[I,J,2]
% DXYZ[I,J,1]=BSL[I,J,1]
% DXYZ[I,J,2]=COS(LAT) BSL[I,J,0]+SIN(LAT) BSL[I,J,2]
% ENDFOR
% ENDFOR
% H=0.
% DEC=0.
% DEC=DEC/180. !PI
% H=H/180. !PI
% FOR I=0,Nant-1 DO BEGIN
% FOR J=0,Nant-1 DO BEGIN
% U[I,J]=DXYZ[I,J,0] SIN(H)+DXYZ[I,J,1] COS(H)
% V[I,J]=DXYZ[I,J,0] (-
% SIN(DEC)) Cos(H)+DXYZ[I,J,1] SIN(DEC) Sin(H)+COS(DEC) DXYZ[I,J,2]
% ENDFOR
% ENDFOR
% freq=12.0d   ;Ghz
% lamda=0.3/freq
% u=u/lamda
% v=v/lamda;%表示在特定的频率下，从源的角度看过去。
% END



%% 读图并设定图片格式部分代码：
pro read_map,source,filename;;%把原图读成数组文件存下来
loadct,3;
device,decomposed=0
file=['ifa141107.png','ifa140817.png','ifa140824.png','ifa160321.png','ifa160905.png']
name=file[3] ;% 可选择上面不同的图像
cd,'sources'
filename=strmid(name,0,9)
map1=read_png(name,R,G,B);

source=dblarr(512,512)
source[0:510,40:471]=map1[0:510,40:471];
histo=histogram(source);%直方图变换
factor=(where(histo[20:n_elements(histo)-1] eq max(histo[20:n_elements(histo)-
    1]))+20)
Flux_max=(256d/factor)+3.
print,flux_max,factor,10^(max(source)/255. ((4.8d)-3)+3)
window,xs=512,ys=512
tvscl,source;%画图，把最大值变成最亮，最小值最暗
write_jpeg,filename+'_source.jpg',TVRD(TRUE=3),TRUE=3;%把窗口图写成文件
cd,'..'
end

%% 密度加权和脏图脏束部分代码：
pro observation
read_map,source,filename
uvgen,U,V
w=1   ;%0  为自然权（natural），1 为均一权（uniform）
pix=512
bpix=pix 2;%脏束大小是脏图的 2 倍；
SRC_ft=fft(source,/CENTER);%原图进行傅里叶变换

Fov=41.9D/60/180 !pi;把 UV 值对应到 512 512 的数组里
maxUV=Pix/Fov
Cor_data=fltarr(60,60,2)
uv_sample=complexarr(pix,pix)
uv_beam=fltarr(bpix,bpix)
WEIGHT=FLTARR(pix,pix)+1.

FOR I=0,60-1 DO BEGIN
FOR J=0,60-1 DO BEGIN
if i ne j then begin
    X=U[I,J]/MAXuv PIX+PIX/2
    Y=V[I,J]/MAXuv PIX+PIX/2
    IF X GE PIX-1 THEN x=0
    IF X LE -PIX+1 THEN x=0
    IF Y GE PIX-1 THEN y=0
    IF Y LE -PIX+1 THEN y=0
    cor_data[i,j,0]=real_part(src_ft[FIX(X),FIX(Y)])
    cor_data[i,j,1]=imaginary(src_ft[FIX(X),FIX(Y)])
    uv_sample[fix(x),fix(y)]=complex(cor_data[i,j,0],cor_data[i,j,1])
    uv_beam[fix(x),fix(y)]=1.0
    if w eq 1 then weight[fix(x),fix(y)]=weight[fix(x),fix(y)]+1.0 ;
        ENDIF
        ENDFOR
        ENDFOR
        print,max(weight)
        dimage=abs(fft(uv_sample/weight))
        dbeam=abs(fft(uv_beam))
        dbeam=shift(dbeam,pix,pix)
        dimage=rotate(dimage,2)
        dbeam=rotate(dbeam,2)
        WINDOW,1,XS=512,YS=512
        TVSCL,dimage
        write_jpeg,filename+'_dimage_'+strcompress(w,/remove_all)+'.jpg',TVRD(TRUE=3),T
        RUE=3
        window,2,xs=1024,ys=1024
        tvscl,dbeam
        write_jpeg,filename+'_dbeam_'+strcompress(w,/remove_all)+'.jpg',TVRD(TRUE=3),T
        RUE=3
        save,dimage,dbeam,filename,w,filename=filename+'_data_'+strcompress(w,/remove_all
        )+'.sav';%脏图脏束都有了的情况下，存到 SAV 中
    end
    
    
    
    %% Högbom CLEAN  部分代码：
    pro HogbomClean;,beam,dirtyimage,image,residual
    w=1
    RESTORE,'ifa141107_data_'+strcompress(w,/remove_all)+'.SAV';
    
    DEVICE, DECOMPOSED = 0
    LOADCT,3
    PIX=512
    BM=dBEAM/MAX(dBEAM)
    residual=dimage/MAX(dimage)
    
    CEN_PIX=30
    CEN_BM=BM[PIX-CEN_PIX/2:PIX+CEN_PIX/2,PIX-
        CEN_PIX/2:PIX+CEN_PIX/2]
    CLEAN_BM=GAUSS2DFIT(CEN_BM,A);%做二维高斯拟合，拟合脏束的宽度，利用这个宽度生成洁束
    PRINT,A[2] 3600/PIX,' Arcsec',A[3] 3600/PIX,' Arcsec',A[0],A[6]
    X=FINDGEN(PIX)-PIX/2
    Y=FINDGEN(PIX)-PIX/2
    UU=FLTARR(PIX,PIX)
    
    XX=X COS(A[6])-Y SIN(A[6])
    YY=X SIN(A[6])-Y COS(A[6])
    
    FOR I=0,PIX-1 DO BEGIN
    FOR J=0,PIX-1 DO BEGIN
    UU[I,J]=(XX[I]/A[2] 3)^2+(YY[J]/A[3] 3)^2
    ENDFOR
    ENDFOR
    
    C_BM=A[1] EXP(-UU/2); % 生成最终的洁束
    
    MODEL=FLTARR(PIX,PIX);  %记录最强值位置和强度的数组
    region=intarr(pix,pix)+1.;  %在哪些区域找最大值
    
    ITERATION=0;%循环次数，这里没用该方式
    GAIN=0.03
    stop_criterion=2.0;%看残图的方差来决定，看是否该停止
    
    while (max(residual region) gt stop_criterion stddev(residual region)) do begin;%观察脏图中的最大值，若大于某个值就继续找，小于就停止
        LOC=WHERE(residual EQ MAX(residual region));%在区域中找最大值，记录下位置
        IF (N_ELEMENTS(LOC) GT 1) THEN BEGIN
        LOCXY=LOC[0];%若返回多个值，取第一个即可
        ENDIF ELSE BEGIN
        LOCXY=LOC
        ENDELSE
        
        LOCY=LOCXY/PIX;%确定返回值对应的 X，Y 具体是多少
        LOCX=LOCXY MOD PIX
        GAIN_residual=MAX(residual) GAIN;%剪掉剩余值的 gain
        MODEL[LOCX,LOCY]=MODEL[LOCX,LOCY]+GAIN_residual;
        residual=residual-BM[pix/2-1-LOCX+pix/2:pix/2-1+pix-1-LOCX+pix/2,pix/2-1-
            LOCY+pix/2:pix/2-1+pix-1-LOCY+pix/2] GAIN_residual;    %从残图中减去对应的beam 的位置
        ITERATION=ITERATION+1;%记录循环次数
        ENDWHILE
        PRINT,ITERATION,MAX(residual),MIN(residual)
        
        WINDOW,0,XS=PIX,YSIZE=PIX
        
        my_convol,model,c_bm,c_img
        c_img=c_img/max(c_img)+residual; 0.25;+residual 0.25;%model 和洁束做卷积，再加回残图
        tvscl,c_img
        write_jpeg,filename+'CleanImage_'+strcompress(w,/remove_all)+'.jpg',tvrd(true=3),true=3;
        save,c_img,filename=filename+'CleanImage_'+strcompress(w,/remove_all)+'.sav'
    end
    
    pro my_convol,inputimage,funct,outputimage
    nsize=size(inputimage)
    pix=nsize[1]
    
    in_ft=fft(inputimage,/center)
    funct_ft=fft(funct,/center)
    
    out_ft=in_ft funct_ft
    outputimage=real_part(fft(out_ft,/inverse,/center))
    outputimage=shift(outputimage,pix/2,pix/2) pix pix
end






