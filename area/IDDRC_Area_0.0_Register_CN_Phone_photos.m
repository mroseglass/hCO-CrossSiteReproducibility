% all are saved as RGB 8-bit, CN is nearly double the pixel and you will need to downsample

%try something, both seem equivalent but "nearest" is supposed to be faster
%CNresize_halfnearest = imresize(CN1P1CA,0.5,"nearest");
%CNresize_halfbox = imresize(CN1P1CA,0.5,"box");
%imwrite(CNresize_halfbox,"CNresize_halfbox.jpg");
%imwrite(CNresize_halfnearest,"CNresize_halfnearest.jpg");

tform = load("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/PiecewiseLinearTransformation2D_CN_44ControlPoint.mat")
t = tform.t;
fixed  = imread('/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/RawImageFiles/FullHema_CNtoRegister-01-01-01-02.png');
fixed  = im2gray(fixed);
fixed = imresize(fixed, 0.2409);
%cd("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/RawImageFiles/CNRawPhonePhotos/CN R1 Day7")
%DirCNR1toR4 = pwd;
%CN_orignal = dir(fullfile(DirCNR1toR4,"*.jpg"));

rootdir = '/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/RawImageFiles/IDDRC_D35/CN R4 Day35';
filelist = dir(fullfile(rootdir, '**/*.*'));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]); 

for k = 1:length(filelist)
    Org_file = filelist(k).name;
    Resized = imread(strcat(filelist(k).folder,"/",filelist(k).name));
    %Resized = imresize(ResizeReady,0.5,"nearest");
    Rfixed = imref2d(size(fixed));
    registered = imwarp(Resized,t,OutputView=Rfixed);
    NewName = strcat("CN_Registered_",extractAfter(filelist(k).folder,70),"_",Org_file);
    imwrite(Resized,strcat("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/1221CNRegisteredImages/CN_scRNAseq_D84/",NewName))
end

Resized = rgb2gray(Resized);
registered = rgb2gray(registered);
imshowpair(Resized,registered,"blend")
imshowpair(Resized,registered,"diff")
imshowpair(A,B,"diff")