folder_list=[""];
dir_folder="G:\ShuttleDrive\Nlx2NRD\E1"; % Directory which parameter file saved
fs=32000;
for ii=1:length(folder_list)
    noise=[

    ];
    ADC=[
        46,45,44,43, ... %TT1%
        42,41,38,37, ... %TT2%
        36,35,34,33, ... %TT3%
        32,48,49,50, ... %TT4%
        51,52,53,54, ... %TT5%
        55,56,40,39, ... %TT6%
        57,58,59,60, ... %TT7%
        61,62,63,47      %TT8%
    ];
    map=[1:32];
    nrdfilename=dir(fullfile(folder_list(ii),'*.nrd'));
    path1=folder_list(ii);
    path1=char(path1)
    nrdfilename1=nrdfilename.name;
    csv='\geom.csv';
    json='\params.json';
    prb='\probe.prb';
    mda='\raw.mda';

    %% Reset 'Samples' array
    [Samples1,Header]=Nlx2MatNRD(strcat(path1,nrdfilename1),1,[0,1],1,1,[]);
    reclen=length(Samples1);
    recduration=(reclen/fs)/60
    Samples=zeros(length(map),reclen);

    %% Export csv
    M=csvread(strcat(dir_folder,'\geom_sd.csv'));
    M=M(map,:);
    csvwrite(strcat(path1,csv),M);

    %% Copy json file
    status=copyfile(strcat(dir_folder,'\params.json'),strcat(path1,json));

    %% Export prb
    prb_output=fopen(strcat(path1,prb),'w');
    tt=[1:8];
    csc=0;
    m_index=1;
    map_index=1;
    fprintf(prb_output,"channel_groups={\n");
    for i=1:length(tt)
        fprintf(prb_output,"%d:\n{\n'channels':[",i-1);
        for n=1:4
            if n~=4
                fprintf(prb_output,'%d,',csc);
                csc=csc+1;
            else
                fprintf(prb_output,"%d],\n'geometry':[",csc);
                csc=csc+1
            end
        end
        for n=1:4
            if n~=4
                fprintf(prb_output,'[%d,%d]',M(m_index,1),M(m_index,2));
                m_index=m_index+1;
            else
                fprintf(prb_output,"[%d,%d]],\n'label':[",M(m_index,1),M(m_index,2));
                m_index=m_index+1;
            end
        end
        for n=1:4
            if n~=4
                fprintf(prb_output,"ch_%d',",map(map_index));
                map_index=map_index+1;
            else
                fprintf(prb_output,"ch_%d']\n},\n",map(map_index));
                map_index=map_index+1;
            end
        end
    end
    fprintf(prb_output,"}");
    fclose(prb_output)

    %% Export NRD
    h=waitbar(0,'Please wait...!');
    for i=1:length(map)
        a=ADC(map(i));
        if ismember(a,noise)
            channel=i
            CSC=map(i)
            Samples(i,:)=zeros;
            waitbar(i/length(map))
        else
            channel=i
            CSC=map(i)
            [Samples1,Header]=Nlx2MatNRD(strcat(path1,nrdfilename1),a,[0,1],1,1,[]);
            ADBitVolts=split(char(Header(14)));
            ADBitVolts=cell2mat(ADBitVolts(2));
            Samples1=1000000*Samples1*str2num(ADBitVolts);
            Samples(i,:)=Samples1;
            waitbar(i/length(map))
        end
    end
    close(h)

    writemda(Samples,strcat(path1,mda),'int16');
    
    % sound(sin(1:2800*3));
end