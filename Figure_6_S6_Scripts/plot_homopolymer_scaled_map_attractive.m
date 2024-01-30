clear all;

constructs={'NUCKS','NUCKSdelK','POLR1F'};
mytitles={'NUCKS','NUCKSdelK','POLR1F'};

reps=[5];


types{1}={'RK'};
types{2}={'DE'};
mycolor3=[[24 154 214]/255; [225 57 55]/255];

colorhex={'#008181','#FEE011','#B4DCB6','#E51B4C','#F9BEBE','#25266B','#863D97','#808080','#4E64AF','#73CDD8','#BCD631','#F58230','#7F8133','#D9BDDB','#BB55A0','#9A6427'};
for i=1:length(colorhex)
    mycolor(i,:)=sscanf(colorhex{i}(2:end),'%2x%2x%2x',[1 3])/255;
end
mycolor2=mycolor(end:-1:1,:);

[hs, seqs]=fastaread('../Figure_6_S6_Data/nucleolar_idrs.fasta');


%figure;
for c=1:length(constructs)
    pos=find(strcmp(hs,constructs{c})==1);
    myseq=seqs{pos};

    chargepos=zeros(1,length(seqs{1}));
    for i=1:length(types)
        tmp=types{i};
        tmp=tmp{1};
        for j=1:length(tmp)
            pos=strfind(myseq,tmp(j));
            chargepos(pos)=i;
        end
    end

    count2=0;
    for r=[1,2,3,4,5]%1:reps(c)
        count2=count2+1;
        smapr(:,:,count2)=importdata(['../Figure_6_S6_Data/' constructs{c} '/340/' num2str(r) '/ana/homopolymer_scaled_map.csv']);   
    end

    smap=mean(smapr,3);

    for i=1:size(smap,1)
        for j=i:size(smap,1)
            smap(j,i)=smap(i,j);
        end
    end
    
    figure;
    subplot(4,4,[6 7 8 10 11 12 14 15 16])
    matones=ones(length(myseq),length(myseq));
    %imagesc(smap-matones); hold on; 
    %colormap(bluewhitered)
    %caxis([-0.5 0.4])
    imagesc(smap); hold on; 
    colormap(jet)
    caxis([0.5 1.5])
    %colorbar
    title(constructs{c})

    for i=1:length(types)
        tmp=types{i};
        tmp=tmp{1};
        if i==1
            for j=1:length(tmp)
                pos=strfind(myseq,tmp(j));
                plot(size(smap,2),pos,'o','markerfacecolor',mycolor3(i,:),'markeredgecolor','k','markersize',6); hold on;
                plot(pos,size(smap,2),'o','markerfacecolor',mycolor3(i,:),'markeredgecolor','k','markersize',6); hold on;
            end
        elseif i==2
            for j=1:length(tmp)
                pos=strfind(myseq,tmp(j));
                plot(pos,1,'o','markerfacecolor',mycolor3(i,:),'markeredgecolor','k','markersize',6); hold on;
                plot(1,pos,'o','markerfacecolor',mycolor3(i,:),'markeredgecolor','k','markersize',6); hold on;
            end
        end
    end	

    %% plot net charge per residue
    bsz=5;
    for i=1:length(myseq)-bsz+1
        ncpr(i)=(length(find(chargepos(i:i+bsz-1)==1))-length(find(chargepos(i:i+bsz-1)==2)))/bsz;
    end
    for i=1:length(myseq)
        if i>=bsz & i<=length(myseq)-bsz
            mncpr(i)=mean(ncpr(i-bsz+1:i));
        elseif i==1
            mncpr(i)=ncpr(i);
        elseif i<bsz
            mncpr(i)=mean(ncpr(1:i));
        elseif i>length(myseq)-bsz
            tmp=length(myseq)-i;
            mncpr(i)=mean(ncpr(end-tmp:end));
        end
    end

    subplot(4,4,2:4)
    plot(1:1:length(myseq),mncpr);
    % Extract positive and negative part
    yp = (mncpr + abs(mncpr))/2;
    yn = (mncpr - abs(mncpr))/2;
    % Plot the data using area function
    area(1:1:length(myseq),yp,'FaceColor',[24 154 214]/255)
    hold on
    area(1:1:length(myseq),yn,'FaceColor',[225 57 55]/255)
    xlim([1 length(myseq)]);  
    ylim([-1 1]); 

    subplot(4,4,[5 9 13])
    plot(1:1:length(myseq),mncpr);
    % Plot the data using area function
    area(1:1:length(myseq),yp,'FaceColor',[24 154 214]/255)
    hold on
    area(1:1:length(myseq),yn,'FaceColor',[225 57 55]/255)
    xlim([1 length(myseq)]); 
    ylim([-1 1]); 
    view(90,90)


    %% add where blocks are
    %alphaval=0.3
    %if c==1
    %    posblock=[[34, 54], [66, 74], [89, 96]];
    %    negblock=[[75, 79], [98, 104], [115, 120], [123, 130]];
    %    for i=1:length(posblock)/2
    %        rectangle('Position',[posblock(2*(i-1)+1), 0, posblock(2*(i-1)+2)-posblock(2*(i-1)+1), length(myseq)+1],'EdgeColor',[0,0,1,alphaval],'FaceColor',[0,0,1,alphaval]);
    %    end
    %    for i=1:length(negblock)/2
    %        rectangle('Position',[0, negblock(2*(i-1)+1), length(myseq)+1, negblock(2*(i-1)+2)-negblock(2*(i-1)+1)],'EdgeColor',[1,0,0,alphaval],'FaceColor',[1,0,0,alphaval])
    %    end
    %end
    
    mncpr
    %if c==2
    %    return
    %end

    %print -painters -depsc 

    clear smapr; clear smap; clear ncpr; clear mncpr; clear yp; clear yn; clear chargepos; 
end

%print -painters -depsc 'dye_weighted_homopolymer_scaled_distance_map.eps'

