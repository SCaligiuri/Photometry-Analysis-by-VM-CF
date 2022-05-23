clear
close all
data= csvread('Gpr151Cre_Oxy_singleInj_M3_270520_2.csv', 2, 0); % import data and skip first 2 header rows

%DORIC ampling rate : 12kS = 12000 samples/sec
%Thor sampling rate : 1017.252625 
%factor between them : 11.79648

%get rid of the data before the demodulation (not a number before
%approx.80ms
startdemodul=data(:,2)>0;
datademodul=data(startdemodul,:);
recording=datademodul;
clear 'datademodul';
samplingrate=12000

rawtime=recording(:,1); %1 table for time
raw405=recording(:,2);%1 table for 405 channel
raw465=recording(:,3);%1 table for 465 channel
rawall=recording(:,4);%1 table for all inputs (specific to DORIC, if > 5 V, bad recording)

%downsampling
drecording = downsample(recording,1200);

dtime=drecording(:,1); %1 table for time
d405=drecording(:,2);%1 table for 405 channel
d465=drecording(:,3);%1 table for 465 channel
dall=drecording(:,4);%1 table for all inputs (specific to DORIC, if > 5 V, bad recording)

sec=12190
dsec=12190/1200
dsecr=round(dsec)


% to plot the raw signal (all the channels)        
figure(999); hold on; grid on;
plot(rawtime, rawall, 'm')
title('rawall all session')
legend('all/doric')

%is the recording ok (<5V for raw all)
limit=5;
if any(rawall > limit);
   
    disp('wrong recording');
    f = msgbox('wrong recording, detector saturated!');
else;
    disp('good recording');
    f = msgbox('good recording, keep analyzing!');
end;

%to calculate a simple df and Df/f without polyfit on channel 405
df =(raw465-raw405);%delta F (465-405) without any correction
dff = df./raw405;%delta F (465-405)/405 without any correction

%down sample
downdf =(d465-d405);%delta F (465-405) without any correction
downdff = downdf./d405;%delta F (465-405)/405 without any correction
%to normalize using Zscore 
zdff=zscore(dff);

figure(1); hold on; grid on;
subplot(221)
plot(rawtime, raw405, 'r:')
title('raw405 all session')
legend('405')
subplot(222)
plot(rawtime, raw465, 'g')
title('raw465 all session')
legend('465')
subplot(223)
plot(rawtime, df, 'g')  
title('deltaF without time correction')
legend('deltaF')
subplot(224)
plot(rawtime, dff, 'm');  
title('deltaF without correction');
legend('deltaF/F');
hline(0,'k:');hold on;

%downsample
figure(998); hold on; grid on;
subplot(221)
plot(dtime, d405, 'r:')
title('raw405 all session down')
legend('405')
subplot(222)
plot(dtime, d465, 'g')
title('raw465 all session down')
legend('465')
subplot(223)
plot(dtime,downdff, 'g')
title('deltaf/F')
legend('deltaf/F')

%si qqchose s'est passe (mettre le temps en sec ici, debile)

secsaline=187;
saline=round(secsaline*dsec);
secnic=323;
nic=round(secnic*dsec);
endsec=1247;
endpts=round(endsec*dsec);

%dff alone
plot(dtime,downdff, 'g')
title('deltaf/F')
legend('deltaf/F')
vline(secsaline,'b','saline');
vline(sec,'r', '2252 20 mg/kg');
vline(endsec,'r', 'end');

sizerecord=size(recording(:,1));
lenghtrecord=sizerecord/sec
sizedrecordd=size(drecording(:,1));
lenghtrecordd=sizedrecordd/dsec

%%
%all recording     
     d405c=controlFit(d405, d465 );

dfc_all = (d465-d405c);              
dffc_all=(dfc_all./ d405c); 
zdffall=zscore(dffc_all);

plot(dtime,dffc_all, 'g')
title('deltaf/F session')
legend('deltaf/F 1 session')
vline(secsaline,'b','saline');
vline(secnic,'r', 'nic low dose 0.125mg.kg-1');
vline(endsec,'r', 'end');

%focus on different parts of this recording
%1 min before saline with the down sample
befsaline= dtime(saline-60*dsecr:saline, :);   
d405befsaline = d405(saline-60*dsecr:saline, :);  
     d465befsaline = d465(saline-60*dsecr:saline, :);% extract GCamp signal for individual epoch'
     
     d405cbefsaline=controlFit(d405befsaline, d465befsaline );

dfc_befsaline = (d465befsaline-d405cbefsaline);              
dffc_befsaline=(dfc_befsaline./d405cbefsaline); 
zscore_befsaline=zdffall(saline-60*dsecr:saline, :);

plot(befsaline,dffc_befsaline, 'g')
title('deltaf/F 1 min before saline')
legend('deltaf/F 1 min before saline')

%event 1 after saline min
aftsaline= dtime(saline:saline+60*dsecr, :);   
d405aftsaline = d405(saline:saline+60*dsecr,:);  
     d465aftsaline = d465(saline:saline+60*dsecr, :);% extract GCamp signal for individual epoch'
     
     d405caftsaline=controlFit(d405aftsaline, d465aftsaline );
 
dfc_aftsaline = (d465aftsaline-d405caftsaline);              
dffc_aftsaline=(dfc_aftsaline./d405caftsaline); 
zscore_aftsaline=zdffall(saline:saline+60*dsecr,:);

plot(aftsaline,dffc_aftsaline, 'g')
title('deltaf/F 1 min after saline')
legend('deltaf/F 1 min after saline')

%1 min before nic with the down sample
befnic= dtime(nic-60*dsecr:nic, :);   
d405befnic = d405(nic-60*dsecr:nic, :);  
     d465befnic = d465(nic-60*dsecr:nic, :);% extract GCamp signal for individual epoch'
     
     d405cbefnic=controlFit(d405befnic, d465befnic );
 
dfc_befnic = (d465befnic-d405cbefnic);              
dffc_befnic=(dfc_befnic./d405cbefnic); 
zscore_befnic=zdffall(nic-60*dsecr:nic, :);
 
plot(befnic,dffc_befnic, 'g')
title('deltaf/F 1 min before nic 0.125mg.kg-1')
legend('deltaf/F 1 min before nic 0.125mg.kg-1')

%1 min after nic with the down sample
aftnic= dtime(nic:nic+60*dsecr, :);   
d405aftnic = d405(nic:nic+60*dsecr, :);  
     d465aftnic = d465(nic:nic+60*dsecr, :);% extract GCamp signal for individual epoch'
     
     d405caftnic=controlFit(d405aftnic, d465aftnic );
 
dfc_aftnic = (d465aftnic-d405caftnic);              
dffc_aftnic=(dfc_aftnic./d405caftnic); 
zscore_aftnic=zdffall(nic:nic+60*dsecr, :); 
 
plot(aftnic,dffc_aftnic, 'g')
title('deltaf/F 1 min after nic 0.125mg.kg-1')
legend('deltaf/F 1 min after nic 0.125mg.kg-1')


%5 min after nic with the down sample
aft5nic= dtime(nic:nic+300*dsecr, :);   
d405aft5nic = d405(nic:nic+300*dsecr, :);  
     d465aft5nic = d465(nic:nic+300*dsecr, :);% extract GCamp signal for individual epoch'
     
     d405caft5nic=controlFit(d405aft5nic, d465aft5nic );
 
dfc_aft5nic = (d465aft5nic-d405caft5nic);              
dffc_aft5nic=(dfc_aft5nic./d405caft5nic); 
zscore_aft5nic=zdffall(nic:nic+300*dsecr, :);  
 
plot(aft5nic,dffc_aft5nic, 'g')
title('deltaf/F 1 min aft5min nic 0.125mg.kg-1')
legend('deltaf/F 1 min aft5min nic 0.125mg.kg-1')

% min 5 a 6 after nic with the down sample
aft5_6nic= dtime(nic+300*dsecr:nic+360*dsecr, :);   
d405aft5_6nic = d405(nic+300*dsecr:nic+360*dsecr, :);  
     d465aft5_6nic = d465(nic+300*dsecr:nic+360*dsecr, :);% extract GCamp signal for individual epoch'
     
     d405caft5_6nic=controlFit(d405aft5_6nic, d465aft5_6nic );
 
dfc_aft5_6nic = (d465aft5_6nic-d405caft5_6nic);              
dffc_aft5_6nic=(dfc_aft5_6nic./d405caft5_6nic); 
zscore_aft5_6nic=zdffall(nic+300*dsecr:nic+360*dsecr, :);   
 
plot(aft5_6nic,dffc_aft5_6nic, 'g')
title('deltaf/F 1 min aft5_6ore nic 0.125 mg.kg-1')
legend('deltaf/F 1 min aft5_6ore nic 0.125 mg.kg-1')

% min 10 a 11 after nic with the down sample
aft10_11nic= dtime(nic+600*dsecr:nic+660*dsecr, :);   
d405aft10_11nic = d405(nic+600*dsecr:nic+660*dsecr, :);  
     d465aft10_11nic = d465(nic+600*dsecr:nic+660*dsecr, :);% extract GCamp signal for individual epoch'
     
     d405caft10_11nic=controlFit(d405aft10_11nic, d465aft10_11nic );
 
dfc_aft10_11nic = (d465aft10_11nic-d405caft10_11nic);              
dffc_aft10_11nic=(dfc_aft10_11nic./d405caft10_11nic); 
zscore_aft10_11nic=zdffall(nic+600*dsecr:nic+660*dsecr, :);   
 
plot(aft10_11nic,dffc_aft10_11nic, 'g')
title('deltaf/F 1 min aft10_11ore nic 0.125 mg.kg-1')
legend('deltaf/F 1 min aft10_11ore nic 0.125 mg.kg-1')

% min 14 a 15 after nic with the down sample (name doesnt fit purpose)
aft15_16nic= dtime(nic+840*dsecr:nic+900*dsecr, :);   
d405aft15_16nic = d405(nic+840*dsecr:nic+900*dsecr, :);  
     d465aft15_16nic = d465(nic+840*dsecr:nic+900*dsecr, :);% extract GCamp signal for individual epoch'
     
     d405caft15_16nic=controlFit(d405aft15_16nic, d465aft15_16nic );
 
dfc_aft15_16nic = (d465aft15_16nic-d405caft15_16nic);              
dffc_aft15_16nic=(dfc_aft15_16nic./d405caft15_16nic); 
zscore_aft15_16nic=zdffall(nic+840*dsecr:nic+900*dsecr, :);   
 
plot(aft15_16nic,dffc_aft15_16nic, 'g')
title('deltaf/F 1 min aft15_16ore nic 0.125 mg.kg-1')
legend('deltaf/F 1 min aft15_16ore nic 0.125 mg.kg-1')


% min 5 a 10 after nic with the down sample
aft5_10nic= dtime(nic+300*dsecr:nic+600*dsecr, :);   
d405aft5_10nic = d405(nic+300*dsecr:nic+600*dsecr, :);  
     d465aft5_10nic = d465(nic+300*dsecr:nic+600*dsecr, :);% extract GCamp signal for individual epoch'
     
     d405caft5_10nic=controlFit(d405aft5_10nic, d465aft5_10nic );
 
dfc_aft5_10nic = (d465aft5_10nic-d405caft5_10nic);              
dffc_aft5_10nic=(dfc_aft5_10nic./d405caft5_10nic); 
zscore_aft5_10nic=zdffall(nic+300*dsecr:nic+600*dsecr, :);    
 
plot(aft5_10nic,dffc_aft5_10nic, 'g')
title('deltaf/F 1 min aft5_10ore nic 0.125 mg.kg-1')
legend('deltaf/F 1 min aft5_10ore nic 0.125 mg.kg-1')

% min 10 a 15 after nic with the down sample
aft10_15nic= dtime(nic+600*dsecr:nic+900*dsecr, :);   
d405aft10_15nic = d405(nic+600*dsecr:nic+900*dsecr, :);  
     d465aft10_15nic = d465(nic+600*dsecr:nic+900*dsecr, :);% extract GCamp signal for individual epoch'
     
     d405caft10_15nic=controlFit(d405aft10_15nic, d465aft10_15nic );
 
dfc_aft10_15nic = (d465aft10_15nic-d405caft10_15nic);              
dffc_aft10_15nic=(dfc_aft10_15nic./d405caft10_15nic); 
zscore_aft10_15nic=zdffall(nic+600*dsecr:nic+900*dsecr, :);  
 
plot(aft10_15nic,dffc_aft10_15nic, 'g')
title('deltaf/F 1 min aft10_15ore nic 0.125 mg.kg-1')
legend('deltaf/F 1 min aft10_15ore nic 0.125 mg.kg-1')



%%
%generate table all session
resultsall_session=table(dtime,dffc_all, zdffall, 'VariableNames', {'time_all_session', 'deltaF_F_session', 'z_score'});
 
resultsbymin=table(befsaline,dffc_befsaline,zscore_befsaline, aftsaline,dffc_aftsaline,zscore_aftsaline,befnic,dffc_befnic,zscore_befnic,aftnic,dffc_aftnic,zscore_aftnic,aft5_6nic,dffc_aft5_6nic, zscore_aft5_6nic,aft10_11nic,dffc_aft10_11nic, zscore_aft10_11nic,aft15_16nic,dffc_aft15_16nic, zscore_aft15_16nic, 'VariableNames', {'before_saline', 'deltaF_F_before_saline', 'zscore_before_saline','after_saline', 'deltaF_F_after_saline','zscore_after_saline','before_nic', 'deltaF_F_before_nic','zscore_before_nic','after_nic', 'deltaF_F_after_nic','zscore_after_nic','after5_6_nic', 'deltaF_F_after5_6_nic','zscore_after5_6_nic','after10_11_nic', 'deltaF_F_after10_11_nic','zscore_after10_11_nic','T15_16','deltaF_15_16', 'zscore_15_16'});

 
resultsby5min=table(aft5nic,dffc_aft5nic,zscore_aft5nic,aft5_10nic,dffc_aft5_10nic,zscore_aft5_10nic,aft10_15nic,dffc_aft10_15nic,zscore_aft10_15nic,'VariableNames', {'time_nic0_5', 'deltaF_F_0_5min_nic','zscore_0_5min_nic','time_nic5_10', 'deltaF_F_5_10min_nic','zscore_5_10min_nic','time_nic10_15', 'deltaF_F_10_15min_nic','zscore_10_15min_nic'});


