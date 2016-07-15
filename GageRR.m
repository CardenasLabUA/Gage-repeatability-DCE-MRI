clear all; close all; clc;

%% Set up time components and determine AIF with Cp10 function
time=0/60:6.4/60:5.4; time=time';
Cp=Cp10(time);

%% Simulate 30 values for ktrans and ve
ktrans=.1*rand(30,1)+0.25;
ve=.1*rand(30,1)+0.4;

%% Generate enhancement curves
for j=1:30
    Ct_tumor(:,j)=myToftsCt2(ktrans(j),ktrans(j)/ve(j),time,Cp)';
end

%% Generate enhancement curve for reference region and add white gaussian noise
S_RR=myToftsCt2(.1,1,time,Cp)';
S_RR=awgn(S_RR,40,'measured');

%% Define guesses for non-linear analysis
lb=zeros(3,1); lb=lb';
ub=ones(3,1)*10; ub=ub';

%% Set initial SNR

%% Set up info for semi-quant analysis
NumBaselineImages=10;
    
%% Loop through analysis changing SNR by 1 each time
mySNR=5:1:40;
for SNR=1:length(mySNR)   
    %% Define number of reptitions to test 30 curves for 3 days by gage analysis
    for Reps=1:1000
        %%Create 3 simulated enhancement curves by adding white gaussian noise to 
        %%original enhancement curve
        for q=1:3    
            for r=1:30
                S_toi=awgn(Ct_tumor(:,r),mySNR(SNR),'measured');
                [maxValue,index] = max(S_toi);
                MER(r,q) = maxValue;
                TTP(r,q) = time(index);
                iauc64(r,q) = trapz(time(1:11),S_toi(1:11));
                slope(r,q) = MER(r,q)./TTP(r,q);
                B=fitdcemri(S_toi,S_RR,time,'nonneg');
                RKtrans1_linear_nonNeg(r,q)=B(1);
                x0=abs(randn(3,1));
                %x0=B(1:3); Uncomment if you want to use initial guesses
                            %for NRRM from LRRM to make NRRM*
                B = fitdcemri(S_toi,S_RR,time,x0,lb,ub,'RRM');
                RKtrans_nonlinear_nonNeg(r,q)=B(1);
                B=fitdcemri(S_toi,S_RR,time,'lsq');
                RKtrans2_linear_lsq(r,q)=B(1);
            end
        end
        
        m=1:30; m=m';
        Parts=[m,m,m]; Parts=[Parts(:)];
        operators=[ones(1,90)]';
        
        %%Perform gage analysis for semi-quant analysis
        observations_MER=[MER(:)];
        observations_TTP=[TTP(:)];
        observations_iauc64=[iauc64(:)];
        observations_slope=[slope(:)];
        T=gagerr( observations_MER,{Parts,operators},'printtable','off','printgraph','off'); 
        MER_vector(Reps,SNR)=T(2,2);
        T=gagerr( observations_TTP,{Parts,operators},'printtable','off','printgraph','off'); 
        TTP_vector(Reps,SNR)=T(2,2);
        T=gagerr( observations_iauc64,{Parts,operators},'printtable','off','printgraph','off'); 
        iauc64_vector(Reps,SNR)=T(2,2);
        T=gagerr( observations_slope,{Parts,operators},'printtable','off','printgraph','off'); 
        slope_vector(Reps,SNR)=T(2,2);
        
        %%Perform gage analysis for linear methods
        observations=[RKtrans1_linear_nonNeg(:);RKtrans2_linear_lsq(:)];
        T=gagerr( observations(1:90),{Parts(1:90),operators(1:90)},'printtable','off','printgraph','off'); 
        Nonneg(Reps,SNR)=T(2,2);
        T=gagerr( observations(91:end),{Parts(91:end),operators(91:end)},'printtable','off','printgraph','off');
        Lsq(Reps,SNR)=T(2,2);
   
        %%Perform gage analysis for nonlinear method
        observations2=[RKtrans_nonlinear_nonNeg(:)];
        T=gagerr( observations2(1:end),{Parts(1:end),operators(1:end)},'printtable','off','printgraph','off');
        Nonlinear(Reps,SNR)=T(2,2);
    end
end