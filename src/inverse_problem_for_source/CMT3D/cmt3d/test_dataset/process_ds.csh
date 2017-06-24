echo 'Pre-processing data and/or syn ...'
process_trinet_syn_new.pl -S -m CMTSOLUTION -a station_file -s 20 -p -d syn_PP syn/*H?*
process_cal_data.pl -m CMTSOLUTION -s 20 -p -d data_PP data/*H?.sac
echo 'Padding zeros to data and syn ...'
pad_zeros.pl -s -l -26 data_PP/*H?.sac >/dev/null
pad_zeros.pl -s -l -26 syn_PP/*H?* >/dev/null
echo "Cutting ... 24 data files to line up with synthetics"
sac <<EOF
cut -25.005 n 4086
r data_PP/GSC.CI.BHE.sac syn_PP/GSC.CI.BHE syn_PP/GSC.CI.BHE.depth syn_PP/GSC.CI.BHE.latitude syn_PP/GSC.CI.BHE.longitude syn_PP/GSC.CI.BHE.Mpp syn_PP/GSC.CI.BHE.Mrp syn_PP/GSC.CI.BHE.Mrr syn_PP/GSC.CI.BHE.Mrt syn_PP/GSC.CI.BHE.Mtp syn_PP/GSC.CI.BHE.Mtt
cut off
w over
cut -25.005 n 4086
r data_PP/GSC.CI.BHN.sac syn_PP/GSC.CI.BHN syn_PP/GSC.CI.BHN.depth syn_PP/GSC.CI.BHN.latitude syn_PP/GSC.CI.BHN.longitude syn_PP/GSC.CI.BHN.Mpp syn_PP/GSC.CI.BHN.Mrp syn_PP/GSC.CI.BHN.Mrr syn_PP/GSC.CI.BHN.Mrt syn_PP/GSC.CI.BHN.Mtp syn_PP/GSC.CI.BHN.Mtt
cut off
w over
cut -25.005 n 4086
r data_PP/GSC.CI.BHZ.sac syn_PP/GSC.CI.BHZ syn_PP/GSC.CI.BHZ.depth syn_PP/GSC.CI.BHZ.latitude syn_PP/GSC.CI.BHZ.longitude syn_PP/GSC.CI.BHZ.Mpp syn_PP/GSC.CI.BHZ.Mrp syn_PP/GSC.CI.BHZ.Mrr syn_PP/GSC.CI.BHZ.Mrt syn_PP/GSC.CI.BHZ.Mtp syn_PP/GSC.CI.BHZ.Mtt
cut off
w over
cut -25.005 n 4088
r data_PP/LAF.CI.BHE.sac syn_PP/LAF.CI.BHE syn_PP/LAF.CI.BHE.depth syn_PP/LAF.CI.BHE.latitude syn_PP/LAF.CI.BHE.longitude syn_PP/LAF.CI.BHE.Mpp syn_PP/LAF.CI.BHE.Mrp syn_PP/LAF.CI.BHE.Mrr syn_PP/LAF.CI.BHE.Mrt syn_PP/LAF.CI.BHE.Mtp syn_PP/LAF.CI.BHE.Mtt
cut off
w over
cut -25.005 n 4083
r data_PP/LAF.CI.BHN.sac syn_PP/LAF.CI.BHN syn_PP/LAF.CI.BHN.depth syn_PP/LAF.CI.BHN.latitude syn_PP/LAF.CI.BHN.longitude syn_PP/LAF.CI.BHN.Mpp syn_PP/LAF.CI.BHN.Mrp syn_PP/LAF.CI.BHN.Mrr syn_PP/LAF.CI.BHN.Mrt syn_PP/LAF.CI.BHN.Mtp syn_PP/LAF.CI.BHN.Mtt
cut off
w over
cut -25.005 n 4086
r data_PP/LAF.CI.BHZ.sac syn_PP/LAF.CI.BHZ syn_PP/LAF.CI.BHZ.depth syn_PP/LAF.CI.BHZ.latitude syn_PP/LAF.CI.BHZ.longitude syn_PP/LAF.CI.BHZ.Mpp syn_PP/LAF.CI.BHZ.Mrp syn_PP/LAF.CI.BHZ.Mrr syn_PP/LAF.CI.BHZ.Mrt syn_PP/LAF.CI.BHZ.Mtp syn_PP/LAF.CI.BHZ.Mtt
cut off
w over
cut -25.005 n 4088
r data_PP/LFP.CI.BHE.sac syn_PP/LFP.CI.BHE syn_PP/LFP.CI.BHE.depth syn_PP/LFP.CI.BHE.latitude syn_PP/LFP.CI.BHE.longitude syn_PP/LFP.CI.BHE.Mpp syn_PP/LFP.CI.BHE.Mrp syn_PP/LFP.CI.BHE.Mrr syn_PP/LFP.CI.BHE.Mrt syn_PP/LFP.CI.BHE.Mtp syn_PP/LFP.CI.BHE.Mtt
cut off
w over
cut -25.955 n 4096
r data_PP/LFP.CI.BHN.sac syn_PP/LFP.CI.BHN syn_PP/LFP.CI.BHN.depth syn_PP/LFP.CI.BHN.latitude syn_PP/LFP.CI.BHN.longitude syn_PP/LFP.CI.BHN.Mpp syn_PP/LFP.CI.BHN.Mrp syn_PP/LFP.CI.BHN.Mrr syn_PP/LFP.CI.BHN.Mrt syn_PP/LFP.CI.BHN.Mtp syn_PP/LFP.CI.BHN.Mtt
cut off
w over
cut -25.955 n 4106
r data_PP/LFP.CI.BHZ.sac syn_PP/LFP.CI.BHZ syn_PP/LFP.CI.BHZ.depth syn_PP/LFP.CI.BHZ.latitude syn_PP/LFP.CI.BHZ.longitude syn_PP/LFP.CI.BHZ.Mpp syn_PP/LFP.CI.BHZ.Mrp syn_PP/LFP.CI.BHZ.Mrr syn_PP/LFP.CI.BHZ.Mrt syn_PP/LFP.CI.BHZ.Mtp syn_PP/LFP.CI.BHZ.Mtt
cut off
w over
cut -25.955 n 4100
r data_PP/MOP.CI.BHE.sac syn_PP/MOP.CI.BHE syn_PP/MOP.CI.BHE.depth syn_PP/MOP.CI.BHE.latitude syn_PP/MOP.CI.BHE.longitude syn_PP/MOP.CI.BHE.Mpp syn_PP/MOP.CI.BHE.Mrp syn_PP/MOP.CI.BHE.Mrr syn_PP/MOP.CI.BHE.Mrt syn_PP/MOP.CI.BHE.Mtp syn_PP/MOP.CI.BHE.Mtt
cut off
w over
cut -25.955 n 4107
r data_PP/MOP.CI.BHN.sac syn_PP/MOP.CI.BHN syn_PP/MOP.CI.BHN.depth syn_PP/MOP.CI.BHN.latitude syn_PP/MOP.CI.BHN.longitude syn_PP/MOP.CI.BHN.Mpp syn_PP/MOP.CI.BHN.Mrp syn_PP/MOP.CI.BHN.Mrr syn_PP/MOP.CI.BHN.Mrt syn_PP/MOP.CI.BHN.Mtp syn_PP/MOP.CI.BHN.Mtt
cut off
w over
cut -25.955 n 3943
r data_PP/MOP.CI.BHZ.sac syn_PP/MOP.CI.BHZ syn_PP/MOP.CI.BHZ.depth syn_PP/MOP.CI.BHZ.latitude syn_PP/MOP.CI.BHZ.longitude syn_PP/MOP.CI.BHZ.Mpp syn_PP/MOP.CI.BHZ.Mrp syn_PP/MOP.CI.BHZ.Mrr syn_PP/MOP.CI.BHZ.Mrt syn_PP/MOP.CI.BHZ.Mtp syn_PP/MOP.CI.BHZ.Mtt
cut off
w over
cut -25.955 n 4090
r data_PP/PAS.CI.BHE.sac syn_PP/PAS.CI.BHE syn_PP/PAS.CI.BHE.depth syn_PP/PAS.CI.BHE.latitude syn_PP/PAS.CI.BHE.longitude syn_PP/PAS.CI.BHE.Mpp syn_PP/PAS.CI.BHE.Mrp syn_PP/PAS.CI.BHE.Mrr syn_PP/PAS.CI.BHE.Mrt syn_PP/PAS.CI.BHE.Mtp syn_PP/PAS.CI.BHE.Mtt
cut off
w over
cut -25.955 n 4107
r data_PP/PAS.CI.BHN.sac syn_PP/PAS.CI.BHN syn_PP/PAS.CI.BHN.depth syn_PP/PAS.CI.BHN.latitude syn_PP/PAS.CI.BHN.longitude syn_PP/PAS.CI.BHN.Mpp syn_PP/PAS.CI.BHN.Mrp syn_PP/PAS.CI.BHN.Mrr syn_PP/PAS.CI.BHN.Mrt syn_PP/PAS.CI.BHN.Mtp syn_PP/PAS.CI.BHN.Mtt
cut off
w over
cut -25.955 n 4106
r data_PP/PAS.CI.BHZ.sac syn_PP/PAS.CI.BHZ syn_PP/PAS.CI.BHZ.depth syn_PP/PAS.CI.BHZ.latitude syn_PP/PAS.CI.BHZ.longitude syn_PP/PAS.CI.BHZ.Mpp syn_PP/PAS.CI.BHZ.Mrp syn_PP/PAS.CI.BHZ.Mrr syn_PP/PAS.CI.BHZ.Mrt syn_PP/PAS.CI.BHZ.Mtp syn_PP/PAS.CI.BHZ.Mtt
cut off
w over
cut -25.955 n 4082
r data_PP/PDE.CI.BHE.sac syn_PP/PDE.CI.BHE syn_PP/PDE.CI.BHE.depth syn_PP/PDE.CI.BHE.latitude syn_PP/PDE.CI.BHE.longitude syn_PP/PDE.CI.BHE.Mpp syn_PP/PDE.CI.BHE.Mrp syn_PP/PDE.CI.BHE.Mrr syn_PP/PDE.CI.BHE.Mrt syn_PP/PDE.CI.BHE.Mtp syn_PP/PDE.CI.BHE.Mtt
cut off
w over
cut -25.955 n 4107
r data_PP/PDE.CI.BHN.sac syn_PP/PDE.CI.BHN syn_PP/PDE.CI.BHN.depth syn_PP/PDE.CI.BHN.latitude syn_PP/PDE.CI.BHN.longitude syn_PP/PDE.CI.BHN.Mpp syn_PP/PDE.CI.BHN.Mrp syn_PP/PDE.CI.BHN.Mrr syn_PP/PDE.CI.BHN.Mrt syn_PP/PDE.CI.BHN.Mtp syn_PP/PDE.CI.BHN.Mtt
cut off
w over
cut -25.955 n 4106
r data_PP/PDE.CI.BHZ.sac syn_PP/PDE.CI.BHZ syn_PP/PDE.CI.BHZ.depth syn_PP/PDE.CI.BHZ.latitude syn_PP/PDE.CI.BHZ.longitude syn_PP/PDE.CI.BHZ.Mpp syn_PP/PDE.CI.BHZ.Mrp syn_PP/PDE.CI.BHZ.Mrr syn_PP/PDE.CI.BHZ.Mrt syn_PP/PDE.CI.BHZ.Mtp syn_PP/PDE.CI.BHZ.Mtt
cut off
w over
cut -25.955 n 4107
r data_PP/RPV.CI.BHE.sac syn_PP/RPV.CI.BHE syn_PP/RPV.CI.BHE.depth syn_PP/RPV.CI.BHE.latitude syn_PP/RPV.CI.BHE.longitude syn_PP/RPV.CI.BHE.Mpp syn_PP/RPV.CI.BHE.Mrp syn_PP/RPV.CI.BHE.Mrr syn_PP/RPV.CI.BHE.Mrt syn_PP/RPV.CI.BHE.Mtp syn_PP/RPV.CI.BHE.Mtt
cut off
w over
cut -25.955 n 4107
r data_PP/RPV.CI.BHN.sac syn_PP/RPV.CI.BHN syn_PP/RPV.CI.BHN.depth syn_PP/RPV.CI.BHN.latitude syn_PP/RPV.CI.BHN.longitude syn_PP/RPV.CI.BHN.Mpp syn_PP/RPV.CI.BHN.Mrp syn_PP/RPV.CI.BHN.Mrr syn_PP/RPV.CI.BHN.Mrt syn_PP/RPV.CI.BHN.Mtp syn_PP/RPV.CI.BHN.Mtt
cut off
w over
cut -25.955 n 4106
r data_PP/RPV.CI.BHZ.sac syn_PP/RPV.CI.BHZ syn_PP/RPV.CI.BHZ.depth syn_PP/RPV.CI.BHZ.latitude syn_PP/RPV.CI.BHZ.longitude syn_PP/RPV.CI.BHZ.Mpp syn_PP/RPV.CI.BHZ.Mrp syn_PP/RPV.CI.BHZ.Mrr syn_PP/RPV.CI.BHZ.Mrt syn_PP/RPV.CI.BHZ.Mtp syn_PP/RPV.CI.BHZ.Mtt
cut off
w over
cut -25.955 n 4076
r data_PP/TOV.CI.BHE.sac syn_PP/TOV.CI.BHE syn_PP/TOV.CI.BHE.depth syn_PP/TOV.CI.BHE.latitude syn_PP/TOV.CI.BHE.longitude syn_PP/TOV.CI.BHE.Mpp syn_PP/TOV.CI.BHE.Mrp syn_PP/TOV.CI.BHE.Mrr syn_PP/TOV.CI.BHE.Mrt syn_PP/TOV.CI.BHE.Mtp syn_PP/TOV.CI.BHE.Mtt
cut off
w over
cut -25.955 n 4107
r data_PP/TOV.CI.BHN.sac syn_PP/TOV.CI.BHN syn_PP/TOV.CI.BHN.depth syn_PP/TOV.CI.BHN.latitude syn_PP/TOV.CI.BHN.longitude syn_PP/TOV.CI.BHN.Mpp syn_PP/TOV.CI.BHN.Mrp syn_PP/TOV.CI.BHN.Mrr syn_PP/TOV.CI.BHN.Mrt syn_PP/TOV.CI.BHN.Mtp syn_PP/TOV.CI.BHN.Mtt
cut off
w over
cut -25.955 n 4047
r data_PP/TOV.CI.BHZ.sac syn_PP/TOV.CI.BHZ syn_PP/TOV.CI.BHZ.depth syn_PP/TOV.CI.BHZ.latitude syn_PP/TOV.CI.BHZ.longitude syn_PP/TOV.CI.BHZ.Mpp syn_PP/TOV.CI.BHZ.Mrp syn_PP/TOV.CI.BHZ.Mrr syn_PP/TOV.CI.BHZ.Mrt syn_PP/TOV.CI.BHZ.Mtp syn_PP/TOV.CI.BHZ.Mtt
cut off
w over
quit
EOF
rm -f tmpPadZeroToSyn.sac
echo 'Bandpassing data and/or synthetics ...'
process_trinet_syn_new.pl -S -t 6/30 -d syn_T006_T030 syn_PP/*.?H?*
rotate.pl syn_T006_T030/*.?H[E1]*
process_cal_data.pl -t 6/30 -d data_T006_T030 -x d data_PP/*.?H?.sac
rotate.pl data_T006_T030/*.?HE.sac.d
