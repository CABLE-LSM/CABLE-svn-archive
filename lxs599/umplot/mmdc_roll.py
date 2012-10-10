import cdms2, sys
import numpy as np

cfile=cdms2.open('MeanMnthDailyCycles_%syrs.nc' % sys.argv[1])
tas = cfile['tas']
qh  = cfile['hfss']
qle = cfile['hfls']
apr = cfile['field5226']
rnt = cfile['field3333']

# Manaus-------------------------------------------------------------------------------------------

tas_ma=np.roll(tas[:,:,0],-7,axis=1)
qh_ma =np.roll(qh[:,:,0],-7,axis=1)
qle_ma=np.roll(qle[:,:,0],-7,axis=1)
apr_ma=np.roll(apr[:,:,0],-7,axis=1)
rnt_ma=np.roll(rnt[:,:,0],-7,axis=1) # changed from -8 to -7

cfout5=cdms2.createDataset('mmdc_manaus_roll_%syrs.nc' % sys.argv[1])
cfout5.write(tas_ma,id=tas.id)
cfout5.write(qh_ma,id=qh.id)
cfout5.write(qle_ma,id=qle.id)
cfout5.write(apr_ma,id=apr.id)
cfout5.write(rnt_ma,id=rnt.id)
cfout5.sync()

# Hay-----------------------------------------------------------------------------------------------

tas_ha=np.roll(tas[:,:,1],20,axis=1)
qh_ha=np.roll(qh[:,:,1],20,axis=1)
qle_ha=np.roll(qle[:,:,1],20,axis=1)
apr_ha=np.roll(apr[:,:,1],20,axis=1)
rnt_ha=np.roll(rnt[:,:,1],20,axis=1)
#=np.roll([:,:,1],20,axis=1)

cfout6=cdms2.createDataset('mmdc_hay_roll_%syrs.nc' % sys.argv[1])
cfout6.write(tas_ha,id=tas.id)
cfout6.write(qh_ha,id=qh.id)
cfout6.write(qle_ha,id=qle.id)
cfout6.write(apr_ha,id=apr.id)
cfout6.write(rnt_ha,id=rnt.id)
cfout6.sync()

# Daly Waters-------------------------------------------------------------------------------------------

tas_dw=np.roll(tas[:,:,2],19,axis=1)
qh_dw=np.roll(qh[:,:,2],19,axis=1)
qle_dw=np.roll(qle[:,:,2],19,axis=1)
apr_dw=np.roll(apr[:,:,2],19,axis=1)
rnt_dw=np.roll(rnt[:,:,2],19,axis=1)
#=np.roll([:,:,2],19,axis=1)

cfout10=cdms2.createDataset('mmdc_dalywate_roll_%syrs.nc' % sys.argv[1])
cfout10.write(tas_dw,id=tas.id)
cfout10.write(qh_dw,id=qh.id)
cfout10.write(qle_dw,id=qle.id)
cfout10.write(apr_dw,id=apr.id)
cfout10.write(rnt_dw,id=rnt.id)
cfout10.sync()

# BATS-------------------------------------------------------------------------------------------

tas_ba=np.roll(tas[:,:,3],2,axis=1)
qh_ba=np.roll(qh[:,:,3],2,axis=1)
qle_ba=np.roll(qle[:,:,3],2,axis=1)
apr_ba=np.roll(apr[:,:,3],2,axis=1)
rnt_ba=np.roll(rnt[:,:,3],2,axis=1)
#=np.roll([:,:,3],2,axis=1)

cfout11=cdms2.createDataset('mmdc_hapexbat_roll_%syrs.nc' % sys.argv[1])
cfout11.write(tas_ba,id=tas.id)
cfout11.write(qh_ba,id=qh.id)
cfout11.write(qle_ba,id=qle.id)
cfout11.write(apr_ba,id=apr.id)
cfout11.write(rnt_ba,id=rnt.id)
cfout11.sync()

# India-------------------------------------------------------------------------------------------

tas_in=np.roll(tas[:,:,4],11,axis=1)
qh_in=np.roll(qh[:,:,4],11,axis=1)
qle_in=np.roll(qle[:,:,4],11,axis=1)
apr_in=np.roll(apr[:,:,4],11,axis=1)
rnt_in=np.roll(rnt[:,:,4],11,axis=1)
#=np.roll([:,:,4],11,axis=1)

cfout7=cdms2.createDataset('mmdc_india_roll_%syrs.nc' % sys.argv[1])
cfout7.write(tas_in,id=tas.id)
cfout7.write(qh_in,id=qh.id)
cfout7.write(qle_in,id=qle.id)
cfout7.write(apr_in,id=apr.id)
cfout7.write(rnt_in,id=rnt.id)
cfout7.sync()

# Hyytiala------------------------------------------------------------------------------------------------

tas_hy=np.roll(tas[:,:,5],4,axis=1)
qh_hy=np.roll(qh[:,:,5],4,axis=1)
qle_hy=np.roll(qle[:,:,5],4,axis=1)
apr_hy=np.roll(apr[:,:,5],4,axis=1)
rnt_hy=np.roll(rnt[:,:,5],4,axis=1)
#=np.roll([:,:,5],4,axis=1)

cfout0=cdms2.createDataset('mmdc_hyytiala_roll_%syrs.nc' % sys.argv[1])
cfout0.write(tas_hy,id=tas.id)
cfout0.write(qh_hy,id=qh.id)
cfout0.write(qle_hy,id=qle.id)
cfout0.write(apr_hy,id=apr.id)
cfout0.write(rnt_hy,id=rnt.id)
cfout0.sync()

# Afirca N---------------------------------------------------------------------------------------------

tas_an=np.roll(tas[:,:,6],2,axis=1)
qh_an  =np.roll(qh[:,:,6],2,axis=1)
qle_an=np.roll(qle[:,:,6],2,axis=1)
apr_an=np.roll(apr[:,:,6],2,axis=1)
rnt_an=np.roll(rnt[:,:,6],2,axis=1)
#=np.roll([:,:,6],2,axis=1)

cfout1=cdms2.createDataset('mmdc_africaN_roll_%syrs.nc' % sys.argv[1])
cfout1.write(tas_an,id=tas.id)
cfout1.write(qh_an,id=qh.id)
cfout1.write(qle_an,id=qle.id)
cfout1.write(apr_an,id=apr.id)
cfout1.write(rnt_an,id=rnt.id)
cfout1.sync()

# Africa S---------------------------------------------------------------------------------------------

tas_as=np.roll(tas[:,:,7],2,axis=1)
qh_as  =np.roll(qh[:,:,7],2,axis=1)
qle_as=np.roll(qle[:,:,7],2,axis=1)
apr_as=np.roll(apr[:,:,7],2,axis=1)
rnt_as=np.roll(rnt[:,:,7],2,axis=1)
#=np.roll([:,:,7],2,axis=1)

cfout21=cdms2.createDataset('mmdc_africaS_roll_%syrs.nc' % sys.argv[1])
cfout21.write(tas_as,id=tas.id)
cfout21.write(qh_as,id=qh.id)
cfout21.write(qle_as,id=qle.id)
cfout21.write(apr_as,id=apr.id)
cfout21.write(rnt_as,id=rnt.id)
cfout21.sync()

# Tumbarumba-------------------------------------------------------------------------------------------

tas_tu=np.roll(tas[:,:,8],20,axis=1)
qh_tu=np.roll(qh[:,:,8],20,axis=1)
qle_tu=np.roll(qle[:,:,8],20,axis=1)
apr_tu=np.roll(apr[:,:,8],20,axis=1)
rnt_tu=np.roll(rnt[:,:,8],20,axis=1)
#=np.roll([:,:,8],20,axis=1)

cfout2=cdms2.createDataset('mmdc_tumbarum_roll_%syrs.nc' % sys.argv[1])
cfout2.write(tas_tu,id=tas.id)
cfout2.write(qh_tu,id=qh.id)
cfout2.write(qle_tu,id=qle.id)
cfout2.write(apr_tu,id=apr.id)
cfout2.write(rnt_tu,id=rnt.id)
cfout2.sync()

# Tharandt-------------------------------------------------------------------------------------------

tas_th=np.roll(tas[:,:,9],2,axis=1)
qh_th=np.roll(qh[:,:,9],2,axis=1)
qle_th=np.roll(qle[:,:,9],2,axis=1)
apr_th=np.roll(apr[:,:,9],2,axis=1)
rnt_th=np.roll(rnt[:,:,9],2,axis=1)
#=np.roll([:,:,9],2,axis=1)

cfout4=cdms2.createDataset('mmdc_tharandt_roll_%syrs.nc' % sys.argv[1])
cfout4.write(tas_th,id=tas.id)
cfout4.write(qh_th,id=qh.id)
cfout4.write(qle_th,id=qle.id)
cfout4.write(apr_th,id=apr.id)
cfout4.write(rnt_th,id=rnt.id)
cfout4.sync()

# Bondville-------------------------------------------------------------------------------------------

tas_bo=np.roll(tas[:,:,10],-9,axis=1)
qh_bo=np.roll(qh[:,:,10],-9,axis=1)
qle_bo=np.roll(qle[:,:,10],-9,axis=1)
apr_bo=np.roll(apr[:,:,10],-9,axis=1)
rnt_bo=np.roll(rnt[:,:,10],-9,axis=1) # changed from -12 to -10 to -9
#=np.roll([:,:,10],-12,axis=1)

cfout3=cdms2.createDataset('mmdc_bondvill_roll_%syrs.nc' % sys.argv[1])
cfout3.write(tas_bo,id=tas.id)
cfout3.write(qh_bo,id=qh.id)
cfout3.write(qle_bo,id=qle.id)
cfout3.write(apr_bo,id=apr.id)
cfout3.write(rnt_bo,id=rnt.id)
cfout3.sync()

# Dinghushan-------------------------------------------------------------------------------------------

tas_di=np.roll(tas[:,:,11],15,axis=1)
qh_di=np.roll(qh[:,:,11],15,axis=1)
qle_di=np.roll(qle[:,:,11],15,axis=1)
apr_di=np.roll(apr[:,:,11],15,axis=1)
rnt_di=np.roll(rnt[:,:,11],15,axis=1)
#=np.roll([:,:,11],15,axis=1)

cfout9=cdms2.createDataset('mmdc_dinghush_roll_%syrs.nc' % sys.argv[1])
cfout9.write(tas_di,id=tas.id)
cfout9.write(qh_di,id=qh.id)
cfout9.write(qle_di,id=qle.id)
cfout9.write(apr_di,id=apr.id)
cfout9.write(rnt_di,id=rnt.id)
cfout9.sync()

# Little Washita-------------------------------------------------------------------------------------------

tas_lw=np.roll(tas[:,:,12],-10,axis=1)
qh_lw=np.roll(qh[:,:,12],-10,axis=1)
qle_lw=np.roll(qle[:,:,12],-10,axis=1)
apr_lw=np.roll(apr[:,:,12],-10,axis=1)
rnt_lw=np.roll(rnt[:,:,12],-10,axis=1)
#=np.roll([:,:,12],-13,axis=1)

cfout8=cdms2.createDataset('mmdc_litwashi_roll_%syrs.nc' % sys.argv[1])
cfout8.write(tas_lw,id=tas.id)
cfout8.write(qh_lw,id=qh.id)
cfout8.write(qle_lw,id=qle.id)
cfout8.write(apr_lw,id=apr.id)
cfout8.write(rnt_lw,id=rnt.id)
cfout8.sync()

# Walker Branch-------------------------------------------------------------------------------------------

tas_wb=np.roll(tas[:,:,13],-8,axis=1)
qh_wb=np.roll(qh[:,:,13],-8,axis=1)
qle_wb=np.roll(qle[:,:,13],-8,axis=1)
apr_wb=np.roll(apr[:,:,13],-8,axis=1)
rnt_wb=np.roll(rnt[:,:,13],-8,axis=1)
#=np.roll([:,:,13],-11,axis=1)

cfout12=cdms2.createDataset('mmdc_walkerbr_roll_%syrs.nc' % sys.argv[1])
cfout12.write(tas_wb,id=tas.id)
cfout12.write(qh_wb,id=qh.id)
cfout12.write(qle_wb,id=qle.id)
cfout12.write(apr_wb,id=apr.id)
cfout12.write(rnt_wb,id=rnt.id)
cfout12.sync()

# Loobos-------------------------------------------------------------------------------------------

tas_lo=np.roll(tas[:,:,14],2,axis=1)
qh_lo=np.roll(qh[:,:,14],2,axis=1)
qle_lo=np.roll(qle[:,:,14],2,axis=1)
apr_lo=np.roll(apr[:,:,14],2,axis=1)
rnt_lo=np.roll(rnt[:,:,14],2,axis=1)
#=np.roll([:,:,14],,axis=1)

cfout14=cdms2.createDataset('mmdc_loobos_roll_%syrs.nc' % sys.argv[1])
cfout14.write(tas_lo,id=tas.id)
cfout14.write(qh_lo,id=qh.id)
cfout14.write(qle_lo,id=qle.id)
cfout14.write(apr_lo,id=apr.id)
cfout14.write(rnt_lo,id=rnt.id)
cfout14.sync()

# Harvard Forest------------------------------------------------------------------------------------------

tas_hv=np.roll(tas[:,:,15],-7,axis=1)
qh_hv=np.roll(qh[:,:,15],-7,axis=1)
qle_hv=np.roll(qle[:,:,15],-7,axis=1)
apr_hv=np.roll(apr[:,:,15],-7,axis=1)
rnt_hv=np.roll(rnt[:,:,15],-7,axis=1)
#=np.roll([:,:,15],,axis=1)

cfout13=cdms2.createDataset('mmdc_harvard_roll_%syrs.nc' % sys.argv[1])
cfout13.write(tas_hv,id=tas.id)
cfout13.write(qh_hv,id=qh.id)
cfout13.write(qle_hv,id=qle.id)
cfout13.write(apr_hv,id=apr.id)
cfout13.write(rnt_hv,id=rnt.id)
cfout13.sync()

# Boreas NSA-------------------------------------------------------------------------------------------

tas_nb=np.roll(tas[:,:,16],-10,axis=1)
qh_nb=np.roll(qh[:,:,16],-10,axis=1)
qle_nb=np.roll(qle[:,:,16],-10,axis=1)
apr_nb=np.roll(apr[:,:,16],-10,axis=1)
rnt_nb=np.roll(rnt[:,:,16],-10,axis=1)
#=np.roll([:,:,16],,axis=1)

cfout16=cdms2.createDataset('mmdc_nsabor_roll_%syrs.nc' % sys.argv[1])
cfout16.write(tas_nb,id=tas.id)
cfout16.write(qh_nb,id=qh.id)
cfout16.write(qle_nb,id=qle.id)
cfout16.write(apr_nb,id=apr.id)
cfout16.write(rnt_nb,id=rnt.id)
cfout16.sync()

# Vielsalm-------------------------------------------------------------------------------------------

tas_vi=np.roll(tas[:,:,17],2,axis=1)
qh_vi=np.roll(qh[:,:,17],2,axis=1)
qle_vi=np.roll(qle[:,:,17],2,axis=1)
apr_vi=np.roll(apr[:,:,17],2,axis=1)
rnt_vi=np.roll(rnt[:,:,17],2,axis=1)
#=np.roll([:,:,17],,axis=1)

cfout15=cdms2.createDataset('mmdc_vielsalm_roll_%syrs.nc' % sys.argv[1])
cfout15.write(tas_vi,id=tas.id)
cfout15.write(qh_vi,id=qh.id)
cfout15.write(qle_vi,id=qle.id)
cfout15.write(apr_vi,id=apr.id)
cfout15.write(rnt_vi,id=rnt.id)
cfout15.sync()


cfile.close()

