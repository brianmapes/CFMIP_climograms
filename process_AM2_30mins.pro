;;; process_AM2_30mins,'../Datasets/AM2/Coar_eee_30mins.nc', 'COARE', [10,11,0,1]
;;; process_AM2_30mins,'../Datasets/AM2/OLDER/Coar_eee.2000.nc', 'COARE', [10,11,0,1]

pro process_AM2_30mins, filename, station_string, months

;;; Read the time-height data all into memory
print, 'Reading file ', filename

@grab_AM2_column.30mins.pro
;@grab_AM2_column.idl

T = temp
TEMP = reform(T(0,0,*,*)) 
TEMP = TEMP - TEMP(*,0) # (1+fltarr(n_elements(time)))

print, 'Now plot it up'

;;;Time is in days

Jan1 = 0
Feb1 = (where(fix(time) eq 31) )(0)
Mar1 = (where(fix(time) eq 31+28) )(0)
Apr1 = (where(fix(time) eq 31+28+31) )(0)
May1 = (where(fix(time) eq 31+28+31+30) )(0)
Jun1 = (where(fix(time) eq 31+28+31+30+31) )(0)
Jul1 = (where(fix(time) eq 31+28+31+30+31+30) )(0)
Aug1 = (where(fix(time) eq 31+28+31+30+31+30+31) )(0)
Sep1 = (where(fix(time) eq 31+28+31+30+31+30+31+31) )(0)
Oct1 = (where(fix(time) eq 31+28+31+30+31+30+31+31+30) )(0)
Nov1 = (where(fix(time) eq 31+28+31+30+31+30+31+31+30+31) )(0)
Dec1 = (where(fix(time) eq 31+28+31+30+31+30+31+31+30+31+30) )(0)
EOY = n_elements(time) -1
dates = [jan1,feb1,mar1,apr1,may1,jun1,jul1,aug1,sep1,oct1,nov1,dec1,EOY]

;for imo = 0,11 do begin
for iimo = 0,n_elements(months)-1 do begin
    IMO = MONTHS[iIMO]

    time1 = dates[imo]
    time2 = dates[imo+1]

;;;;;;;;;;;;;;;;;;;;;;;;;;;; variables for plot: GFDL names

tplot = time(time1:time2)
pplot = pfull ;; hPa
NLEVS = n_elements(pplot)
ntimes = n_elements(tplot)
plevs = [4,6,8,11,14,18] ;;; a nice set for correlations

LWTOA  = reform( olr    (0,0,time1:time2) )
LWTOAC = reform( olr_clr(0,0,time1:time2) )
SWTOA = reform( (swup_toa -swdn_toa) (0,0,time1:time2) ) ;;; positive upward; vlaue negative
SWTOAC= reform( (swup_toa_clr -swdn_toa_clr) (0,0,time1:time2) )
SWUPTOA = reform( swup_toa (0,0,time1:time2) )
SWUPTOAC= reform( swup_toa_clr (0,0,time1:time2) ) ;; small positive number

SWCF = SWTOAC-SWTOA ;;; a negative number, reflected solar
LWCF = LWTOAC-LWTOA ;;; a positive number, reduced OLR
momeanCF = mean(SWCF+LWCF)

QSW = reform( tdt_sw(0,0,*,time1:time2) ) *86400.
QLW = reform( tdt_lw(0,0,*,time1:time2) ) *86400.
;Q1plot = reform( (tdt_sw + tdt_lw + tdt_conv + tdt_ls + tdt_vdif) (0,0,*,time1:time2) ) *86400.
;Q2plot = reform( (qdt_conv + qdt_ls + qdt_vdif) (0,0,*,time1:time2) ) *86400.
Q1plot = reform( (tdt_sw + tdt_lw + tdt_conv + tdt_ls ) (0,0,*,time1:time2) ) *86400.
Q2plot = reform( (qdt_conv + qdt_ls) (0,0,*,time1:time2) ) *86400.

MAXCWC = 0.5 ;; g/kg
CWC = reform( (liq_wat + ice_wat) (0,0,*,time1:time2) ) *1000 ;;; g/kg
CWClevels = [.0033,0.01,0.03,0.09,0.27,0.81] *MAXCWC ;;; g/kg
CWClev_label = 'log contours'

;;; After KW256
CWClevels = [0.002, 0.01*( 1+indgen(10) )*2]
CWClev_label = '2,20,40,60... mg/kg'

RHplot = reform( rh(0,0,*,time1:time2) )
RHlevels = 10*findgen(11)

TPERT = Temp(*,time1:time2)
TPERT = TPERT - TPERT(*,0) # (1+fltarr(time2-time1+1))

OMEG = reform(omega(0,0,*,time1:time2)) * 86400./100. ;; Pa/s -> hPa/d
OM = 300 ;;  
print, 'Omega scale hPa/d', max(abs(OMEG))
OMLEVELS = -OM + OM*findgen(21)/10. 

;;; Divergence from omega
DIV = OMEG
for it=0,n_elements(tplot)-1 do $
  div(*,it) = -deriv( OMEG(*,it) )/deriv(pplot) *1.e6/86400. ;; units 1e-6 s-1

hlp, div
Uplot = reform(ucomp(0,0,*,time1:time2))
Vplot = reform(vcomp(0,0,*,time1:time2))

rain = reform( precip(0,0,time1:time2) ) *3600. ;;; mm/h rather than mm/s
rainLS = reform( prec_ls(0,0,time1:time2) )  *3600. ;;; mm/h rather than mm/s
momeanrain = mean(rain)
cloudfraction = reform( cld_amt(0,0,*,time1:time2) )
ccloudfraction = cloudfraction*0

MODEL_NAME = 'AM2_30'
SMOO = 49
LEV500 = 9 ;;; level where p ~500, for Bony label
;;;;; these need adding for pageplots
VERSION = 'p?'
site = station_string

;;;;;;;;;;; Filename for postscript:

monames = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
mostr = monames(imo)

PSdirstub = '../PSFILES/PS_'+MODEL_NAME +'/'
namestub = MODEL_NAME + '.' + station_string + '.' + monames(imo) 
statstub = '.rain'+str(fix(momeanrain*24)) + 'mmd.CRF' + str(fix(momeanCF)) +'W.OM'+ str(fix(mean(omeg(LEV500,*))))+'hpad'

fname = PSdirstub+namestub+statstub+'.ps'
print, fname

stop

tall
device, file=fname

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; MODEL
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; INDEPENDENT
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;PLOTTING CODE EXTERNAL/SHARED
@plots_page_gfdl.idl
;@plots_page_lowertrop.idl
;@plots_regress_nocurd.idl
Psout

;@rain_composites_only.idl

endfor ;;; IMONTH LOOP

end
