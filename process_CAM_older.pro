;;; Run with 
;;; process_CAM,'../Datasets/CAM/eul64x128_amip-156E_2S.2000.nc', 'COARE', [0,1,10,11]

pro process_CAM, filename, station_string, months

;;; Read the time-height data all into memory
print, 'Reading CAM file ', filename

@grab_CAM_column.idl
halflevels = ilev

;;;;;;; CUrrently, NCAR 20 min output has value1-0-0-value2-0-0-value3-0-0. 
;;;;;;; Fix that so the values are sustained (as used in the model)

;;;;;;; CUrrently, NCAR 20 min output has value1-0-0-value2-0-0-value3-0-0. 
;;;;;;; Fix that so the values are sustained (as used in the model)
index = findgen(n_elements(time)/3)
FLNT = reform(flnt(0,0,*))  
     FLNT [index*3+2] = FLNT [index*3] & FLNT [index*3+1] = FLNT [index*3]
FLNTC= reform(flntc(0,0,*)) 
     FLNTC[index*3+2] = FLNTC[index*3] & FLNTC[index*3+1] = FLNTC[index*3]
fsntoa = reform(fsntoa(0,0,*))  
     fsntoa [index*3+2] = fsntoa [index*3] & fsntoa [index*3+1] = fsntoa [index*3]
fsntoaC= reform(fsntoac(0,0,*)) 
     fsntoaC[index*3+2] = fsntoaC[index*3] & fsntoaC[index*3+1] = fsntoaC[index*3]
cloud = reform(cloud(0,0,*,*))  
     cloud [*,index*3+2] = cloud [*,index*3] & cloud [*,index*3+1] = cloud [*,index*3]
CONCLD = reform(concld(0,0,*,*))  
     CONCLD [*,index*3+2] = CONCLD [*,index*3] & CONCLD [*,index*3+1] = CONCLD [*,index*3]

TEMP = reform(T(0,0,*,*)) 
TEMP = TEMP - TEMP(*,0) # (1+fltarr(n_elements(time)))

;;;; SCAM time is in hours I think: fix it
time = time - 122. ;; CAM dates start from September

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;; variables for plot: NCAR names
;!x.style = 1

tplot = time(time1:time2)
ntimes = n_elements(tplot)
pplot = lev
phalf = halflevels
NLEVS = n_elements(pplot)
plevs = [11,15,18,20,22,24] ;; nice pressure levels for correlations

LWTOA = flnt(time1:time2)
LWTOAC= flntC(time1:time2)
SWTOA = -fsntoa(time1:time2)
SWTOAC= -fsntoaC(time1:time2)

SWUPTOA = -reform( fsU(0,0,0,time1:time2) )
SWUPTOAC= -reform( fsUC(0,0,0,time1:time2)) 

SWCF = SWTOAC-SWTOA ;;; a negative number, reflected solar
LWCF = LWTOAC-LWTOA ;;; a positive number, reduced OLR
momeanCF = mean(SWCF+LWCF)

LWCUR = reform( FLU(0,0,*,time1:time2) - FLUC(0,0,*,time1:time2) )
LWCURD = LWCUR*0 ;; dummy arrray
for it = 0, time2-time1-1 do LWCURD(*,it) = deriv(LWCUR(*,it))/deriv(phalf)

SWCUR = reform( FSU(0,0,*,time1:time2) - FSUC(0,0,*,time1:time2) )
SWCURD = SWCUR*0 ;; dummy arrray
for it = 0, time2-time1-1 do SWCURD(*,it) = -deriv(SWCUR(*,it))/deriv(phalf)

QSW = reform( QRS(0,0,*,time1:time2) ) *86400. ;;; K/d
QLW = reform( QRL(0,0,*,time1:time2) ) *86400. ;;; K/d

Q1plot = reform( (DTV+DTCOND)(0,0,*,time1:time2) ) *86400. ;;; K/d
Q2plot = reform( -(VD01 +DCQ)(0,0,*,time1:time2) )    *86400. *2.5e6/1007. ;;; kg/kg/s -> K/d

;;;;; Cloud

MAXCWC = 0.5 ;; g/kg
CWC = reform( (CLDICE + CLDLIQ) (0,0,*,time1:time2) ) *1000 ;;; g/kg
CWClevels = [.0033,0.01,0.03,0.09,0.27,0.81] *MAXCWC ;;; g/kg

RHplot = reform( relhum(0,0,*,time1:time2) ) 
RHlevels = 10*findgen(15)
;;; recall this is a pert already TEMP = TEMP - TEMP(*,0) # (1+fltarr(n_elements(time)))
TPERT = Temp(*,time1:time2)

OMEG = reform(omega(0,0,*,time1:time2)) * 86400./100. ;; Pa/s -> hPa/d
OM = 300 ;;  
print, 'Omega scale hPa/d', max(abs(OMEG))
OMLEVELS = -OM + OM*findgen(21)/10. 

;;; Divergence from omega
DIV = OMEG
for it=0,n_elements(tplot)-1 do $
  div(*,it) = -deriv( OMEG(*,it) )/deriv(pplot) *1.e6/86400. ;; units 1e-6 s-1

;;;; Winds
Uplot = reform( U(0,0,*,time1:time2) )
Vplot = reform( V(0,0,*,time1:time2) )

;;;;; Rain
rain = reform( prect(0,0,time1:time2) ) *3600. *1000. ;;; mm/h rather than m/s!
rainLS = reform( precL(0,0,time1:time2) )  *3600. *1000. ;;; mm/h rather than m/s!
momeanrain = mean(rain)

cloudfraction = CLOUD(*,time1:time2) 
ccloudfraction = CONCLD(*,time1:time2) 

MODEL_NAME = 'CAM'
SMOO = 73
LEV500 = 18 ;; level where p~500, for Bony discrimination

;;;;; these need adding for pageplots
VERSION = 'rio'
site = station_string
CWClev_label = 'log contours'

;;;;;;;;;;; Filename for postscript AND HISTOGRAM SAVE FILES:
;;;;;;;;;;; model.'omeg'+omeg500.longitude.month.ps
monames = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
mostr = monames(imo)

fname = '../PSFILES/PS_'+MODEL_NAME +'/'+MODEL_NAME + '.' + $
  station_string + '.' + monames(imo) + '.' + $
  'rain'+str(fix(momeanrain*24)) + '.CRF' + str(fix(momeanCF)) + $
  '.ps'

tall
device, file=fname

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; MODEL
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; INDEPENDENT
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;PLOTTING CODE EXTERNAL/SHARED

@plots_page_gfdl.idl
@plots_regress.idl

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

psout
endfor ;;; IMONTH LOOP

end

