;;; Run with 
;;; process_CAM,'../Datasets/CAM/eul64x128_amip-156E_2S.2000.nc', 'COARE', [0,1,10,11]
;;; process_CAM,'../Datasets/CAM/eul64x128_amip-85W_20S.2000.nc', 'WHOI_85W', [0,3,7,10]

pro process_CAM_stratus, filename, station_string, months  ;;; prange 1000 to 600

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
fsu = reform(fsu(0,0,*,*))  
     fsu [*,index*3+2] = fsu [*,index*3] & fsu [*,index*3+1] = fsu [*,index*3]
fsuC= reform(fsuc(0,0,*,*)) 
     fsuC[*,index*3+2] = fsuC[*,index*3] & fsuC[*,index*3+1] = fsuC[*,index*3]
flu = reform(flu(0,0,*,*))  
     flu [*,index*3+2] = flu [*,index*3] & flu [*,index*3+1] = flu [*,index*3]
fluC= reform(fluc(0,0,*,*)) 
     fluC[*,index*3+2] = fluC[*,index*3] & fluC[*,index*3+1] = fluC[*,index*3]

qrs= reform(qrs(0,0,*,*)) 
     qrs[*,index*3+2] = qrs[*,index*3] & qrs[*,index*3+1] = qrs[*,index*3]
qrl= reform(qrl(0,0,*,*)) 
     qrl[*,index*3+2] = qrl[*,index*3] & qrl[*,index*3+1] = qrl[*,index*3]

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
phalf = halflevels ;; mb
NLEVS = n_elements(pplot)
plevs = [11,15,18,20,22,24] ;; nice pressure levels for correlations

LWTOA = flnt(time1:time2)
LWTOAC= flntC(time1:time2)
SWTOA = -fsntoa(time1:time2)
SWTOAC= -fsntoaC(time1:time2)

SWUPTOA = -fsU(0,time1:time2) 
SWUPTOAC= -fsUC(0,time1:time2)

SWCF = SWTOAC-SWTOA ;;; a negative number, reflected solar
LWCF = LWTOAC-LWTOA ;;; a positive number, reduced OLR
momeanCF = mean(SWCF+LWCF)

SWDC= fsdsc(time1:time2) ;;; down at surface, clear sky ~ incoming, positive values
NSWCF = SWCF / (SWDC > 1) ;; normalized SWCF, will be 0 at night, negative values like SWCF

LWCUR = FLU(*,time1:time2) - FLUC(*,time1:time2) 
LWCURD = LWCUR*0 ;; dummy arrray
for it = 0, time2-time1-1 do LWCURD(*,it) = -deriv(LWCUR(*,it))/deriv(phalf) ;;; Wm-2 per mb

SWCUR =  FSU(*,time1:time2) - FSUC(*,time1:time2) 
SWCURD = SWCUR*0 ;; dummy arrray
for it = 0, time2-time1-1 do SWCURD(*,it) = deriv(SWCUR(*,it))/deriv(phalf)

;;;; Normalize by incoming clear-sky radiation
NSWCURD = SWCUR*0 ;; dummy arrray
for ip = 0, n_elements(phalf)-1 do NSWCURD(ip,*) = SWCURD(ip,*)/(SWDC >1) ;; avoid /0

QSW = QRS(*,time1:time2)  *86400. ;;; K/d
QLW = QRL(*,time1:time2)  *86400. ;;; K/d

Q1plot = QSW+QLW+ reform( (DTCOND+DTV)(0,0,*,time1:time2) ) *86400. ;;; K/d
Q2plot = reform( -(VD01 +DCQ)(0,0,*,time1:time2) )    *86400. *2.5e6/1007. ;;; kg/kg/s -> K/d

;;; include DTV (diffusion) if available
Q1plot = QSW+QLW+ reform( (DTV+DTCOND)(0,0,*,time1:time2) )*86400. 

;;;;; Cloud

MAXCWC = 0.5 ;; g/kg
CWC = reform( (CLDICE + CLDLIQ) (0,0,*,time1:time2) ) *1000 ;;; g/kg
CWClevels = [.0033,0.01,0.03,0.09,0.27,0.81] *MAXCWC ;;; g/kg

;;;; after KW256 but with a low one too
CWClevels = [0.002, 0.01*( 1+indgen(10) )*2]
CWClev_label = '2,20,40,60... mg/kg'

RHplot = reform( relhum(0,0,*,time1:time2) ) 
RHlevels = 10*findgen(11)
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
rainLS2 = reform( precSH(0,0,time1:time2) )  *3600. *1000. ;;; Shallow (Hack) rain
momeanrain = mean(rain)

cloudfraction = CLOUD(*,time1:time2) 
ccloudfraction = CONCLD(*,time1:time2) 

MODEL_NAME = 'CAM'
SMOO = 73
LEV500 = 18 ;; level where p~500, for Bony discrimination

;;;;; these need adding for pageplots
VERSION = 'rio'
site = station_string

;;;;;;;;;;; Filename for postscript AND HISTOGRAM SAVE FILES:
;;;;;;;;;;; model.'omeg'+omeg500.longitude.month.ps
monames = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
mostr = monames(imo)

PSdirstub = '../PSFILES/PS_'+MODEL_NAME +'/'
namestub = MODEL_NAME + '.' + station_string + '.' + mostr
statstub = '.rain'+str(fix(momeanrain*24)) + 'mmd.CRF' + str(fix(momeanCF)) +'W.OM'+ str(fix(mean(omeg(LEV500,*))))+'hpad'

fname = PSdirstub+namestub+statstub+'.stratus.ps'
print, fname

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; MODEL
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; INDEPENDENT
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;PLOTTING CODE EXTERNAL/SHARED

;;;; Now main processing
tall
device, file=fname

@plots_page_lowertrop.idl
;@plots_regress.idl
@plots_regress_TOAbase.idl

psout

;@rain_composites_only.idl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

endfor ;;; IMONTH LOOP
;stop

end

