;PRO

; track the system time
print, systime(/utc)

; total wavenumber, generate Spherical Harmonics Functions from n=NN1 to n=NN2
NN1 = 1
NN2 = 360 
totalN = NN2-NN1+1

; read the data and coordinates
pfname = '../input/RestartLatLon_ne60.cam.h1.0001-01-02-00000.nc'
fid = ncdf_open(pfname)
ncdf_varget, fid, 'lat', lat
ncdf_varget, fid, 'lon', lon
ncdf_varget, fid, 'area', area
ncdf_close, fid
ncol = n_elements(lon)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; geographical parameters

; change lat and lon to colatitude and longitude in radian
r_lon = lon*!dpi/180.d
colat = lat + 90.d
r_colat = colat*!dpi/180.d
; x is cosine of colatitude
x = cos(r_colat)
; y is sine of colatitude
y = sin(r_colat)

tol = 100.d*(machar(/double)).xmin

; specifying the parameters at certain angles increase precision

;   DO NOT replace any numbers that are not at 0 degree, 90 degrees, 180
; degrees, and 270 degrees. It messes up the signals in the spectrum tail. 

cos0 = 1.d
sin0 = 0.d
cos30 = sqrt(3.d)/2.d
sin30 = 0.5d
cos45 = 1.d/sqrt(2.d)
sin45 = 1.d/sqrt(2.d)
cos60 = 0.5d
sin60 = sqrt(3.d)/2.d
cos90 = 0.d
sin90 = 1.d

i_cos = where(abs(r_colat) lt tol)
x[i_cos] = cos0
;i_cos = where(abs(r_colat-!dpi/6.d) lt tol)
;x[i_cos] = cos30
;i_cos = where(abs(r_colat-!dpi/4.d) lt tol)
;x[i_cos] = cos45
;i_cos = where(abs(r_colat-!dpi/3.d) lt tol)
;x[i_cos] = cos60
i_cos = where(abs(r_colat-!dpi/2.d) lt tol)
x[i_cos] = cos90
;i_cos = where(abs(r_colat-!dpi*2.d/3.d) lt tol)
;x[i_cos] = -cos60
;i_cos = where(abs(r_colat-!dpi*3.d/4.d) lt tol)
;x[i_cos] = -cos45
;i_cos = where(abs(r_colat-!dpi*5.d/6.d) lt tol)
;x[i_cos] = -cos30
i_cos = where(abs(r_colat-!dpi) lt tol)
x[i_cos] = -cos0

i_xeq1 = where(abs(x-1.d) lt tol or abs(x+1.d) lt tol)

j_sin = where(abs(r_colat) lt tol)
y[j_sin] = sin0
;j_sin = where(abs(r_colat-!dpi/6.d) lt tol)
;y[j_sin] = sin30
;j_sin = where(abs(r_colat-!dpi/4.d) lt tol)
;y[j_sin] = sin45
;j_sin = where(abs(r_colat-!dpi/3.d) lt tol)
;y[j_sin] = sin60
j_sin = where(abs(r_colat-!dpi/2.d) lt tol)
y[j_sin] = sin90
;j_sin = where(abs(r_colat-!dpi*2.d/3.d) lt tol)
;y[j_sin] = sin60
;j_sin = where(abs(r_colat-!dpi*3.d/4.d) lt tol)
;y[j_sin] = sin45
;j_sin = where(abs(r_colat-!dpi*5.d/6.d) lt tol)
;y[j_sin] = sin30
j_sin = where(abs(r_colat-!dpi) lt tol)
y[j_sin] = sin0

; geographical parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; recurrence parameters
; exp(i * lambda)

; specifying the parameters at certain angles increase precision
eil = exp(dcomplex(0.d,r_lon))

;   DO NOT replace any numbers that are not at 0 degree, 90 degrees, 180
; degrees, and 270 degrees. It messes up the signals in the spectrum tail. 

k = where(abs(r_lon) lt tol)
eil[k] = dcomplex(cos0,sin0)
;k = where(abs(r_lon-!dpi/6.d) lt tol)
;eil[k] = dcomplex(cos30,sin30)
;k = where(abs(r_lon-!dpi/4.d) lt tol)
;eil[k] = dcomplex(cos45,sin45)
;k = where(abs(r_lon-!dpi/3.d) lt tol)
;eil[k] = dcomplex(cos60,sin60)
k = where(abs(r_lon-!dpi/2.d) lt tol)
eil[k] = dcomplex(cos90,sin90)
;k = where(abs(r_lon-!dpi*2.d/3.d) lt tol)
;eil[k] = dcomplex(-cos60,sin60)
;k = where(abs(r_lon-!dpi*3.d/4.d) lt tol)
;eil[k] = dcomplex(-cos45,sin45)
;k = where(abs(r_lon-!dpi*5.d/6.d) lt tol)
;eil[k] = dcomplex(-cos30,sin30)
k = where(abs(r_lon-!dpi) lt tol)
eil[k] = dcomplex(-cos0,sin0)
;k = where(abs(r_lon-!dpi*7.d/6.d) lt tol)
;eil[k] = dcomplex(-cos30,-sin30)
;k = where(abs(r_lon-!dpi*5.d/4.d) lt tol)
;eil[k] = dcomplex(-cos45,-sin45)
;k = where(abs(r_lon-!dpi*4.d/3.d) lt tol)
;eil[k] = dcomplex(-cos60,-sin60)
k = where(abs(r_lon-!dpi*3.d/2.d) lt tol)
eil[k] = dcomplex(cos90,-sin90)
;k = where(abs(r_lon-!dpi*5.d/3.d) lt tol)
;eil[k] = dcomplex(cos60,-sin60)
;k = where(abs(r_lon-!dpi*7.d/4.d) lt tol)
;eil[k] = dcomplex(cos45,-sin45)
;k = where(abs(r_lon-!dpi*11.d/6.d) lt tol)
;eil[k] = dcomplex(cos30,-sin30)
k = where(abs(r_lon-!dpi*2.d) lt tol)
eil[k] = dcomplex(cos0,sin0)

Factor1 = y * eil
Factor2 = x * eil
; recurrence parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Start New from n=0

if(NN1 le 3) then begin

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; for n=0 calculate Spherical Harmonics Functions, only for Y

  n = 0
  m = 0

  ;;;;;;;;;;;;;;;;;;;;;;
  ; Y
  print, 'n, m = ', n, m
  ; Note that rlon is a ncol-element array, so is the return value Ynm
  Ynm = spher_harm(r_colat,r_lon,n,m,/double)
  ; Y
  ;;;;;;;;;;;;;;;;;;;;;

  pfname_out = '../output/ne60/Y_n=0000.nc'
  spawn, 'rm '+pfname_out+' '+pfname_out+'.gz'
  fid = ncdf_create(pfname_out,/clobber,/netcdf4)
  gpid = ncdf_dimdef(fid,"GridPoint#",ncol)
  mid = ncdf_dimdef(fid,"m",n+1)
  RIid = ncdf_dimdef(fid,"ReIm",2)
  lonid = ncdf_vardef(fid,"lon",[gpid],/double)
  latid = ncdf_vardef(fid,"lat",[gpid],/double)
  Yid = ncdf_vardef(fid,"Y",[gpid,mid,RIid],/double)
  ncdf_varput, fid, "lon", lon
  ncdf_varput, fid, "lat", lat
  temp_array = dblarr(ncol,n+1,2)
  temp_array[*,*,0] = real_part(Ynm[*,*])
  temp_array[*,*,1] = Imaginary(Ynm[*,*])
  ncdf_varput, fid, "Y", temp_array
  ncdf_close, fid
  ; use temporary function to delete temp_array
  tempvar = size(temporary(temp_array))
;  spawn, 'gzip '+pfname_out
  ; for n=0 calculate Spherical Harmonics Functions, only for Y
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; for n=1 calculate Spherical Harmonics Functions

  n = 1
  Y_nM2_m = dcomplexarr(ncol,n+1)
  vPSI_nM2_m = dcomplexarr(ncol,n+1,2)

  ;;;;;;;;;;;;;;;;;;;;;;
  ; Y
  for m = 0, n do begin
    print, 'n, m = ', n, m
    ; Note that rlon is a ncol-element array, so is the return value Ynm
    Ynm = spher_harm(r_colat,r_lon,n,m,/double)
    Y_nM2_m[*,m] = Ynm
  endfor
  ; Y
  ;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;
  ; PSI
  m = 0
  vPSInm_theta = - sqrt(3.d/4.d/!dpi) * dcomplex(y,0.d)
  vPSInm_lambda = dcomplexarr(ncol)
  vPSI_nM2_m[*,m,0] = vPSInm_theta
  vPSI_nM2_m[*,m,1] = vPSInm_lambda
  ;
  m = 1
  vPSInm_theta = - sqrt(3.d/8.d/!dpi) * eil * x
  vPSInm_lambda = - sqrt(3.d/8.d/!dpi) * eil * dcomplex(0.d,1.d)
  vPSI_nM2_m[*,m,0] = vPSInm_theta
  vPSI_nM2_m[*,m,1] = vPSInm_lambda
  ; PSI
  ;;;;;;;;;;;;;;;;;;;;;

  pfname_out = '../output/ne60/Y_n=0001.nc'
  spawn, 'rm '+pfname_out+' '+pfname_out+'.gz'
  fid = ncdf_create(pfname_out,/clobber,/netcdf4)
  gpid = ncdf_dimdef(fid,"GridPoint#",ncol)
  mid = ncdf_dimdef(fid,"m",n+1)
  RIid = ncdf_dimdef(fid,"ReIm",2)
  Yid = ncdf_vardef(fid,"Y",[gpid,mid,RIid],/double)
  temp_array = dblarr(ncol,n+1,2)
  temp_array[*,*,0] = real_part(Y_nM2_m[*,*])
  temp_array[*,*,1] = Imaginary(Y_nM2_m[*,*])
  ncdf_varput, fid, "Y", temp_array
  ncdf_close, fid
;  spawn, 'gzip '+pfname_out

  pfname_out = '../output/ne60/PSI_n=0001.nc'
  spawn, 'rm '+pfname_out+' '+pfname_out+'.gz'
  fid = ncdf_create(pfname_out,/clobber,/netcdf4)
  gpid = ncdf_dimdef(fid,"GridPoint#",ncol)
  mid = ncdf_dimdef(fid,"m",n+1)
  RIid = ncdf_dimdef(fid,"ReIm",2)
  PSI_thetaid = ncdf_vardef(fid,"PSI_theta",[gpid,mid,RIid],/double)
  PSI_lambdaid = ncdf_vardef(fid,"PSI_lambda",[gpid,mid,RIid],/double)
  temp_array[*,*,0] = real_part(vPSI_nM2_m[*,*,0])
  temp_array[*,*,1] = Imaginary(vPSI_nM2_m[*,*,0])
  ncdf_varput, fid, "PSI_theta", temp_array
  temp_array[*,*,0] = real_part(vPSI_nM2_m[*,*,1])
  temp_array[*,*,1] = Imaginary(vPSI_nM2_m[*,*,1])
  ncdf_varput, fid, "PSI_lambda", temp_array
  ncdf_close, fid
  ; use temporary function to delete temp_array
  tempvar = size(temporary(temp_array))
  ; use temporary function to delete vPSI_nM2_m
  tempvar = size(temporary(vPSI_nM2_m))
;  spawn, 'gzip '+pfname_out

  ; for n=1 calculate Spherical Harmonics Functions
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; for n=2 calculate Spherical Harmonics Functions

  n = 2
  Y_nM1_m = dcomplexarr(ncol,n+1)
  vPSI_nM1_m = dcomplexarr(ncol,n+1,2)

  ;;;;;;;;;;;
  ; Y
  for m = 0, n do begin
    print, 'n, m = ', n, m
    ; Note that rlon is a ncol-element array, so is the return value Ynm
    Ynm = spher_harm(r_colat,r_lon,n,m,/double)
    Y_nM1_m[*,m] = Ynm
  endfor
  ; Y
  ;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;
  ; PSI
  m = 0
  vPSInm_theta = - 3.d/2.d * sqrt(5.d/!dpi) * dcomplex(x*y,0.d)
  vPSInm_lambda = dcomplexarr(ncol)
  vPSI_nM1_m[*,m,0] = vPSInm_theta
  vPSI_nM1_m[*,m,1] = vPSInm_lambda
  ;
  m = 1
  vPSInm_theta = - sqrt(15.d/8.d/!dpi) * eil * cos(2.d*r_colat)
  vPSInm_lambda = - sqrt(15.d/8.d/!dpi) * eil * dcomplex(0.d,x)
  vPSI_nM1_m[*,m,0] = vPSInm_theta
  vPSI_nM1_m[*,m,1] = vPSInm_lambda
  ;
  m = 2
  vPSInm_theta = sqrt(15.d/8.d/!dpi) * y * eil*eil * x
  vPSInm_lambda = sqrt(15.d/8.d/!dpi) * y * eil*eil * dcomplex(0.d,1.d)
  vPSI_nM1_m[*,m,0] = vPSInm_theta
  vPSI_nM1_m[*,m,1] = vPSInm_lambda
  ; PSI
  ;;;;;;;;;;;;;;;;;;;

  pfname_out = '../output/ne60/Y_n=0002.nc'
  spawn, 'rm '+pfname_out+' '+pfname_out+'.gz'
  fid = ncdf_create(pfname_out,/clobber,/netcdf4)
  gpid = ncdf_dimdef(fid,"GridPoint#",ncol)
  mid = ncdf_dimdef(fid,"m",n+1)
  RIid = ncdf_dimdef(fid,"ReIm",2)
  Yid = ncdf_vardef(fid,"Y",[gpid,mid,RIid],/double)
  temp_array = dblarr(ncol,n+1,2)
  temp_array[*,*,0] = real_part(Y_nM1_m[*,*])
  temp_array[*,*,1] = Imaginary(Y_nM1_m[*,*])
  ncdf_varput, fid, "Y", temp_array
  ncdf_close, fid
;  spawn, 'gzip '+pfname_out

  pfname_out = '../output/ne60/PSI_n=0002.nc'
  spawn, 'rm '+pfname_out+' '+pfname_out+'.gz'
  fid = ncdf_create(pfname_out,/clobber,/netcdf4)
  gpid = ncdf_dimdef(fid,"GridPoint#",ncol)
  mid = ncdf_dimdef(fid,"m",n+1)
  RIid = ncdf_dimdef(fid,"ReIm",2)
  PSI_thetaid = ncdf_vardef(fid,"PSI_theta",[gpid,mid,RIid],/double)
  PSI_lambdaid = ncdf_vardef(fid,"PSI_lambda",[gpid,mid,RIid],/double)
  temp_array[*,*,0] = real_part(vPSI_nM1_m[*,*,0])
  temp_array[*,*,1] = Imaginary(vPSI_nM1_m[*,*,0])
  ncdf_varput, fid, "PSI_theta", temp_array
  temp_array[*,*,0] = real_part(vPSI_nM1_m[*,*,1])
  temp_array[*,*,1] = Imaginary(vPSI_nM1_m[*,*,1])
  ncdf_varput, fid, "PSI_lambda", temp_array
  ncdf_close, fid
  ; use temporary function to delete temp_array
  tempvar = size(temporary(temp_array))
;  spawn, 'gzip '+pfname_out

  ; for n=2 calculate Spherical Harmonics Functions
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

endif else begin

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; for n=NN1-2 read in Spherical Harmonics Functions

  n = NN1-2
  Y_nM2_m = dcomplexarr(ncol,n+1)
  vPSI_nM2_m = dcomplexarr(ncol,n+1,2)
  temp_array = dblarr(ncol,n+1,2)

  fm = '(i4.4)'


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Y
  pfname_in = '../output/ne60/Y_n='+string(n,form=fm)+'.nc'
  spawn, 'gunzip '+pfname_in+'.gz'
  fid_in = ncdf_open(pfname_in)
  ncdf_varget, fid_in, 'Y', temp_array
  for icol = 0, ncol-1 do begin
  for m = 0, n do begin
    Y_nM2_m[icol,m] = dcomplex(temp_array[icol,m,0],temp_array[icol,m,1])
  endfor
  endfor
  ncdf_close, fid_in
;  spawn, 'gzip '+pfname_in
  ; Y
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

  ;;;;;;;;;;;;;;;;;;;;;
  ; PSI
  pfname_in = '../output/ne60/PSI_n='+string(n,form=fm)+'.nc'
  spawn, 'gunzip '+pfname_in+'.gz'
  fid_in = ncdf_open(pfname_in)
  ncdf_varget, fid_in, 'PSI_theta', temp_array
  for icol = 0, ncol-1 do begin
  for m = 0, n do begin
    vPSI_nM2_m[icol,m,0] = dcomplex(temp_array[icol,m,0],temp_array[icol,m,1])
  endfor
  endfor
  ncdf_varget, fid_in, 'PSI_lambda', temp_array
  for icol = 0, ncol-1 do begin
  for m = 0, n do begin
    vPSI_nM2_m[icol,m,1] = dcomplex(temp_array[icol,m,0],temp_array[icol,m,1])
  endfor
  endfor
  ncdf_close, fid_in
;  spawn, 'gzip '+pfname_in
  ; PSI
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; use temporary function to delete temp_array
  tempvar = size(temporary(temp_array))

  ; for n=NN1-2 read in Spherical Harmonics Functions
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; for n=NN1-1 read in Spherical Harmonics Functions

  n = NN1-1
  Y_nM1_m = dcomplexarr(ncol,n+1)
  vPSI_nM1_m = dcomplexarr(ncol,n+1,2)
  temp_array = dblarr(ncol,n+1,2)

  fm = '(i4.4)'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Y
  pfname_in = '../output/ne60/Y_n='+string(n,form=fm)+'.nc'
  spawn, 'gunzip '+pfname_in+'.gz'
  fid_in = ncdf_open(pfname_in)
  ncdf_varget, fid_in, 'Y', temp_array
  for icol = 0, ncol-1 do begin
  for m = 0, n do begin
    Y_nM1_m[icol,m] = dcomplex(temp_array[icol,m,0],temp_array[icol,m,1])
  endfor
  endfor
  ncdf_close, fid_in
;  spawn, 'gzip '+pfname_in
  ; Y
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

  ;;;;;;;;;;;;;;;;;;;;;
  ; PSI
  pfname_in = '../output/ne60/PSI_n='+string(n,form=fm)+'.nc'
  spawn, 'gunzip '+pfname_in+'.gz'
  fid_in = ncdf_open(pfname_in)
  ncdf_varget, fid_in, 'PSI_theta', temp_array
  for icol = 0, ncol-1 do begin
  for m = 0, n do begin
    vPSI_nM1_m[icol,m,0] = dcomplex(temp_array[icol,m,0],temp_array[icol,m,1])
  endfor
  endfor
  ncdf_varget, fid_in, 'PSI_lambda', temp_array
  for icol = 0, ncol-1 do begin
  for m = 0, n do begin
    vPSI_nM1_m[icol,m,1] = dcomplex(temp_array[icol,m,0],temp_array[icol,m,1])
  endfor
  endfor
  ncdf_close, fid_in
;  spawn, 'gzip '+pfname_in
  ; PSI
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; use temporary function to delete temp_array
  tempvar = size(temporary(temp_array))

  ; for n=NN1-1 read in Spherical Harmonics Functions
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




endelse


;;;;;;;;;;;;;;;
; the iteration for computing Spherical Harmonics Functions from n=NN1 to n=NN2
for n = max([NN1,3]), NN2 do begin

  Y_n_m = dcomplexarr(ncol,n+1)
  vPSI_n_m = dcomplexarr(ncol,n+1,2)
  SpherHarm_recur2, $
    r_lon, x, y, eil, Factor1, Factor2, area, n, i_xeq1, $
    Y_nM2_m, Y_nM1_m, $
    reform(vPSI_nM1_m[*,n-1,*]), $
    Y_n_m, vPSI_n_m

  Y_nM2_m = Y_nM1_m
  Y_nM1_m = Y_n_m
  vPSI_nM1_m = vPSI_n_m

  fm = '(i4.4)'

  pfname_out = '../output/ne60/Y_n='+string(n,form=fm)+'.nc'
  spawn, 'rm '+pfname_out+' '+pfname_out+'.gz'
  fid = ncdf_create(pfname_out,/clobber,/netcdf4)
  gpid = ncdf_dimdef(fid,"GridPoint#",ncol)
  mid = ncdf_dimdef(fid,"m",n+1)
  RIid = ncdf_dimdef(fid,"ReIm",2)
  Yid = ncdf_vardef(fid,"Y",[gpid,mid,RIid],/double)
  temp_array = dblarr(ncol,n+1,2)
  temp_array[*,*,0] = real_part(Y_n_m[*,*])
  temp_array[*,*,1] = Imaginary(Y_n_m[*,*])
  ncdf_varput, fid, "Y", temp_array
  ncdf_close, fid
;  spawn, 'gzip '+pfname_out

  pfname_out = '../output/ne60/PSI_n='+string(n,form=fm)+'.nc'
  spawn, 'rm '+pfname_out+' '+pfname_out+'.gz'
  fid = ncdf_create(pfname_out,/clobber,/netcdf4)
  gpid = ncdf_dimdef(fid,"GridPoint#",ncol)
  mid = ncdf_dimdef(fid,"m",n+1)
  RIid = ncdf_dimdef(fid,"ReIm",2)
  PSI_thetaid = ncdf_vardef(fid,"PSI_theta",[gpid,mid,RIid],/double)
  PSI_lambdaid = ncdf_vardef(fid,"PSI_lambda",[gpid,mid,RIid],/double)
  temp_array[*,*,0] = real_part(vPSI_nM1_m[*,*,0])
  temp_array[*,*,1] = Imaginary(vPSI_nM1_m[*,*,0])
  ncdf_varput, fid, "PSI_theta", temp_array
  temp_array[*,*,0] = real_part(vPSI_nM1_m[*,*,1])
  temp_array[*,*,1] = Imaginary(vPSI_nM1_m[*,*,1])
  ncdf_varput, fid, "PSI_lambda", temp_array
  ncdf_close, fid
  ; use temporary function to delete temp_array
  tempvar = size(temporary(temp_array))
;  spawn, 'gzip '+pfname_out

endfor


END
