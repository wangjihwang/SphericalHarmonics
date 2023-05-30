PRO SpherHarm_recur2, $
  r_lon, x, y, eil, Factor1, Factor2, area, n, i_xeq1, $
  Y_nM2_m, Y_nM1_m, $
  vPSInm_last, $
  Y_n_m, vPSI_n_m
;;;;;;;;;;;;;;;;;;;;;;;;;;
; Author: Jih-Wang (Aaron) Wang
; Date: 2016/7/14
;;;;;;;;;;;;;;;;;;;;;;;;;;
; Purpose:
;   This program takes scalar and vector spherical harmonics functions for
; last two total wavenumbers (n-2 and n-1) and use recurrence relations to
; calculate the scalar and vector spherical harmonics functions for next total
; wavenumber (n).
;;;;;;;;;;;;;;;;;;;;;;;;;;
; Note:
;   Only the orders greater than or equal to 0 (m = 0,1,2,...,n) are
; considered, because the functions for orders less than 0
; (m = -1,-2,...,-n) are simply either conjugates (even order) or negative
; conjugates (odd order) of the functions for orders greater than 0 (i.e.,
; f_n_-m = (-1)^m * conj(f_n_m)). The algorithm in use is chosen so that the
; errors resulting from recurrence do not accumulate too quickly.
;;;;;;;;;;;;;;;;;;;;;;;;;;
; INPUT
; ncol: data point counts
; lon (1-D float/double array; [ncol]): longitudes of the data points on the
;   sphere
; lat (1-D float/double array; [ncol]): latitudes of the data points on the
;   sphere
; area (1-D double array; [ncol]): area of each grid cell
;   total of area array should fully cover the sphere (4.*!pi)
; n (integer): degree of associated legendre polynomials or total wavenumber of
;   output
; nM2: total wavenumber = n-2
; nM1: total wavenumber = n-1
; m: order of associated legendre polynomials (in certain orienation equal to
;   zonal wavenumber)
; Y_nM2_m (2-D complex float/double array; [ncol,n-1]): scalar spherical
;   harmonics functions for total wavenumber=n-2
; Y_nM1_m (2-D complex float/double array; [ncol,n]): scalar spherical
;   harmonics functions for total wavenumber=n-1
; vPSInm_last (2-D complex float/double array; [ncol,2]): last vector
;   spherical harmonics functions for m=n-1
;;;;;;;;;;;;;;;;;;;;;;
; OUTPUT
; Ynm (2-D complex double array; [ncol,n+1]): scalar spherical harmonics
;   functions for total wavenumber=n
; vPSInm (3-D complex double array; [ncol,n+1,2]); first vector spherical
;   harmonics functions for total wavenumber=n
;;;;;;;;;;;;;;;;;;;;;;

ncol = n_elements(area)
Ynm = dcomplexarr(ncol,n+1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; m = 0
m = 0
print, 'n, m = ', n, m

a = - sqrt((n-1.d)*(n-1.d)/(2.d*n-3.d)/(2.d*n-1.d))
b = sqrt((2.d*n-1.d)*(2.d*n+1.d)/double(n)/double(n))
xb = x*b
ab = a*b
Ynm[*,0] = xb*Y_nM1_m[*,0]+ab*Y_nM2_m[*,0]
; special case for x=1 or -1, and m=0
Ynm[i_xeq1,0] = sqrt((2.d*n+1.d)/4.d/!dpi)
Y_n_m[*,0] = Ynm[*,0]
; m = 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; m = 1, 2, ..., n-2
for m = 1, n-2 do begin
  print, 'n, m = ', n, m

  a = - sqrt((n-m-1.d)*(n+m-1.d)/(2.d*n-3.d)/(2.d*n-1.d))
  b = sqrt((2.d*n-1.d)*(2.d*n+1.d)/double(n-m)/double(n+m))
  xb = x*b
  ab = a*b

  Ynm[*,m] = xb*Y_nM1_m[*,m]+ab*Y_nM2_m[*,m]

  ; special case when x=1 or -1 and m/=0
  Ynm[i_xeq1,m] = dcomplex(0.d,0.d)

  ; DO NOT replace any numbers that are not at poles. It will mess up
  ; the signals in spectrum tail.

  Y_n_m[*,m] = Ynm[*,m]
endfor
; m = 1, 2, ..., n-2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; m = n-1
m = n-1
print, 'n, m = ', n, m

g = sqrt(2.d*n+1.d)
gx = g*x
Ynm[*,m] = gx * Y_nM1_m[*,m]

; special case when x=1 or -1 and m/=0
Ynm[i_xeq1,m] = dcomplex(0.d,0.d)

; DO NOT replace any numbers that are not at poles. It will mess up
; the signals in spectrum tail.

Y_n_m[*,m] = Ynm[*,m]
; m = n-1
;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;
; m = n
m = n
print, 'n, m = ', n, m

g = - sqrt((2.d*n+1.d)/(2.d*n))
gFactor1 = g*Factor1
Ynm[*,m] = gFactor1*Y_nM1_m[*,m-1]

; special case when x=1 or -1 and m/=0
Ynm[i_xeq1,m] = dcomplex(0.d,0.d)

; DO NOT replace any numbers that are not at poles. It will mess up
; the signals in spectrum tail.

Y_n_m[*,m] = Ynm[*,m]
; m = n
;;;;;;;;;;;;;;;;;;;;;;;;;

; Y
;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; vPSI

;;;;;;;;;;;;;;;;;;;;;;;;;;
; m = 0, 1, ..., n-1

for m = 0, n-1 do begin

  mx_y = m * x / y
  m_y = m / y

  Ynm = Y_n_m[*,m]
  Ynmp1 = Y_n_m[*,m+1]

  vPSInm_theta = mx_y*Ynm + sqrt(double(n-m)*double(n+m+1))*Ynmp1/eil
  vPSInm_lambda = dcomplex(0.d,1.d)*m_y*Ynm
  vPSInm_theta[i_xeq1] = dcomplex(0.d,0.d)
  vPSInm_lambda[i_xeq1] = dcomplex(0.d,0.d)
  vPSI_n_m[*,m,0] = vPSInm_theta
  vPSI_n_m[*,m,1] = vPSInm_lambda

endfor
; m = 0, 1, ..., n-1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;
; m = n
m = n
print, 'n, m = ', n, m

g = - sqrt((2.d*n+1.d)/(2.d*n))
gFactor1 = g*Factor1
gFactor2 = g*Factor2

vPSInm_theta = g * (Factor1*vPSInm_last[*,0]+Factor2*Y_nM1_m[*,m-1])
vPSInm_lambda = g * (Factor1*vPSInm_last[*,1]+eil*dcomplex(0.d,1.d)*Y_nM1_m[*,m-1])
vPSI_n_m[*,m,0] = vPSInm_theta
vPSI_n_m[*,m,1] = vPSInm_lambda
; m = n
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; vPSI
;;;;;;;;;;;;;;;;;;;;;;;;;;;;


END
